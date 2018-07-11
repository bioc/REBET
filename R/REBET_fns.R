#REBET function:
#Y: vector of Phenotye
#G: n*p matrix of Genotype
#E: sub-region annotation
#X: matrix of covariates
#a,b: parameters of beta distribution

rebet <- function(Y, G, E, outcome=NULL, X=NULL, a=1, b=1, saveMem=FALSE)
{
  #need add some sanity check of parameter of REBET
  if (!is.matrix(G)) stop("ERROR: G must be a matrix with at least two columns")
  n <- nrow(G)
  p <- ncol(G)
  if ((is.null(p)) || (p < 2)) stop("ERROR: G must be a matrix with at least two columns")

  E.unique <- unique(E)
  if (p != length(E)) stop("ERROR: length of E and number of col of G must be the same")
  if (a <= 0) stop("ERROR: a must be positive")
  if (b <= 0) stop("ERROR: b must be positive")
  #if ((p.bound < 0) || (p.bound > 1)) stop("ERROR: p.bound must be between 0 and 1")

  if (n != length(Y)) stop("ERROR: length of Y and number of rows of G must be the same")
  if (!is.null(outcome)) {
    if (!(outcome %in% c("continuous","binary"))) stop("ERROR: outcome type must be continous, binary or NULL")
  } else {
    temp <- is.finite(Y)
    if (all(Y[temp] %in% 0:1)) {
      outcome <- "binary"
    } else {
      outcome <- "continuous"
    }
  }
  if (is.null(X)) X <- matrix(data=1, nrow=n, ncol=1)
  if (is.vector(X)) X <- matrix(X, ncol=1)
  if ((!is.matrix(X)) && (!is.vector(X))) stop("ERROR: X must be a matrix or vector")
  if (n != nrow(X)) stop("ERROR: length of Y and number of rows of X must be the same")

  # Remove missing values
  temp <- (!is.finite(Y)) | (rowSums(!is.finite(G)) > 0) | (rowSums(!is.finite(X)) > 0)  
  if (any(temp)) {
    Y <- Y[!temp]
    G <- G[!temp, , drop=FALSE]
    X <- X[!temp, , drop=FALSE]
    n <- nrow(G)
  }
  if (!n) stop("ERROR: No data to process")

  # Add intercept if needed
  if (all(colSums(X == 1) != n)) X <- as.matrix(cbind(1, X))
  ncov <- ncol(X)

  # Local function to compute diagonal of SIGMA for continuous case
  getSIGMA2_cont <- function() {

    if (saveMem) {
      ret <- double(p)
      tmp <- .C("SIGMA2_cont", as.numeric(G), as.integer(n), as.integer(p), 
                as.numeric(X), as.integer(ncov), as.numeric(inv.XtX),
                as.numeric(tilde.sigma2), ret=ret, PACKAGE="REBET")
      ret <- tmp$ret 
    } else {
      SIGMA <- Gt%*% (diag(1, n) - X %*% inv.XtX %*% t(X)) %*% G * tilde.sigma2  
      ret   <- diag(SIGMA)
    }
    ret

  } # END: getSIGMA2_cont

  # Local function to compute diagonal of SIGMA for binary case
  getSIGMA2_binary <- function() {

    if (saveMem) {
      ret <- double(p)
      tmp <- .C("SIGMA2_binary", as.numeric(G), as.integer(n), as.integer(p), 
                as.numeric(X), as.integer(ncov), as.numeric(inv.XtDX),
                as.numeric(d), ret=ret, PACKAGE="REBET")
      ret <- tmp$ret 
    } else {
      SIGMA <- Gt %*% (diag(d) - (d %o% d) * (X %*% (inv.XtDX) %*% t(X))) %*% G
      ret   <- diag(SIGMA)
    }
    ret

  } # END: getSIGMA2_binary

  # Local function to compute XtDX <- t(X)%*%D%*%X for binary case
  getXtDX <- function() {

    if (saveMem) {
      ret <- double(ncov*ncov)
      tmp <- .C("computeXtDX", as.numeric(X), as.numeric(d), as.integer(n), 
                as.integer(ncov), ret=ret, PACKAGE="REBET")
      ret <- matrix(tmp$ret, nrow=ncov, ncol=ncov)
    } else {
      ret <- t(X)%*%diag(d)%*%X
    }
    ret
  
  } # END: getXtDX

  #step 1:calculate score statistics vector and its variance matrix
  if (outcome == "continuous"){
    
    fit          <- try(lm(Y ~ X - 1))
    if ("try-error" %in% class(fit)) stop("ERROR: Linear regression failed")
    tilde.sigma2 <- summary(fit)$sigma^2
    res          <- fit$resid   
    Gt           <- t(G)
    S            <- Gt %*% res 
    
    XtX          <- t(X)%*%X
    inv.XtX      <- try(solve(XtX))
    if ("try-error" %in% class(inv.XtX)) {
      print(summary(fit))
      stop("ERROR: matrix inversion failed")
    }
    
    SIGMA2       <- getSIGMA2_cont() 

    rm(Gt, XtX, res)
    gc()
  }else if (outcome == "binary"){
    fit      <- try(glm(formula = Y ~ X - 1, family = binomial(link = logit)))
    if (("try-error" %in% class(fit)) || (!fit$converged)) {
      stop("ERROR: Logistic regression failed")
    }
    mu       <- fit$fitted.value
    res      <- Y - mu
    Gt       <- t(G)
    S        <- Gt %*% res 
    d        <- mu * (1 - mu)
    XtDX     <- getXtDX()
    inv.XtDX <- try(solve(XtDX))
    if ("try-error" %in% class(inv.XtDX)) {
      print(summary(fit))
      stop("ERROR: matrix inversion failed")
    }

    SIGMA2   <- getSIGMA2_binary() 

    rm(Gt, XtDX, res, d, mu)
    gc()
  }else{
    stop("ERROR: outcome must be continous or binary")
  }

  # Check SIGMA2
  if (any(SIGMA2 <= 0)) stop("ERROR: values of SIGMA2 <= 0")

  #step 2: specify weight OMEGA
  #G.bar.star <- G.bar
  G.bar       <- (colSums(G))/(n)
  G.bar.star  <- G.bar/2
  OMEGA       <- dbeta(G.bar.star, shape1=a, shape2=b) 
  
  #step 3: creat input for fastASSET.R: PI_k and Z_k
  PI_Z_kGroups <- PI_kGroups(E,S,SIGMA2,OMEGA)
 
  #step 4: run  ASSET
  snps       <- "Gene"
  traits.lab <- paste("Region_", E.unique, sep="")
  omega      <- sqrt(PI_Z_kGroups$PI_k)
  alpha      <- PI_Z_kGroups$Z_k
  beta.hat   <- alpha/omega 
  sigma.hat  <- 1/omega
  cor        <- diag(length(E.unique))
  ncase      <- matrix(1,1,length(E.unique))
  ncntl      <- ncase

  # Use ASSET 
  res <- h.traits(snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=cor,
                cor.numr=FALSE, search=NULL, side=2, meta=T,
                zmax.args=NULL)
  
  return(res)
}


PI_kGroups <- function(E,S,SIGMA2,OMEGA)
{  
  p <- length(E)
  E.unique <- unique(E)
  m <- length(E.unique)
  SIGMA <- sqrt(SIGMA2)
  S.scaled <- S/SIGMA 
  
  Z_k <- rep(NA,m)
  PI_k <- rep(NA,m)
  for(j in 1:(m))
  {
    is.j <- E %in% E.unique[j]   
    S.scaled.j <- S.scaled[is.j]
    SIGMA.j <- SIGMA[is.j]
    OMEGA.j <- OMEGA[is.j]
    
    OS.j <- OMEGA.j*SIGMA.j
    SumOS_j2 <- sum(OS.j^2)
    Z_k[j] <- sum(OS.j/sqrt(SumOS_j2)*S.scaled.j)
    
    PI_k[j] <- SumOS_j2
  }
  PI_k <- PI_k/sum(OMEGA^2*SIGMA2)
  return (list(PI_k=PI_k,Z_k=Z_k))  
}

simultion1 <- function(n,p,m,prop,beta0,beta,E)
{
  N <- round((n/2)*1.2 / (exp(beta0)/(1+exp(beta0)))) #cohorat size
  
  #Simulate design matrix
  GCohort <- matrix(0,N,p)
  
  for(j in 1:p)
  {
    MAF <- prop*j
    GCohort[,j] <- rbinom(N, 2, MAF)
  }
  
  #GCohort[GCohort==2] = 1
  
  eta <- beta0+GCohort %*% beta
  
  pCase <- exp(eta)/(1+exp(eta))
  
  YCohort <- matrix(rbinom(N,1,pCase),N,1)
  
  casePool <- which(YCohort==1)
  controlPool <- which(YCohort==0)
  
  case <- sample(casePool, n/2)
  control <-sample(controlPool,n/2) 
  
  G <- rbind(GCohort[case,],GCohort[control,])
  Y <- YCohort[c(case,control),]
  return(list(Y=Y, G=G))
}