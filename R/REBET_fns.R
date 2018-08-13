#REBET function:
#Y: vector of Phenotye
#G: n*p matrix of Genotype
#E: sub-region annotation
#X: matrix of covariates
#a,b: parameters of beta distribution

rebet <- function(response, genotypes, subRegions, responseType=NULL, 
        covariates=NULL, shape1=1, shape2=1, saveMem=FALSE)
{
    #need add some sanity check of parameter of REBET
    if (!is.matrix(genotypes)) {
        stop("ERROR: genotypes must be a matrix with at least two columns")
    }
    n <- nrow(genotypes)
    p <- ncol(genotypes)
    if ((is.null(p)) || (p < 2)) {
        stop("ERROR: genotypes must be a matrix with at least two columns")
    }
    E.unique <- unique(subRegions)
    if (p != length(subRegions)) {
        stop("ERROR: length of subRegions and number of col of genotypes must be the same")
    }
    if (shape1 <= 0) stop("ERROR: shape1 must be positive")
    if (shape2 <= 0) stop("ERROR: shape2 must be positive")

    if (n != length(response)) {
      stop("ERROR: length of response and number of rows of genotypes must be the same")
    }
    if (!is.null(responseType)) {
        if (!(responseType %in% c("continuous","binary"))) {
            stop("ERROR: responseType type must be continous, binary or NULL")
        }
    } else {
        temp <- is.finite(response)
        if (all(response[temp] %in% 0:1)) {
          responseType <- "binary"
        } else {
          responseType <- "continuous"
        }
    }
    if (is.null(covariates)) covariates <- matrix(data=1, nrow=n, ncol=1)
    if (is.vector(covariates)) covariates <- matrix(covariates, ncol=1)
    if ((!is.matrix(covariates)) && (!is.vector(covariates))) {
        stop("ERROR: covariates must be a matrix or vector")
    }
    if (n != nrow(covariates)) {
        stop("ERROR: length of response and number of rows of covariates must be the same")
    }

    # Remove missing values
    temp <- (!is.finite(response)) | (rowSums(!is.finite(genotypes)) > 0) |
        (rowSums(!is.finite(covariates)) > 0)  
    if (any(temp)) {
        response <- response[!temp]
        genotypes <- genotypes[!temp, , drop=FALSE]
        covariates <- covariates[!temp, , drop=FALSE]
        n <- nrow(genotypes)
    }
    if (!n) stop("ERROR: No data to process")

    # Add intercept if needed
    if (all(colSums(covariates == 1) != n)) covariates <- as.matrix(cbind(1, covariates))
    ncov <- ncol(covariates)

    

    #step 1:calculate score statistics vector and its variance matrix
    if (responseType == "continuous"){
    
        fit          <- try(lm(response ~ covariates - 1))
        if ("try-error" %in% class(fit)) {
            stop("ERROR: Linear regression failed")
        }
        tilde.sigma2 <- summary(fit)$sigma^2
        res          <- fit$resid   
        Gt           <- t(genotypes)
        S            <- Gt %*% res 
    
        XtX          <- t(covariates)%*%covariates
        inv.XtX      <- try(solve(XtX))
        if ("try-error" %in% class(inv.XtX)) {
            print(summary(fit))
            stop("ERROR: matrix inversion failed")
        }
    
        SIGMA2       <- .getSIGMA2_cont(saveMem, p, genotypes, n, covariates, ncov, 
                        inv.XtX, tilde.sigma2, Gt) 

        rm(Gt, XtX, res)
        gc()
    }else if (responseType == "binary"){
        fit <- try(glm(formula = response ~ covariates - 1, family = binomial(link = logit)))
        if (("try-error" %in% class(fit)) || (!fit$converged)) {
            stop("ERROR: Logistic regression failed")
        }
        mu       <- fit$fitted.value
        res      <- response - mu
        Gt       <- t(genotypes)
        S        <- Gt %*% res 
        d        <- mu * (1 - mu)
        XtDX     <- .getXtDX(saveMem, ncov, covariates, d, n)
        inv.XtDX <- try(solve(XtDX))
        if ("try-error" %in% class(inv.XtDX)) {
            print(summary(fit))
            stop("ERROR: matrix inversion failed")
        }

        SIGMA2   <- .getSIGMA2_binary(saveMem, p, genotypes, n, covariates, ncov, 
                    inv.XtDX, d, Gt) 

        rm(Gt, XtDX, res, d, mu)
        gc()
    } else{
        stop("ERROR: responseType must be continous or binary")
    }

    # Check SIGMA2
    if (any(SIGMA2 <= 0)) stop("ERROR: values of SIGMA2 <= 0")

    #step 2: specify weight OMEGA
    #G.bar.star <- G.bar
    G.bar       <- (colSums(genotypes))/(n)
    G.bar.star  <- G.bar/2
    OMEGA       <- dbeta(G.bar.star, shape1=shape1, shape2=shape2) 
  
    #step 3: creat input for fastASSET.R: PI_k and Z_k
    PI_Z_kGroups <- .PI_kGroups(subRegions,S,SIGMA2,OMEGA)
 
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
    res <- h.traits(snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl, 
        cor=cor, cor.numr=FALSE, search=NULL, side=2, meta=TRUE, 
        zmax.args=NULL)
  
    return(res)
}


