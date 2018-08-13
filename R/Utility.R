.PI_kGroups <- function(E,S,SIGMA2,OMEGA)
{  
    p <- length(E)
    E.unique <- unique(E)
    m <- length(E.unique)
    SIGMA <- sqrt(SIGMA2)
    S.scaled <- S/SIGMA 
  
    Z_k <- rep(NA,m)
    PI_k <- rep(NA,m)
    for(j in seq_len(m))
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

# Local function to compute diagonal of SIGMA for continuous case
.getSIGMA2_cont <- function(saveMem, p, G, n, X, ncov, inv.XtX, tilde.sigma2, 
    Gt) {

    if (saveMem) {
        ret <- double(p)
        tmp <- .C("SIGMA2_cont", as.numeric(G), as.integer(n), as.integer(p), 
               as.numeric(X), as.integer(ncov), as.numeric(inv.XtX),
               as.numeric(tilde.sigma2), ret=ret, PACKAGE="REBET")
        ret <- tmp$ret 
    } else {
        SIGMA <- Gt%*% (diag(1, n)-X %*% inv.XtX %*% t(X)) %*% G * tilde.sigma2
        ret   <- diag(SIGMA)
    }
    ret

} # END: getSIGMA2_cont

# Local function to compute diagonal of SIGMA for binary case
.getSIGMA2_binary <- function(saveMem, p, G, n, X, ncov, inv.XtDX, d, Gt) {

        if (saveMem) {
            ret <- double(p)
            tmp <- .C("SIGMA2_binary", as.numeric(G), as.integer(n), 
                as.integer(p), as.numeric(X), as.integer(ncov), 
                as.numeric(inv.XtDX),
                as.numeric(d), ret=ret, PACKAGE="REBET")
            ret <- tmp$ret 
        } else {
            SIGMA <- Gt%*%(diag(d)-(d %o% d) * (X %*% (inv.XtDX) %*% t(X)))%*%G
            ret   <- diag(SIGMA)
        }
        ret

} # END: getSIGMA2_binary

# Local function to compute XtDX <- t(X)%*%D%*%X for binary case
.getXtDX <- function(saveMem, ncov, X, d, n) {

        if (saveMem) {
            ret <- double(ncov*ncov)
            tmp <- .C("computeXtDX", as.numeric(X), as.numeric(d), 
                as.integer(n), as.integer(ncov), ret=ret, PACKAGE="REBET")
            ret <- matrix(tmp$ret, nrow=ncov, ncol=ncov)
        } else {
            ret <- t(X)%*%diag(d)%*%X
        }
        ret
  
} # END: getXtDX