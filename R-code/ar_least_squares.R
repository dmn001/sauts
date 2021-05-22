### forward least squares 

ar_fls <- function(ts, p=1, center=TRUE)
{
    if(center) ts <- ts - mean(ts)
    N <- length(ts)
    Y <- ts[(p+1):N]
    X <- matrix(nrow=N-p,ncol=p)
    for(i in 1:p) X[,i] <- ts[(p+1-i):(N-i)]
    fls_lm <- lm(Y ~ X -1)
    return(list(coeffs=coef(fls_lm),
                innov_var=ls_innov_var_est(fls_lm)))
}

### backward least squares 

ar_bls <- function(ts, p=1, center=TRUE)
{
    if(center) ts <- ts - mean(ts)
    N <- length(ts)
    Y <- ts[1:(N-p)]
    X <- matrix(nrow=N-p,ncol=p)
    for(i in 1:p) X[,i] <- ts[(1+i):(N-p+i)]
    bls_lm <- lm(Y ~ X -1)
    return(list(coeffs=coef(bls_lm),
                innov_var=ls_innov_var_est(bls_lm)))
}

### forward/backward least squares 

ar_fbls <- function(ts, p=1, center=TRUE)
{
    if(center) ts <- ts - mean(ts)
    N <- length(ts)
    Y <- c(ts[(p+1):N],ts[1:(N-p)])
    X <- matrix(nrow=2*(N-p),ncol=p)
    for(i in 1:p) X[,i] <- c(ts[(p+1-i):(N-i)],ts[(1+i):(N-p+i)])
    fbls_lm <- lm(Y ~ X -1)
    return(list(coeffs=coef(fbls_lm),
                innov_var=ls_innov_var_est(fbls_lm)))
}

###

ls_innov_var_est <- function(lm_obj) sum(resid(lm_obj)^2)/lm_obj$df
