### NOTE: the formula in SAPA2e seems to imply that we need to
###       multiply by delta_t rather than divide, but the
###       definition of the autocorrelation of the taper
###       involves a multiplication by delta_t, which becomes
###       a square 

B_H_bar <- function(the_tapers,delta_t=1,weights=NULL)
{
    the_tapers <- as.matrix(the_tapers)
    N <- nrow(the_tapers)
    K <- ncol(the_tapers)
    if(is.null(weights)) weights <- rep(1/K,K)
    for(k in 1:K)
    {
        temp <- Re(fft(abs(fft(c(the_tapers[,k],rep(0,N)))^2)))/(2*N)
        results <- if(k==1) temp*weights[1] else results + temp*weights[k]
    }
    return(1/(sum(results^2)*delta_t))
}
