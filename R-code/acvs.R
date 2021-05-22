### compute sample acvs using FFT algorithm (with option to compute
### acvs corresponding to a direct spectral estimator other than
### the periodogram by specifying a taper)

acvs <- function(ts, center=TRUE, acs=FALSE, unbiased=FALSE, taper=NULL)
{
    if(center) ts <- ts - mean(ts)
    N <- length(ts)
    taper_p <- !is.null(taper)
    if(!taper_p) taper <- rep(1/sqrt(N),N)
    results <- Re(fft(abs(fft(c(taper*ts,rep(0,N))))^2)[1:N])/(2*N)
    if(unbiased) results <- results * N /(N:1)
    if(acs) results <- results/results[1]
    return(list(acvs=results,lags=0:(N-1),center=center,acs=acs,unbiased=unbiased,taper=taper_p))
}
