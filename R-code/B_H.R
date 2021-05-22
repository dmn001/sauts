### standard measure of effective bandwidth of a direct spectral
### estimator (based on autocorrelation width of spectral window)
###
### NOTE: the formula in SAUTS seems to imply that we need to
###       multiply by delta.t rather than divide, but the
###       definition of the autocorrelation of the taper
###       involves a multiplication by delta_t, which becomes
###       a square 

B_H <- function(taper,delta_t=1)
{
    N <- length(taper)
    return(1/(sum((Re(fft(abs(fft(c(taper,rep(0,N)))^2)))/(2*N))^2)*delta_t))
}
