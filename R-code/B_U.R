### standard measure of bandwidth of spectral window
### NOTE: lw and taper (if supplied) should be vectors of length N
###       (the sample size); lw should contain the lag window
###       evaluated at lag 0, 1, ..., N-1

B_U <- function(lw_1_sided,taper=NULL,delta_t=1)
{
    N <- length(lw_1_sided)
    if(is.null(taper)) taper <- default_taper(N)
    h_star_h_sq <- (Re(inverse_dft(abs(dft(c(taper,rep(0,N))))^2))[1:N])^2
    return(1/(delta_t*(1 + 2*sum(lw_1_sided[-1]^2 * h_star_h_sq[-1]))))
}
