### compute SDF for ARMA process

arma_coeffs_to_sdf <- function(ar_coeffs, ma_coeffs, innov_var=1, pad_factor=4, delta_t=1, N_pad=NULL)
{
    p <- length(ar_coeffs)
    q <- length(ma_coeffs)
    if(is.null(N_pad)) N_pad <- 2^(ceiling(log2(max(p+1,q+1)))) * pad_factor
    return(list(freqs=(0:floor(N_pad/2))/(N_pad*delta_t),sdfe=innov_var*delta_t*abs(fft(c(1,-ma_coeffs,rep(0,N_pad-(q+1))))[1:(floor(N_pad/2)+1)])^2/abs(fft(c(1,-ar_coeffs,rep(0,N_pad-(p+1))))[1:(floor(N_pad/2)+1)])^2))
}
