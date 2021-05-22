### 

J_multiple_freq <- function(ts,the_taper,center=TRUE,recenter=FALSE,normalize=FALSE,pad_factor=1,delta_t=1)
{
    tapered_ts <- create_tapered_series(ts,the_taper,center=center,recenter=recenter,normalize=normalize)
    N <- length(tapered_ts)
    N_padded <- round(N * pad_factor)
    ts_padded <- c(tapered_ts,rep(0,N_padded-N))
    M <- floor(N_padded/2) # number of nonzero frequencies
    freqs <- (0:M)/(N_padded*delta_t)
    J_dft <- sqrt(delta_t)*fft(ts_padded)[1:(M+1)]
    return(list(freqs=freqs,
                J_dft=J_dft))
}

### needed?

J_one_freq <- function(freq,ts,the_taper,center=TRUE,recenter=FALSE,normalize=FALSE,delta_t=1)
{
    tapered_ts <- create_tapered_series(ts,the_taper,center=center,recenter=recenter,normalize=normalize)
    N <- length(tapered_ts)
    return(sqrt(delta_t)*sum(exp(complex(real=rep(0,N),imag=(-2*pi*freq*(0:(N-1))*delta_t)))*tapered_ts))
}
