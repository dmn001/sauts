### compute Fejer's kernel

fejer_kernel <- function(N, delta_t=1, N_pad=2048, tiny=10^(-5))
{
    ones_in_zeros <- c(rep(1,N),rep(0,N_pad-N))
    freqs <- (-((N_pad/2)-1):(N_pad/2))/N_pad
    x <- abs(fft(ones_in_zeros)[1:((N_pad/2)+1)])^2/N
    x[x <= tiny] <- tiny
    kernel <- c(rev(x[c(-1,-length(x))]),x)
    return(list(freqs=freqs/delta_t, kernel=kernel*delta_t))
}

###

fejer_kernel_single_freq <- function(f,N,delta_t=1)
{
    if(isTRUE(all.equal(sin(pi*f*delta_t),0)))
        delta_t*N
    else
        delta_t*(sin(N*pi*f*delta_t))^2/(N*(sin(pi*f*delta_t))^2)
}

