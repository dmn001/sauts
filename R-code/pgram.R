### compute periodogram over a grid of frequencies

pgram <- function(ts,pad_factor=1,delta_t=1,center=TRUE,two_sided_p=FALSE,start_with_0_p=FALSE)
{
    N <- length(ts)
    N_padded <- round(N * pad_factor)
    ts_padded <- rep(0,N_padded)
    ts_padded[1:N] <- if(center) ts-mean(ts) else ts
    M <- floor(N_padded/2) # number of nonzero frequencies
    if(two_sided_p)
    {
        if(start_with_0_p)  # from zero to just below twice Nyquist
        {
            freqs <- (0:(N_padded-1))/(N_padded*delta_t)
            pgram <- delta_t*abs(fft(ts_padded))^2/N
        }
        else # from just beyond -Nyquist to Nyquist
        {
            freqs <- c(((M+1)-N_padded):(-1),(0:M))/(N_padded*delta_t)
            temp <- abs(fft(ts_padded))^2
            pgram <- delta_t*c(temp[(M+2):N_padded],temp[1:(M+1)])/N
        }
    }
    else # from zero to Nyquist (or just below Nyquist)
    {
        freqs <- (0:M)/(N_padded*delta_t)
        pgram <- delta_t*abs(fft(ts_padded)[1:(M+1)])^2/N
    }
    return(list(freqs=freqs, sdfe=pgram, cc=list(up=10*log10(2/qchisq(0.025,2)),down=-10*log10(2/qchisq(0.975,2)),edof=2,width=3*N/((2*N^2+1)*delta_t))))
}
