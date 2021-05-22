### create direct spectral estimate
###
### NOTE: requires create_tapered_series, do_crisscross_dse and B_H
###
### NOTE: sum of squares of the_taper is intended to be unity

direct_sdf_est <- function(ts,the_taper,center=TRUE,recenter=FALSE,normalize=FALSE,pad_factor=1,N_pad=NULL,delta_t=1,two_sided_p=FALSE,start_with_0_p=FALSE)
{
    tapered_ts <- create_tapered_series(ts,the_taper,center=center,recenter=recenter,normalize=normalize)
    N <- length(ts)
    if(is.null(N_pad)) N_pad <- round(N * pad_factor)
    ts_padded <- c(tapered_ts,rep(0,N_pad-N))
    M <- floor(N_pad/2)  # number of nonzero frequencies
    if(two_sided_p)
    {
        if(start_with_0_p)  # from zero to just below twice Nyquist
        {
            freqs <- (0:(N_pad-1))/(N_pad*delta_t)
            dse <- delta_t*abs(fft(ts_padded))^2
        }
        else # from just beyond -Nyquist to Nyquist
        {
            freqs <- c(((M+1)-N_pad):(-1),(0:M))/(N_pad*delta_t)
            temp <- abs(fft(ts_padded))^2
            dse <- delta_t*c(temp[(M+2):N_pad],temp[1:(M+1)])
        }
    }
    else # from zero to Nyquist (or just below Nyquist)
    {
        freqs <- (0:M)/(N_pad*delta_t)
        dse <- delta_t*abs(fft(ts_padded)[1:(M+1)])^2
    }
    return(list(freqs=freqs, sdfe=dse, cc=do_crisscross_dse(the_taper,delta_t=delta_t)))
}
