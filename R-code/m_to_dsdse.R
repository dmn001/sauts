###

m_to_dsdse <- function(my_m, my_dft_dse, hanning_taper_p=FALSE, delta_t=1)
{
    N_dse <- length(my_dft_dse)
    ## NOTE: currently limited to default or Hanning tapers
    ##       due to edof calculation, which hasn't been worked out
    ##       in general; note that B_U has been worked
    ##       out for tapers in general
    the_taper <- if(hanning_taper_p) hanning_taper(N_dse) else default_taper(N_dse)
    temp <- dnorm(0:(N_dse/2-1),sd=1/(sqrt(2)*my_m))
    g_raw <- c(temp,0,rev(temp[-1]))
    g_m <- g_raw/sum(g_raw)
    edofs <- if(hanning_taper_p) 2/(sum(g_m^2) + sum(g_m*circular_shift(g_m,1))) else 2/sum(g_m^2)
    dft_g_m <- dft(g_m)
    return(list(freqs=(0:(N_dse-1))/N_dse,
                sdfe=Re(inverse_dft(my_dft_dse*dft_g_m)),
                cc=list(up = 10*log10(edofs/qchisq(0.025,edofs)),
                        down = -10*log10(edofs/qchisq(0.975,edofs)),
                        edof = edofs,
                        width = B_U_ds_dse(g_m,the_taper,delta_t=delta_t))))
}
