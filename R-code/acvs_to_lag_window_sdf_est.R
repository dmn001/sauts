### 

acvs_to_lag_window_sdf_est <- function(acvs_est, m=NULL, taper_parameter=NULL, taper=default_taper, lag_window=default_lag_window, delta_t=1, pad_factor=1, two_sided_p=TRUE)
{
    if(two_sided_p)
    {
        acvs_est_2_sided <- acvs_est  # length is 2*N -1 
        N <- (length(acvs_est_2_sided)+1)/2
    }
    else
    {
        N <- length(acvs_est)
        acvs_est_2_sided <- c(rev(acvs_est[-1]),acvs_est)
    }
    the_taper <- if(is.numeric(taper)) taper else taper(N,taper_parameter)
    Nm1 <- N - 1
    lw <- if(is.numeric(lag_window)) lag_window else sapply((-Nm1):Nm1, lag_window, m)
    N_dft <- round(2*N*pad_factor)
    M <- floor(N_dft/2)  # number of nonzero frequencies
    for_dft <- rep(0,N_dft)
    lw_acvs <- lw * acvs_est_2_sided
    for_dft[1:N] <- lw_acvs[N:(2*N-1)]
    for_dft[(N_dft-N+2):N_dft] <- lw_acvs[1:Nm1]
    list(freqs=(0:M)/(delta_t*N_dft),
         sdfe=delta_t*Re(dft(for_dft))[1:(M+1)],
         cc=do_crisscross_lwe(lw[-(1:Nm1)],the_taper,delta_t=delta_t),
         lw_2_sided=lw,
         acvs=lw_acvs[-(1:(N-1))])
}
