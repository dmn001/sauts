### compute lag window SDF estimate using time serie

ts_to_lag_window_sdf_est <- function(ts, m=NULL, taper_parameter=NULL, taper=default_taper, lag_window=default_lag_window, center=TRUE, recenter=FALSE, normalize=FALSE, delta_t=1, pad_factor=1)
{
    N <- length(ts)
    the_taper <- if(is.numeric(taper)) taper else taper(N,taper_parameter)
    ts_for_dft <- create_tapered_series(ts,the_taper,center=center,recenter=recenter,normalize=normalize)
    acvs_est <- acvs(ts_for_dft,center=FALSE)
    acvs_est_2_sided <- N*c(rev(acvs_est$acvs[-1]),acvs_est$acvs)
    Nm1 <- N - 1
    lw <- if(is.numeric(lag_window)) c(rev(lag_window[-1]),lag_window) else sapply((-Nm1):Nm1, lag_window, m)
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
