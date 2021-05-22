##########################################################################
### adaptation of gcv_dsp to work with Parzen lag window estimates
### (but with just the default data taper) rather than discretely
### smoothed periodogram estimates
###
### NOTE: input is the time series, not the periodogram
##########################################################################

gcv_lw <- function(m,ts,keep_0_p=TRUE,swap_p=FALSE,lag_window=parzen_lag_window,delta_t=1)
{
    N <- length(ts)
    Nyquist_p <- is_even(N)
    ##
    M_chi_2 <- floor((N-1)/2)  # number of Fourier frequencies > 0 and < Nyquist
    M <- M_chi_2 + (if(keep_0_p) 1 else 0) + (if(Nyquist_p) 1 else 0)
    ##
    temp <- pgram(ts,delta_t=delta_t)
    freqs <- temp$freqs
    M_with_0 <- length(freqs)  # M_with_0 is same as M if keep_0_p is TRUE
    pgram_ts <- temp$sdfe
    lw_est_all <- ts_to_lag_window_sdf_est(ts,m=m,center=TRUE,lag_window=lag_window,delta_t=delta_t)
    edofs <- lw_est_all$cc$edof
    B_U <- lw_est_all$cc$width
    g_0 <- mean(lw_est_all$lw_2_sided)
    dofs_model <- (1-2*g_0)^2
    lw_est <- lw_est_all$sdfe[seq(1,by=2,length.out=M_with_0)]
    temp <- if(swap_p) ((lw_est/pgram_ts) - log(lw_est/pgram_ts) - 1)/dofs_model else ((pgram_ts/lw_est) - log(pgram_ts/lw_est) -1 )/dofs_model
    gcv <- (sum(temp[2:(M_chi_2+1)]) +
            if(keep_0_p) 0.5*temp[1] else 0 + if(Nyquist_p) 0.5*temp[N/2+1] else 0)/M_with_0
    return(list(N=N,
                Nyquist_p=Nyquist_p,
                M_chi_2=M_chi_2,
                M=M,
                M_with_0=M_with_0,
                edofs=edofs,
                B_U=B_U,
                dofs_model=dofs_model,
                cc=lw_est_all$cc,
                lw_2_sided=lw_est_all$lw_2_sided,
                full_lw_sdfe=lw_est_all$sdfe,
                full_lw_freqs=lw_est_all$freqs,
                freqs=freqs,
                pgram_ts=pgram_ts,
                sdfe=lw_est,
                gcv=gcv))
}
