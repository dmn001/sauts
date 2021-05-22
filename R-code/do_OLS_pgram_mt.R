do_OLS_pgram_mt <- function(d,N_rep=10000,seed=42,N_ts=1024,low_freq=0,high_freq=1/6,N_chop_mt=7,K=5,use_freqs_p=TRUE)
{
    set.seed(seed)
    fd_setup <- sim_circ_embed_setup(N_ts, gen_FD_ACVS_given_d, d)
    the_sine_tapers <- taper("sine",N_ts,n.taper=K)
    all_freqs <- (0:(N_ts/2))/N_ts
    N_all_freqs <- length(all_freqs)
    pgram_ave <- rep(0,N_all_freqs)
    mt_ave <- rep(0,N_all_freqs)
    gg_indices <- which(all_freqs > low_freq & all_freqs <= high_freq)
    i_range <- gg_indices[1]:gg_indices[length(gg_indices)]
    N_reg_freqs <- length(i_range)
    log_freqs_for_reg <- if(use_freqs_p) log(all_freqs[i_range]) else log(2*(sin(pi*all_freqs[i_range])))
    slope_ests_pgram <- vector("numeric",N_rep)
    slope_ests_mt <- matrix(nrow=N_rep,ncol=N_chop_mt+1)
    for(n in 1:N_rep)
    {
        ts_fd <- sim_circ_embed_from_setup(fd_setup)
        ts_fd_pgram <- pgram(ts_fd)
        pgram_ave <- pgram_ave + ts_fd_pgram$sdfe
        ts_fd_mt <- mt_sdf_est(ts_fd,the_sine_tapers)
        mt_ave <- mt_ave + ts_fd_mt$sdfe
        slope_ests_pgram[n] <- coef(lm(log(ts_fd_pgram$sdfe[i_range]) ~ log_freqs_for_reg))[2]
        slope_ests_mt[n,1] <- coef(lm(log(ts_fd_mt$sdfe[i_range]) ~ log_freqs_for_reg))[2]
        if(N_chop_mt > 0)
            for(k in 2:(N_chop_mt+1))
                slope_ests_mt[n,k] <- coef(lm(log(ts_fd_mt$sdfe[i_range[-(1:(k-1))]]) ~ log_freqs_for_reg[-(1:(k-1))]))[2]
    }
    d_ests_pgram <- -slope_ests_pgram/2
    d_ests_mt <- -slope_ests_mt/2
    return(list(freqs=all_freqs,
                pgram_ave=pgram_ave/N_rep,
                mt_ave=mt_ave/N_rep,
                d_ests_pgram=d_ests_pgram,
                d_ests_mt=d_ests_mt,
                d_mean_pgram=mean(d_ests_pgram),
                d_mean_mt=colMeans(d_ests_mt),
                d_mse_pgram=mean((d_ests_pgram-d)^2),
                d_mse_mt=colMeans((d_ests_mt-d)^2)))
}
