###

reshaped_mt_sdfe <- function(k_reshape,F_test_goodies,mt_est,k_span=NULL)
{
    k_freq <- k_reshape + 1
    the_tapers <- F_test_goodies$the_tapers
    K_tapers <- ncol(the_tapers)
    N_ts     <- nrow(the_tapers)
    N_fft    <- N_ts * F_test_goodies$pad_factor
    delta_t <- F_test_goodies$delta_t
    H_ks <- sapply(1:K_tapers,function(k) delta_t*fft(c(the_tapers[,k],rep(0,N_fft-N_ts))))
    J_k <- F_test_goodies$J_k
    C_hat <- F_test_goodies$C_hat
    freqs <- F_test_goodies$freqs
    if(is.null(k_span)) k_span <- max(round(1.25*mt_est$cc$width*N_fft*delta_t/2) - 1,1)
    return(list(freqs=freqs[(k_freq-k_span-1):(k_freq+k_span+1)],
                sdfe=c(mt_est$sdfe[k_freq-k_span-1],
                       sapply((k_freq-k_span):(k_freq+k_span), function(j) sum(sapply(1:K_tapers, function(k) abs(J_k[j,k] - C_hat[k_freq]*(if(j<k_freq) Conj(H_ks[k_freq+1-j,k]) else H_ks[j-k_freq+1,k])/sqrt(delta_t))^2)))/K_tapers,
                       mt_est$sdfe[k_freq+k_span+1]),
                k_span=k_span))
}

