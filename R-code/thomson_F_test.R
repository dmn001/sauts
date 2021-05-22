###

thomson_F_test <- function(ts,the_tapers,pad_factor=1,delta_t=1,...)
{
    K_tapers <- ncol(the_tapers)
    H_at_zeros <- sapply(1:K_tapers,function(k) delta_t*sum(the_tapers[,k]))
    H_at_zeros[seq(2,K_tapers,2)] <- 0
    sum_H_at_zeros_squared <- sum(H_at_zeros^2)
    J_k <- sapply(1:K_tapers,function(k) J_multiple_freq(ts,the_tapers[,k],pad_factor=pad_factor,delta_t=delta_t,...)$J_dft)
    C_hat <- sqrt(delta_t)*rowSums(t(t(J_k)*H_at_zeros))/sum_H_at_zeros_squared
    J_k_hat <- (matrix(C_hat,ncol=1) %*% H_at_zeros)/sqrt(delta_t)  
    denom <- delta_t*rowSums(sapply(1:K_tapers,function(k) abs(J_k[,k] - J_k_hat[,k])^2))
    F_test <- (K_tapers-1)*abs(C_hat)^2*sum_H_at_zeros_squared/denom
    freqs <- get_dse_freqs(length(ts),pad_factor=pad_factor,delta_t=delta_t)
    return(list(F_test=F_test,
                freqs=freqs,
                C_hat=C_hat,
                J_k=J_k,
                J_k_hat=J_k_hat,
                the_tapers=the_tapers,
                delta_t=delta_t,
                pad_factor=pad_factor,
                K_tapers=K_tapers,
                N_ts=nrow(the_tapers)))
}
