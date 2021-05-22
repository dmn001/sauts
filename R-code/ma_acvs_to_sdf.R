### 

ma_acvs_to_sdf <- function(acvs_1_sided, delta_t=1)
{
    acvs_2_sided <- c(rev(acvs_1_sided[-1]),acvs_1_sided)
    N <- (length(acvs_2_sided)+1)/2
    Nm1 <- N - 1
    N_dft <- 2*N
    for_dft <- rep(0,N_dft)
    for_dft[1:N] <- acvs_2_sided[N:(2*N-1)]
    for_dft[(N_dft-N+2):N_dft] <- acvs_2_sided[1:Nm1]
    return(list(freqs=(0:(N_dft/2))/N_dft,
                sdf=Re(dft(for_dft))[1:(N_dft/2 + 1)]))
}
