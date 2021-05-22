### WOSA sdf estimation

wosa_sdf_est <- function(ts,the_taper,n_or_block_starts=round(length(the_taper)/2),center=TRUE,recenter=FALSE,normalize=FALSE,pad_factor=1,delta_t=1,do_crisscross=do_crisscross_wosa)
{
    if(center) ts <- ts - mean(ts)
    N_S <- length(the_taper)
    block_starts <- if(length(n_or_block_starts) == 1)
                        seq(1,length(ts)-N_S+1,n_or_block_starts)
                    else
                        n_or_block_starts
    N_B <- length(block_starts)
    dses <- vector("list",N_B)
    N_padded <- round(N_S * pad_factor)
    M <- floor(N_padded/2)  # number of nonzero frequencies
    freqs <- (0:M)/(N_padded*delta_t)
    for(k in 1:N_B)
    {
        tapered_ts <- create_tapered_series(ts[block_starts[k]:(block_starts[k]+N_S-1)],the_taper,center=recenter)
        ts_padded <- c(tapered_ts,rep(0,N_padded-N_S))
        dse <- delta_t*abs(fft(ts_padded)[1:(M+1)])^2
        dses[[k]] <- structure(list(freqs=freqs, sdfe=dse, cc=do_crisscross_dse(the_taper,delta_t=delta_t)))
        results <- if(k==1) dse else results + dse
    }
    return(list(freqs=freqs,
                sdfe=results/N_B,
                cc=do_crisscross(the_taper,N_B,n_or_block_starts,delta_t=delta_t),
                dses=dses))
}
