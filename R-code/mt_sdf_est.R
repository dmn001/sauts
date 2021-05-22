### NOTE: multitapers are in the columns of the_tapers - each column should
###       have a sum of squares equal to unity

mt_sdf_est <- function(ts,the_tapers,center=TRUE,recenter=FALSE,normalize=FALSE,pad_factor=1,delta_t=1,weights=NULL)
{
    the_tapers <- as.matrix(the_tapers)
    N <- length(ts)  # should be the same as nrow(the_tapers)
    K <- ncol(the_tapers)
    if(is.null(weights)) weights <- rep(1/K,K)
    eigenspectra <- vector("list",K)
    N_padded <- round(N * pad_factor)
    M <- floor(N_padded/2) # number of nonzero frequencies
    freqs <- (0:M)/(N_padded*delta_t)
    for(k in 1:K)
    {
        tapered_ts <- create_tapered_series(ts,the_tapers[,k],center=center,recenter=recenter,normalize=normalize)
        ts_padded <- c(tapered_ts,rep(0,N_padded-N))
        dse <- delta_t*abs(fft(ts_padded)[1:(M+1)])^2
        eigenspectra[[k]] <- list(freqs=freqs,
                                  sdfe=dse,
                                  cc=do_crisscross_dse(the_tapers[,k],delta_t=delta_t))
        results <- if(k==1) dse*weights[1] else results + dse*weights[k]
    }
    return(list(freqs=freqs,
                sdfe=results,
                cc=do_crisscross_mt(the_tapers,delta_t=delta_t),
                eigenspectra=eigenspectra))
}
