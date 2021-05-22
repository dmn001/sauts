### expected value of DCT-based periodogram
###
### acvs should be a vector of length N containing a theoretical ACVS at lags 0, ..., N-1
###
### NOTE: requires ev_lag_window_sdf_estimator

ev_DCTII <- function(the_acvs)
{
    N <- length(the_acvs)
    freqs <- (0:(N-1))/(2*N)
    results <- ev_lag_window_sdf_estimator(the_acvs,rep(1/sqrt(N),N),N_pad=2*N)$sdf_ev[-(N+1)]
    results[1] <- the_acvs[1] + 2*sum(the_acvs[-1]*(1-(1:(N-1))/N))
    for(k in 1:(N-1))
        results[k+1] <- results[k+1] - 2*sum(the_acvs[-1]*sin(pi*(1:(N-1))*k/N))/(sin(pi*k/N)*N)
    return(list(sdf_ev=results,freqs=freqs))
}

ev_lag_window_sdf_estimator <- function(acvs, taper={N <- length(acvs); rep(1/sqrt(N),N)}, lag_window=rep(1,length(acvs)), N_pad=2*length(acvs), delta_t=1)
{
    N_ts <- length(acvs)
    acvs_taper <- Re(fft(abs(fft(c(taper,rep(0,N_pad-N_ts)))^2)))/N_pad
    acvs_wrapped <- c(lag_window*acvs, rep(0,N_pad-2*N_ts+1),lag_window[N_ts:2]*acvs[N_ts:2])
    return(list(sdf_ev=delta_t*Re(fft(acvs_wrapped*acvs_taper))[1:((N_pad/2)+1)],freqs=(0:(N_pad/2))/(N_pad*delta_t)))
}
