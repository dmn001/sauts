### expected value of lag window estimator
###
### acvs should be a vector of length N containing a theoretical ACVS at lags 0, ..., N-1;
### taper should be a vector of length N normalized so that its sum of squares is 1;
### lag_window should be a vector of length N containing a lag window to be applied at lags 0, 1, ..., N-1
### N_pad is the length of the vector handed to the fft - this needs to be at least
### as large as the length of acvs

ev_lag_window_sdf_estimator <- function(acvs, taper={N <- length(acvs); rep(1/sqrt(N),N)}, lag_window=rep(1,length(acvs)), N_pad=2*length(acvs), delta_t=1)
{
    N_ts <- length(acvs)
    acvs_taper <- Re(fft(abs(fft(c(taper,rep(0,N_pad-N_ts)))^2)))/N_pad
    acvs_wrapped <- c(lag_window*acvs, rep(0,N_pad-2*N_ts+1),lag_window[N_ts:2]*acvs[N_ts:2])
    return(list(sdf_ev=delta_t*Re(fft(acvs_wrapped*acvs_taper))[1:((N_pad/2)+1)],
                freqs=(0:(N_pad/2))/(N_pad*delta_t),
                B_U=B_U(lag_window,taper,delta_t)))
}
