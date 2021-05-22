### requires next_power_of_two
### requires dft
### requires beta_W_swb_exact
### requires B_W_swb_exact

lag_window_to_smoothing_window <- function(w,N_pad=2*next_power_of_2(length(w)),delta_t=1)
{
    N_w <- length(w)
    for_dft <- rep(0,N_pad)
    for_dft[1:N_w] <- w
    for_dft[(N_pad-N_w+2):N_pad] <- rev(w[-1])
    N_freqs <- trunc(N_pad/2) + 1
    squared_beta_W <- beta_W_swb_exact_squared(w,delta_t=delta_t)
    return(list(freqs=(0:(N_freqs-1))/(N_pad*delta_t),
                sw=delta_t*Re(dft(for_dft)[1:N_freqs]),
                squared_beta_W=squared_beta_W,
                beta_W=if(squared_beta_W>=0) sqrt(squared_beta_W) else NA,
                B_W=B_W_swb_exact(w,delta_t=delta_t)))
}

