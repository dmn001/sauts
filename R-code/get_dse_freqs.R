### when using sapply to collect several sdfe values from direct_sdf_est,
### freqs values are lost - the following function computes them

get_dse_freqs <- function(N,pad_factor=1,delta_t=1)
{
    N_padded <- round(N * pad_factor)
    M <- floor(N_padded/2) # number of nonzero frequencies
    return((0:M)/(N_padded*delta_t))
}
