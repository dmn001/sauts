### generate crisscross for WOSA with Hanning data taper and 50% overlap

do_crisscross_hanning_50_wosa <- function(the_taper,N_B,n_or_block_starts,delta_t=1)
{
    edof <- 36*N_B^2/(19*N_B-1)
    return(list(up = 10*log10(edof/qchisq(0.025,edof)),
                down = -10*log10(edof/qchisq(0.975,edof)),
                width = B_H(the_taper,delta_t),
                edof=edof))
}
