### create crisscross for WOSA SDF esimate
###
### NOTE: appropriate for n, but an approximation for block starts

do_crisscross_wosa <- function(the_taper,N_B,n_or_block_starts,delta_t=1)
{
    n <- round(mean(n_or_block_starts))
    N_S <- length(the_taper)
    the_denom <- 1
    if(N_B > 1) for(m in 1:(N_B-1))
                {
                    if(m*n > N_S-1) break()
                    the_denom <- the_denom + 2*(1-m/N_B)*N_S*acvs(the_taper,center=FALSE)$acvs[m*n+1]
                }
    edof <- 2*N_B/the_denom
    return(list(up = 10*log10(edof/qchisq(0.025,edof)),
         down = -10*log10(edof/qchisq(0.975,edof)),
         width = B_H(the_taper,delta_t),
         edof=edof))
}
