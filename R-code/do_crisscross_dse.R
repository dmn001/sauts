### create values needed to display bandwidth/CI crisscross
### for a direct spectral estimate
###
### NOTE: requires B_H

do_crisscross_dse <- function(taper,delta_t=1)
{
    edof <- 2
    list(up = 10*log10(edof/qchisq(0.025,edof)),
         down = -10*log10(edof/qchisq(0.975,edof)),
         edof = edof,
         width = B_H(taper,delta_t))
}
