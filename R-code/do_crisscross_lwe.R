### create values needed to display bandwidth/CI crisscross
### for a lag window SDF estimate
###
### NOTE: requires B_U

do_crisscross_lwe <- function(lw_1_sided,taper=NULL,delta_t=1)
{
    N <- length(lw_1_sided)
    if(is.null(taper)) taper <- rep(1/sqrt(N),N)  # default taper
    edof <- 2/(sum(taper^4)*(1+ 2*sum(lw_1_sided[-1]^2)))
    h2_star_h2 <- Re(inverse_dft(abs(dft(c(taper^2,rep(0,N))))^2))[1:N]
    alt_edof <- 2/(h2_star_h2[1]+ 2*sum(lw_1_sided[-1]^2*h2_star_h2[-1]))
    list(up = 10*log10(edof/qchisq(0.025,edof)),
         down = -10*log10(edof/qchisq(0.975,edof)),
         edof = edof,
         width = B_U(lw_1_sided,taper,delta_t),
         alt_up = 10*log10(alt_edof/qchisq(0.025,alt_edof)),
         alt_down = -10*log10(alt_edof/qchisq(0.975,alt_edof)),
         alt_edof = alt_edof)
}
