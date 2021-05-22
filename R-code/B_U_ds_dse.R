###

B_U_ds_dse <- function(g_m,taper,delta_t=1)
{
    N <- length(g_m)
    h_star_h_sq <- (Re(inverse_dft(abs(dft(c(taper,rep(0,N))))^2))[1:N])^2
    return(1/(delta_t*(1 + 2*sum(abs(dft(g_m)[-1])^2 * h_star_h_sq[-1]))))
}
