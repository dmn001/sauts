###

B_U_ds_pgram <- function(g_m,delta_t=1)
{
    N <- length(g_m)
    return(1/(delta_t*(1 + 2*sum(abs(dft(g_m)[-1])^2 * (1 - (1:(N-1))/N)^2))))
}
