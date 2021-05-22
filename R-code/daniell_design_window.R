###

daniell_design_window <- function(m,delta_t=1)
{
    list(freqs=c(0,rep(1/(2*m*delta_t),2)),
         dw=c(rep(m*delta_t,2),0))
}
