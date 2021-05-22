### NOTE: requires step_down_LD_recursions
###       requires sim_ar_process

sim_arma_process <- function(N, phis=NULL, thetas=NULL, sig2=1, LD_stuff=step_down_LD_recursions(phis,sig2,proc="FALSE"))
{
    p <- length(phis)
    q <- length(thetas)
    if(p==0 && q==0) return(rnorm(N,sd=sqrt(sig2)))
    else
    {
        if(p==0) return(filter(rnorm(N+q,sd=sqrt(sig2)),c(1,-thetas),method="conv",sides=1)[-(1:q)])
        else
        {
            ts <- sim_ar_process(N+q,LD_stuff=LD_stuff)
            if(q>0) return(filter(ts,c(1,-thetas),method="conv",sides=1)[-(1:q)])
            else
                return(ts)
        }
    }
}
