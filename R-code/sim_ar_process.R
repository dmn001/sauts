### NOTE: requires step_down_LD_recursions

sim_ar_process <- function(N, phis=NULL, sig2=1, LD_stuff=step_down_LD_recursions(phis,sig2,proc="FALSE"), innovations=rnorm(N))
{
    vars <- LD_stuff$pe_vars
    coeffs <- LD_stuff$coeffs
    p <- length(vars) - 1
    ts <- rep(0,N)
    ts[1] <- innovations[1]*sqrt(vars[1])
    if(p>1)
    {
        for (n in 2:p)
        {
            ts[n] <- as.numeric(ts[(n-1):1] %*% coeffs[[n-1]]) + innovations[n]*sqrt(vars[n])
        }
    }
    innov_sd <- sqrt(vars[p+1])
    ar_coeffs <- coeffs[[p]]
    for (n in (p+1):N)
    {
        ts[n] <- as.numeric(ts[(n-1):(n-p)] %*% ar_coeffs) + innovations[n]*innov_sd
    }
    return(ts)
}

