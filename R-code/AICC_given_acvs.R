### version of AICC function that operates on ACVS ...

AICC_given_acvs <- function(ts_in,acvs,N_parms,center=TRUE)
{
    ## NOTE: ASSUMES ts_in AND acvs HAVE SAME LENGTH
    ts <- if(center) ts_in - mean(ts_in) else ts_in
    N <- length(ts)
    LD_stuff <- LD_recursions(acvs)
    pacs <- LD_stuff$pacs
    rs <- rep(1,N)
    for(k in 1:(N-1)) rs[k] <- 1/prod(1-pacs[k:(N-1)]^2)
    ss <- ts[1]^2/rs[1]
    for(j in 2:N)
    {
        coeffs <- LD_stuff$blpc[[j-1]]
        ss <- ss + (ts[j] - sum(ts[(j-1):1]*coeffs))^2/rs[j]
    }
    return(list(AICC=N + N*log(2*pi/N) + N*log(ss) + sum(log(rs)) + 2*(N_parms+1)*N/(N-N_parms-2),
                rs=rs,
                pacs=pacs,
                coeffs=coeffs,
                LD_stuff=LD_stuff))
}
