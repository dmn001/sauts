###

ar_reduced_likelihood <- function(ts,ar_coeffs,center=TRUE)
{
    if(center)  ts <- ts - mean(ts)
    n <- length(ts)
    var_ts <- sum(ts^2)/n
    p <- length(ar_coeffs)
    LD_stuff <- step_down_LD_recursions(ar_coeffs,var=var_ts)
    pacs <- rep(0,p)
    for(k in 1:p) pacs[k] <- LD_stuff$coeffs[[k]][k]
    rs <- rep(1,n)
    for(k in 1:p) rs[k] <- 1/prod(1-pacs[k:p]^2)
    ss <- ts[1]^2/rs[1]
    for(j in 2:n)
    {
        if(j < p+2) coeffs <- LD_stuff$coeffs[[j-1]]
        n_coeffs <- length(coeffs)
        ss <- ss + (ts[j] - sum(ts[(j-1):(j-n_coeffs)]*coeffs))^2/rs[j]
    }
    return(n + n*log(2*pi/n) + n*log(ss) + sum(log(rs)))
}
