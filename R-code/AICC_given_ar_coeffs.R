###

AICC_given_ar_coeffs <- function(ts_in,ar_coeffs,center=TRUE)
{
    ts <- if(center) ts_in - mean(ts_in) else ts_in
    N <- length(ts)
    var_ts <- sum(ts^2)/N
    p <- length(ar_coeffs)
    LD_stuff <- step_down_LD_recursions(ar_coeffs,var=var_ts)
    pacs <- rep(0,p)
    for(k in 1:p) pacs[k] <- LD_stuff$coeffs[[k]][k]
    rs <- rep(1,N)
    for(k in 1:p) rs[k] <- 1/prod(1-pacs[k:p]^2)
    ss <- ts[1]^2/rs[1]
    for(j in 2:N)
    {
        if(j < p+2) coeffs <- LD_stuff$coeffs[[j-1]]
        N_coeffs <- length(coeffs)
        ss <- ss + (ts[j] - sum(ts[(j-1):(j-N_coeffs)]*coeffs))^2/rs[j]
    }
    return(list(AICC=N + N*log(2*pi/N) + N*log(ss) + sum(log(rs)) + 2*(p+1)*N/(N-p-2),
                rs=rs,
                pacs=pacs,
                N_coeffs=N_coeffs,
                coeffs=coeffs,
                LD_stuff=LD_stuff))
}
