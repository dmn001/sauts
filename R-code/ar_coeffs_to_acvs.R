### NOTE: requires step_down_LD_recursions

ar_coeffs_to_acvs <- function(coeffs, max_lag, variance=1, process_variance_p=TRUE, list_of_lower_order_phi=NULL, list_of_lower_order_pev=NULL)
{
    p <- length(coeffs)
    temp <- step_down_LD_recursions(coeffs, variance=variance, process_variance_p=process_variance_p)
    acvs <- rep(0,max_lag+1)
    acvs[1] <- temp$pe_vars[1]
    if (max_lag > 0)
    {
        for (k in 1:(min(p, max_lag)))
        {
            acvs[k+1] <- temp$pe_vars[k]*temp$coeffs[[k]][k]
            if(k>1)
                acvs[k+1] <- acvs[k+1] + as.vector(temp$coeffs[[k-1]]%*%acvs[k:2])
        }
        if (max_lag > p)
        {
            tmp <- temp$coeffs[[p]]
            for (k in (p+1):max_lag)
                acvs[k+1] <- as.vector(tmp %*% acvs[k:(k-p+1)])
        }
    }
    return(acvs)
}
