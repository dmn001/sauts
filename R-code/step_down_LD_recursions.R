### step-down Levinson-Durbin recursions

step_down_LD_recursions <- function(coeffs, variance=1, process_variance_p=TRUE)
{
    p <- length(coeffs)
    if(p < 1) stop("phis (AR coefficients) must be not specified")
    pm1 <- p - 1
    k <- p
    temp <- vector(mode="list", length=p)
    temp[[p]] <- coeffs
    coeffs_k <- coeffs
    if(p>1)
    {
        for (i in 1:pm1)
        {
            km1 <- k - 1
            den <- 1 - coeffs_k[k]^2
            js <- 1:km1
            temp[[km1]] <-(coeffs_k[js] + coeffs_k[k]*coeffs_k[k-js])/den
            coeffs_k <- temp[[km1]]
            k <- k - 1
        }
    }
    pe_vars <- rep(0,p+1)
    if(process_variance_p)
    {
        pe_vars[1] <- variance
        for (i in 1:p)
            pe_vars[i+1] <- pe_vars[i] * (1 - (temp[[i]][i])^2)
    }
    else
    {
        pe_vars[p+1] <- variance
        for (i in p:1)
            pe_vars[i] <- pe_vars[i+1] / (1 - (temp[[i]][i])^2)
    } 
    return(list(coeffs=temp,pe_vars=pe_vars))
}
