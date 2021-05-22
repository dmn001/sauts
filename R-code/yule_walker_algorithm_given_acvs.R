### Yule-Walker method

yule_walker_algorithm_given_acvs <- function(acvs, p=length(acvs)-1)
{
    blpc <- vector(mode="list", length=p)
    phis <- acvs[2]/acvs[1]
    pev  <- rep(acvs[1],p+1)
    blpc[[1]] <- phis
    pacs <- rep(phis,p)
    pev[2] <- pev[1]*(1-phis^2)
    if(p > 1)
    {
        for(k in 2:p)
        {
            old_phis <- phis
            phis <- rep(0,k)
            ## compute kth order pacs (reflection coefficient)
            phis[k] <- (acvs[k+1] - sum(old_phis*acvs[k:2]))/pev[k]
            phis[1:(k-1)] <- old_phis - phis[k]*rev(old_phis)
            blpc[[k]] <- phis
            pacs[k]  <- phis[k]
            pev[k+1] <- pev[k]*(1-phis[k]^2)
        }
    }
    return(list(coeffs=phis,
                innov_var=pev[p+1],
                pev=pev,pacs=pacs,
                blpc=blpc))
}
