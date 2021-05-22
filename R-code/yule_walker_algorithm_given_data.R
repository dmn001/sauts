### 

yule_walker_algorithm_given_data <- function(ts, p=1, center=TRUE)
{
    if(center) ts <- ts - mean(ts)
    e_F <- c(ts,rep(0,p))
    e_B <- c(0,ts,rep(0,p-1))
    N <- length(ts)
    pev  <- rep(sum(ts^2)/N,p+1)
    pacs <- rep(0,p)
    blpc <- vector(mode="list", length=p)
    fpe  <- vector(mode="list", length=p)
    bpe  <- vector(mode="list", length=p)
    for(k in 1:p)
    {
        phis <- rep(0,k)
        ## compute kth order pacs (reflection coefficient)
        phis[k] <- sum(e_F*e_B)/sum(e_F^2)
        if(k>1) phis[1:(k-1)] <- old_phis - phis[k]*rev(old_phis)
        blpc[[k]] <- phis
        pacs[k]  <- phis[k]
        pev[k+1] <- pev[k]*(1-phis[k]^2)
        old_phis <- phis
        old_e_F <- e_F
        old_e_B <- e_B
        e_F <- old_e_F - phis[k]*old_e_B
        e_B <- circular_shift(old_e_B - phis[k]*old_e_F,1)
        fpe[[k]] <- e_F[(k+1):N]
        bpe[[k]] <- e_B[(k+2):(N+1)]
    }
    return(list(coeffs=phis,
                innov_var=pev[p+1],
                pev=pev,
                pacs=pacs,
                blpc=blpc,
                forward_pe=fpe,
                backward_pe=bpe,
                N=N))
}

yule_walker <- yule_walker_algorithm_given_data
