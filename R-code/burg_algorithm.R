### Burg's algorithm

burg_algorithm <- function(ts, p=1, center=TRUE)
{
    if(center) ts <- ts - mean(ts)
    e_F <- ts
    e_B <- ts
    N <- length(ts)
    pev  <- rep(sum(ts^2)/N,p+1)
    pacs <- rep(0,p)
    blpc <- vector(mode="list", length=p)
    fpe  <- vector(mode="list", length=p)
    bpe  <- vector(mode="list", length=p)
    M <- N
    for(k in 1:p)
    {
        phis <- rep(0,k)
        ## compute kth order pacs (reflection coefficient)
        phis[k] <- 2*sum(e_F[-1]*e_B[-M])/sum(e_F[-1]^2 + e_B[-M]^2)
        if(k>1) phis[1:(k-1)] <- old_phis - phis[k]*rev(old_phis)
        blpc[[k]] <- phis
        pacs[k]  <- phis[k]
        pev[k+1] <- pev[k]*(1-phis[k]^2)
        old_phis <- phis
        old_e_F <- e_F
        old_e_B <- e_B
        e_F <- old_e_F[-1] - phis[k]*old_e_B[-M]
        e_B <- old_e_B[-M] - phis[k]*old_e_F[-1]
        fpe[[k]] <- e_F
        bpe[[k]] <- e_B
        M <- M - 1
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

burg <- burg_algorithm
