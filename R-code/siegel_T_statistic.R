###

siegel_c_and_beta <- function(M,lam=0.6,...)
{
    g_F <- g_F_exact_via_M(M,...)
    E_T <- (1 - lam*g_F)^M 
    E_T_squared <- (2*(1 - lam*g_F)^(M+1) + (M-1)*max(0,(1 - 2*lam*g_F))^(M+1))/(M+1)
    my_c <- (E_T_squared - E_T^2)/(4*E_T)
    list(c=my_c, beta=E_T/my_c, E_T=E_T, E_T_squared=E_T_squared)
}

###

siegel_t_distribution_via_M <- function(x,M,lam=0.6,...)
{
    g_F <- g_F_exact_via_M(M,...)
    Q_j_k <- function(j,k) lgamma(M+1) - lgamma(M-j+1) - lgamma(j+1) +
                               lgamma(j) - lgamma(j-k) - lgamma(k+1) +
                               lgamma(M) - lgamma(M-k) - lgamma(k+1) +
                               k*log(x) + (M-k-1)*log(max(1-j*lam*g_F-x,0))
    return(x^(M-1) + sum(sapply(1:M,function(j) sum(sapply(0:min(j-1,M-2),function(k) (-1)^(j+k+1) * exp(Q_j_k(j,k)))))))
}

###

siegel_t_exact_via_M <- function(M,interval,lam=0.6,alpha=0.05,low_fiddle=1,high_fiddle=0.995,my_tol=.Machine$double.eps^0.5)
{
    return(uniroot(function(x) siegel_t_distribution_via_M(x,M,lam=lam,alpha=alpha,low_fiddle=low_fiddle,high_fiddle=high_fiddle,my_tol=my_tol)-alpha,interval,tol=my_tol)$root)
}

###

siegel_t_approximate_lam_0p6_via_M <- function(M,alpha=0.05)
{
    return(if(alpha==0.05) 1.033 * M^(-0.72356) else if(alpha==0.01) 1.4987 * M^(-0.79695) else stop(paste("alpha (", alpha, ") must be set to 0.1 or 0.05\n",sep="")))
}

