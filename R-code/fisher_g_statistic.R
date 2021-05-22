### 

fisher_g_statistic <- function(ts, exact_p=FALSE, ...)
{
    N <- length(ts)
    M <- floor((N-1)/2)
    if(M <= 1)
    {
        g_statistic <- 1
        g_F <- 1
        reject_null_hypothesis <- FALSE 
    }
    else
    {
        pgram_no_0_no_Nyquist <- pgram(ts)$sdfe[2:(M+1)]
        g_statistic <- max(pgram_no_0_no_Nyquist)/sum(pgram_no_0_no_Nyquist)
        g_F <- if(exact_p) g_F_exact_via_M(M,...)
               else
               {
                   if(N>3000) g_F_hat(N,...) else g_F_tilde(M,...)
               }
        reject_null_hypothesis <- g_statistic > g_F
    }
    return(list(g_statistic=g_statistic, g_F=g_F, reject_null_hypothesis=reject_null_hypothesis))
}

###

g_F_tilde <- function(N,alpha=0.05)
{
    M <- floor((N-1)/2)
    return(1 - (alpha/M)^(1/(M-1)))
}

###

g_F_tilde_via_M <- function(M,alpha=0.05) return(1 - (alpha/M)^(1/(M-1)))

###

g_F_hat <- function(N,alpha=0.05) -2*log(-2*log(1-alpha)/N)/N

###

g_F_exact_via_M <- function(M,alpha=0.05,low_fiddle=1.1,high_fiddle=0.995,my_tol=.Machine$double.eps^0.5,table=NA)
{
    if(is.na(table) || (!is.na(table) && !is.element(M,table[,1])))
    {
            temp <- g_F_tilde_via_M(M,alpha=alpha)
            return(uniroot(function(x) g_F_distribution_via_M(x,M)-alpha,c(low_fiddle*temp,high_fiddle*temp),tol=my_tol)$root)
        }
    else  # if table is supplied, get g_F from it rather than by a (potentially lengthy) computation
        return(as.vector(table[which(M == table[,1]),2]))
}

###

g_F_distribution_via_N <- function(g_F,N)
{
    M <- floor((N-1)/2)
    if(g_F >= 1) return(0)
    else
    {
        if(g_F <= 1/M) return(1)
        else
            return(sum(sapply(1:min(floor(1/g_F),M),function(j) (-1)^(j-1) * exp(lgamma(M+1)-lgamma(M-j+1)-lgamma(j+1)) * (1-j*g_F)^(M-1))))
    }
}

###

g_F_distribution_via_M <- function(g_F,M)
{
    if(g_F >= 1) return(0)
    else
    {
        if(g_F <= 1/M) return(1)
        else
            return(sum(sapply(1:min(floor(1/g_F),M),function(j) (-1)^(j-1) * exp(lgamma(M+1)-lgamma(M-j+1)-lgamma(j+1)) * (1-j*g_F)^(M-1))))
    }
}
