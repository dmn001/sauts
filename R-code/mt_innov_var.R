###

mt_innov_var <- function(ts,the_tapers,center=TRUE)
{
    if(center) ts <- ts - mean(ts)
    N <- length(ts)
    M_mt <- floor(N/2)
    sdfe_mt <- mt_sdf_est(ts,the_tapers=the_tapers,center=center)$sdfe[2:(M_mt+1)]
    K <- ncol(the_tapers)
    digammas <- rep(digamma(K),M_mt)
    if(M_mt == N/2) digammas[M_mt] <- digamma(K/2)
    return(K*exp(mean(log(sdfe_mt) - digammas - trigamma(K)/(2*M_mt))))
}
