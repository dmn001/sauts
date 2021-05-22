### function for creating the Monte Carlo study that is the focus
### of Figure 406

do_it_innov_var <- function(ts,Ls=c(2,3,4),Ks=c(2,3,4),the_tapers=NULL)
{
    N <- length(ts)
    ## Davis-Jones _..
    M_dj <- floor((N-1)/2)
    euler <- -digamma(1)
    the_pgram <- pgram(ts,center=FALSE)$sdfe[2:(M_dj+1)]
    iv_dj <- exp(mean(log(the_pgram)+euler))
    ## Hannan-Nicholls ...
    N_Ls <- length(Ls)
    iv_hn <- rep(NA,N_Ls)
    for(l in 1:N_Ls)
    {
        L <- Ls[l]
        M_hn <- floor((N-1)/(2*L))
        iv_hn[l] <- L*exp(mean(log(colMeans(matrix(the_pgram[1:(L*M_hn)],nrow=L)))-digamma(L)))
    }
    ## Pukkilla-Nyquist ...
    M_pn <- M_dj
    the_PN_taper_1 <- cosine_taper(N,12.5/N)
    the_PN_taper_2 <- cosine_taper(N,13/N)
    the_PN_taper_3 <- trapezoidal_taper(N)
    iv_pn <- rep(NA,3)
    iv_pn[1] <- exp(mean(log(direct_sdf_est(ts,the_taper=the_PN_taper_1,center=FALSE)$sdfe[2:(M_pn+1)]) + euler - pi^2/(12*M_pn)))
    iv_pn[2] <- exp(mean(log(direct_sdf_est(ts,the_taper=the_PN_taper_2,center=FALSE)$sdfe[2:(M_pn+1)]) + euler - pi^2/(12*M_pn)))
    iv_pn[3] <- exp(mean(log(direct_sdf_est(ts,the_taper=the_PN_taper_3,center=FALSE)$sdfe[2:(M_pn+1)]) + euler - pi^2/(12*M_pn)))
    ## multitaper ...
    M_mt <- floor(N/2)
    N_Ks <- length(Ks)
    iv_mt <- rep(NA,N_Ks)
    if(is.null(the_tapers)) the_tapers <- t(as.matrix(taper("sine",N,n.taper=max(Ks))))
    for(k in 1:N_Ks)
    {
        K <- Ks[k]
        digammas <- rep(digamma(K),M_mt)
        if(M_mt == N/2) digammas[M_mt] <- digamma(K/2)
        iv_mt[k] <- K*exp(mean(log(mt_sdf_est(ts,the_tapers=the_tapers[,1:K],center=FALSE)$sdfe[2:(M_mt+1)]) - digammas - trigamma(K)/(2*M_mt)))
    }
    return(list(dj=iv_dj,hn=iv_hn,pn=iv_pn,mt=iv_mt))
}

