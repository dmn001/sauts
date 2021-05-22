### NOTE: multitapers are in the columns of the_tapers - each column should
###       have a sum of squares equal to unity

amt_sdf_est <- function(ts,the_tapers,NW,center=TRUE,recenter=FALSE,normalize=FALSE,pad_factor=1,delta_t=1,max_it=100,tol=0.05,percent_CI=0.95)
{
    N <- length(ts)
    K <- ncol(the_tapers)  # K is assumed to be 2 or more
    eigenspectra <- vector("list",K)
    N_padded <- round(N * pad_factor)
    M <- floor(N_padded/2) # number of nonzero frequencies
    freqs <- (0:M)/(N_padded*delta_t)
    n_its <- rep(NA,M+1)
    edofs <- rep(0,M+1)
    mu <- if(center) mean(ts) else 0
    sig2_delta_t <- delta_t * sum((ts - mu)^2)/N
    ## compute K eigenspectra
    for(k in 1:K)
    {
        a_taper <- the_tapers[,k]
        tapered_ts <- create_tapered_series(ts,a_taper,center=center,recenter=recenter,normalize=normalize)
        ts_padded <- c(tapered_ts,rep(0,N_padded-N))
        dse <- delta_t*abs(fft(ts_padded)[1:(M+1)])^2
        eigenspectra[[k]] <- structure(list(freqs=freqs, sdfe=dse, cc=do_crisscross_dse(a_taper,delta_t=delta_t)))
    }
    ## compute eigenvalues
    eigenvalues <- sapply(1:K,function(k) compute_slepian_eigenvalue(as.vector(the_tapers[,k]),NW))
    ## initial SDF estimate
    mt_est_initial <- (eigenvalues[1]*eigenspectra[[1]]$sdfe + eigenvalues[2]*eigenspectra[[2]]$sdfe)/sum(eigenvalues[1:2])
    mt_est_final <- mt_est_initial
    ## note: need to figure out how to get around these `for' loops!!!
    for(j_freq in 1:(M+1))  # frequency-by-frequency estimate ...
    {
        mt_est_current <- mt_est_initial[j_freq]
        for(i in 1:max_it) 
        {
            b_squared <- (mt_est_current/(eigenvalues*mt_est_current + (1-eigenvalues)*sig2_delta_t))^2
            d_num <- b_squared * eigenvalues
            d <- d_num/sum(d_num)
            edofs[j_freq] <- 2*(sum(b_squared*eigenvalues)^2)/sum((b_squared*eigenvalues)^2)
            mt_est_final[j_freq] <- sum(d*sapply(1:K,function(k) eigenspectra[[k]]$sdfe[j_freq]))
            ## do enough iterations so that mt_est_current is not the same as mt_est_initial
            if(i > 1 && (abs(mt_est_final[j_freq]-mt_est_current)/mt_est_current) < tol)
            {
                n_its[j_freq] <- i
                break()
            }
            mt_est_current <- mt_est_final[j_freq]
        }
    }
    p <- (1-percent_CI)/2
    return(list(freqs=freqs,sdfe=mt_est_final,edofs=edofs,ci_lower=edofs*mt_est_final/qchisq(1-p,edofs),ci_upper=edofs*mt_est_final/qchisq(p,edofs),eigenspectra=eigenspectra,p=p,percent_CI=percent_CI,sig2=sig2_delta_t/delta_t,eigenvalues=eigenvalues,n_its=n_its))
}
