### multitaper-based CIs for SDF based on two types of jackknifing

thomson_chave_miller_jackknife <- function(mt_est)
{
    N_freq <- length(mt_est$freqs)
    lower_CI_tc <- rep(0,N_freq)
    upper_CI_tc <- rep(0,N_freq)
    lower_CI_miller <- rep(0,N_freq)
    upper_CI_miller <- rep(0,N_freq)
    Stmt <- rep(0,N_freq)
    K <- length(mt_est$eigenspectra)
    everyone <- 1:K
    t_lower <- qt(0.025,K-1)  # same as -qt(0.975,K-1)
    t_upper <- qt(0.975,K-1)  # same as -qt(0.025,K-1)
    for(k in 1:N_freq)
    {
        ## Thomson-Chave jackknife
        Shmt <- mt_est$sdfe[k]
        log_Shmt <- log(Shmt)
        Shmtk <- unlist(lapply(mt_est$eigenspectra,function(x) x$sdfe[k]))
        log_Shmt_minus_j <- sapply(1:K, function(m) log(sum(Shmtk[everyone[-m]])/(K-1)))
        log_Shmt_minus <- mean(log_Shmt_minus_j)
        s_hat <- sqrt( (K-1)*sum((log_Shmt_minus_j - log_Shmt_minus)^2)/K )
        lower_CI_tc[k] <- exp(s_hat*t_lower) * Shmt
        upper_CI_tc[k] <- exp(s_hat*t_upper) * Shmt
        ## Miller jackknife
        log_Stmt <- K*log_Shmt - (K-1)*log_Shmt_minus
        Stmt[k] <- exp(log_Stmt)
        lower_CI_miller[k] <- exp(s_hat*t_lower) * Stmt[k]
        upper_CI_miller[k] <- exp(s_hat*t_upper) * Stmt[k]
    }
    return(list(freqs=mt_est$freqs,
                sdfe=mt_est$sdfe,
                debiased_sdfe=Stmt,
                tc_lci=lower_CI_tc,
                tc_uci=upper_CI_tc,
                m_lci=lower_CI_miller,
                m_uci=upper_CI_miller,
                tc_dB_range=range(dB(upper_CI_tc/lower_CI_tc)),
                m_dB_range=range(dB(upper_CI_miller/lower_CI_miller)),
                tc_dB_ave=sum(dB(upper_CI_tc/lower_CI_tc))/N_freq,
                m_dB_ave=sum(dB(upper_CI_miller/lower_CI_miller))/N_freq,
                tc_dB_median=median(dB(upper_CI_tc)) - median(dB(lower_CI_tc)),
                m_dB_median=median(dB(upper_CI_miller)) - median(dB(lower_CI_miller))))
}
