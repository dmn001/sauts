### compute expected value of biased estimator of ACVS

ev_shp <- function(tau,N,the_acvs)
{
    a_tau <- abs(tau)
    s_func <- function(tau) return(the_acvs[abs(tau)+1])
    double_sums <- sum(sapply(0:(N-a_tau-1),function(s) sum(sapply(0:(N-1),function(u) s_func(s-u+a_tau)+s_func(s-u)))))
    single_sum <- s_func(0) + 2*sum(sapply(1:(N-1),function(u) (1-u/N)*s_func(u)))
    return((N-a_tau)*s_func(tau)/N - double_sums/N^2 + (N-a_tau)*single_sum/N^2)
}
