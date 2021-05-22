### compute expected value of square of biased estimator of ACVS
### NOTE: even for N=64, it takes a LONG time to execute this function!!!

ev_shp_squared <- function(tau,N,the_acvs)
{
    a_tau <- abs(tau)
    s_func <- function(tau) return(the_acvs[abs(tau)+1])
    A_func <- function(j,k,l,m)
    {
        s_func(j-k)*s_func(l-m)+s_func(j-l)*s_func(k-m)+s_func(j-m)*s_func(k-l)
    }
    B_func <- function(j,k,l)
    {
        (s_func(j-k)*sum(sapply(0:(N-1),function(u) s_func(l-u)))+
            s_func(j-l)*sum(sapply(0:(N-1),function(u) s_func(k-u)))+
            s_func(k-l)*sum(sapply(0:(N-1),function(u) s_func(j-u))))/N
    }
    C_func <- function(j,k)
    {
        (s_func(j-k)*sum(sapply(0:(N-1),function(u) sum(sapply(0:(N-1),function(v) s_func(u-v)))))+
            2*sum(sapply(0:(N-1),function(u) s_func(j-u)))*sum(sapply(0:(N-1),function(u) s_func(k-u))))/N^2
    }
    D_func <- function(j)
    {
        3*sum(sapply(0:(N-1),function(u) s_func(j-u)))*sum(sapply(0:(N-1),function(u) sum(sapply(0:(N-1),function(v) s_func(u-v)))))/N^3
    }
    E_func <- function()
    {
        3*(sum(sapply(0:(N-1),function(u) sum(sapply(0:(N-1),function(v) s_func(u-v))))))^2/N^4
    }
    double_sapply <- sum(sapply(0:(N-a_tau-1),function(t) sum(sapply(0:(N-a_tau-1),function(u) A_func(t+a_tau,t,u+a_tau,u)-B_func(t+a_tau,t,u+a_tau)-B_func(t+a_tau,t,u)-B_func(t+a_tau,u+a_tau,u)-B_func(t,u+a_tau,u)+C_func(t+a_tau,u+a_tau)+C_func(t+a_tau,u)+C_func(t,u+a_tau)+C_func(t,u)))))
    single_sapply <- 2*(N-a_tau)*sum(sapply(0:(N-a_tau-1),function(t) C_func(t+a_tau,t)-D_func(t+a_tau)-D_func(t)))
    no_sapply <- (N-a_tau)^2*E_func()
    return((double_sapply+single_sapply+no_sapply)/N^2)
}
