###

compute_cD_and_Im <- function(C_D,N_prime,max_m,lag_window=parzen_lag_window,delta_t=1)
{
    N_L <- trunc((N_prime-1)/2) + 1
    N_U <- N_prime - N_L + 1
    ## NOTE: C_D should have length N_U (SHOULD CHECK THIS!!!)
    C_D_two_sided <- c(C_D,rev(C_D[2:N_L]))
    c_D <- Re(inverse_dft(C_D_two_sided))/delta_t
    lambda <- dB(exp(1))
    gamma <- -digamma(1)
    for_first_sum <- c_D^2 - (lambda*pi)^2/(6*N_prime*delta_t^2)
    taus <- c(0:(N_U-1),(N_L-1):1)
    C_1 <- (lambda*gamma)^2/delta_t
    C_2 <- (lambda*pi)^2/(6*N_prime*delta_t)
    even_p <- (N_prime%%2) == 0
    get_Im <- function(m)
    {
        p_lw_m <- sapply(taus,lag_window,m)
        C_1 + delta_t*sum(for_first_sum*(1-p_lw_m)^2) + C_2*(sum(p_lw_m^2) + 1 + if(even_p) p_lw_m[N_U]^2 else 0)
    }
    return(list(I_m=sapply(1:max_m,get_Im),
                c_D=c_D))
}

