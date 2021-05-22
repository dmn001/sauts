### smoothing window bandwidths ...

beta_W_swb_exact <- function(lw,delta_t=1)
  {
    N_m_1 <- length(lw) - 1
    taus_odd <- rev(seq(1,N_m_1,2))
    taus_even <- rev(seq(2,N_m_1,2))
    sqrt(1 + (12/pi^2) * (sum(lw[taus_even+1]/taus_even^2) - sum(lw[taus_odd+1]/taus_odd^2)))/delta_t
  }

### useful in cases where beta is imaginary ...

beta_W_swb_exact_squared <- function(lw,delta_t=1)
  {
    N_m_1 <- length(lw) - 1
    taus_odd <- rev(seq(1,N_m_1,2))
    taus_even <- rev(seq(2,N_m_1,2))
    (1 + (12/pi^2) * (sum(lw[taus_even+1]/taus_even^2) - sum(lw[taus_odd+1]/taus_odd^2)))/delta_t
  }
