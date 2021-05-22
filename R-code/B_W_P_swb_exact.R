### smoothing window bandwidth due to Parzen

B_W_P_swb_exact <- function(lw,delta_t=1)
  {
    1/(delta_t*(1+ 2*sum(lw[-1])))
  }
