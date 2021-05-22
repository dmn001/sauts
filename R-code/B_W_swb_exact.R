### smoothing window bandwith due to Jenkins

B_W_swb_exact <- function(lw,delta_t=1) 1/(delta_t*(1+ 2*sum(lw[-1]^2)))
