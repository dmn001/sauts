###

bartlett_priestley_design_window <- function(m,N_freq=1025,delta_t=1)
{
    freqs <- (0:(N_freq-1))/(2*N_freq)
    dw <- sapply(freqs,function(f) if(abs(f) > 1/(2*m)) 0 else 1.5*m*(1-(2*f*m)^2)) 
    list(freqs=freqs/delta_t,
         dw=dw*delta_t)
}

quadratic_design_window    <- bartlett_priestley_design_window
epanechnikov_design_window <- bartlett_priestley_design_window

