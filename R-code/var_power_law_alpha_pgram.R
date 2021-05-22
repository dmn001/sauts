var_power_law_alpha_pgram <- function(N,f_H=1/3)
{
    nonzero_Fourier_freqs <- (1:floor(N/2))/N
    xs <- log(nonzero_Fourier_freqs[which(nonzero_Fourier_freqs <= f_H)])
    x_bar <- mean(xs)
    return(pi^2/(6*sum((xs - x_bar)^2)))
}

### for backward compatability (ZAP EVENTUALLY!)

var_alpha_pgram <- var_power_law_alpha_pgram
