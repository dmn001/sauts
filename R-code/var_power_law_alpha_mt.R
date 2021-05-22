var_power_law_alpha_mt <- function(N,K=5,f_H=1/3)
    {
        nonzero_Fourier_freqs <- (1:floor(N/2))/N
        M <- sum(nonzero_Fourier_freqs <= f_H)
        M_0 <- ceiling((K+1)/2)
        gg <- (nonzero_Fourier_freqs >= M_0/N) & (nonzero_Fourier_freqs <= f_H)
        xs <- log(nonzero_Fourier_freqs)
        xs_centered <- xs - mean(xs[gg])
        temp <- 0
        for(j in 0:(M-M_0)) temp <- temp + xs_centered[M_0+j]*sum(sapply(max(-j,-K):min(K,M-M_0-j),function(k) xs_centered[k+M_0+j]*(1-abs(k)/(K+1))))
        return(list(M=M,
                    M_0=M_0,
                    N_xs=sum(gg),
                    var_alpha=trigamma(K)*temp/(sum(xs_centered[gg]^2))^2))
    }

### for backward compatability (ZAP EVENTUALLY!)

mt_var_alpha <- var_power_law_alpha_mt
