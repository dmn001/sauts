### 

cumulative_pgram <- function(ts, alpha=0.05, delta_t=1)
{
    N <- length(ts)
    M <- round(if(is_odd(N)) (N-1)/2 else (N-2)/2)
    freqs <- (1:(M-1))/(N*delta_t)
    if(M < 2) stop(paste("length of ts (N=",N,") is too small: N must be greater than 4",sep=""))
    else
    {
        pgram_no_0_no_Nyquist <- pgram(ts,center=FALSE)$sdfe[2:(M+1)]
        c_pgram <- cumsum(pgram_no_0_no_Nyquist[-M])/sum(pgram_no_0_no_Nyquist)
        ##
        temp     <- (1:(M-1))/(M-1) - c_pgram
        i_temp   <- which.max(temp)
        D_plus   <- temp[i_temp]
        f_plus   <- freqs[i_temp]
        ##
        temp     <- c_pgram - (0:(M-2))/(M-1)
        i_temp   <- which.max(temp)
        D_minus  <- temp[i_temp]
        f_minus  <- freqs[i_temp]
        ##
        temp     <- c(D_plus, D_minus)
        i_temp   <- which.max(temp)
        D_stat   <- temp[i_temp]
        f_D_stat <- c(f_plus,f_minus)[i_temp]
        ##
        tilde_D_alpha <- tilde_D(M,alpha)
        reject_null_hypothesis <- D_stat > tilde_D_alpha
        L_u <- cbind(c(0,((M-1)*(1-tilde_D_alpha)+1)/(N*delta_t)),c(tilde_D_alpha-1/(M-1),1))
        L_l <- cbind(c((M-1)*tilde_D_alpha/(N*delta_t),1/(2*delta_t)),c(0,N*0.5/(M-1)-tilde_D_alpha))
        return(list(D_stat=D_stat, tilde_D_alpha=tilde_D_alpha, reject_null_hypothesis=reject_null_hypothesis, alpha=alpha, c_pgram=c_pgram, freqs=freqs, L_u=L_u, L_l=L_l, D_plus=D_plus, f_plus=f_plus, D_minus=D_minus, f_minus=f_minus, f_D_stat=f_D_stat))
    }
}

### values from p. 732 of Stephens (1974), second line of Table 1A
### percentages (times 100):  0.15  0.1   0.05  0.025 0.01
### percentage points:        1.138 1.224 1.358 1.480 1.628

tilde_D <- function(M,alpha=0.05)
    {
        C_alpha <-  if(isTRUE(all.equal(alpha,0.15))) 1.138 else if(isTRUE(all.equal(alpha,0.1))) 1.224 else if(isTRUE(all.equal(alpha,0.05))) 1.358 else if(isTRUE(all.equal(alpha,0.025))) 1.480 else if(isTRUE(all.equal(alpha,0.01))) 1.628 else stop(paste("alpha (", alpha, ") must be set to 0.15, 0.1, 0.05, 0.025 or 0.01\n",sep=""))
        temp <- sqrt(M-1)
        return(C_alpha/(temp + 0.12 +0.11/temp))
    }
