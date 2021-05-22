### compute periodogram at one or more particular frequencies

pgram_single_freqs <- function(f,ts,delta_t=1,center=TRUE)
{
    N <- length(ts)
    if(center) ts <- ts-mean(ts)
    N_f <- length(f)
    temp <- rep(0,N_f)
    my_sum <- if(N_f == 1) sum else rowSums
    for(n in 1:N) temp <- temp + my_sum(sapply(1:N,function(m) ts[n]*ts[m]*cos(2*pi*f*(m-n)*delta_t)))
    pgram <- delta_t*temp/N
    return(list(freqs=f, sdfe=pgram, cc=list(up=10*log10(2/qchisq(0.025,2)),down=-10*log10(2/qchisq(0.975,2)),edof=2,width=3*N/((2*N^2+1)*delta_t))))
}

### deprecated version that works with just one frequency
### 
### pgram_single_freq <- function(f,ts,delta_t=1,center=TRUE)
### {
###     N <- length(ts)
###     if(center) ts <- ts-mean(ts)
###     temp <- 0
###     for(n in 1:N) temp <- temp + sum(sapply(1:N,function(m) ts[n]*ts[m]*cos(2*pi*f*(m-n)*delta_t)))
###     pgram <- delta_t*temp/N
###     return(list(freqs=f, sdfe=pgram, cc=list(up=10*log10(2/qchisq(0.025,2)),down=-10*log10(2/qchisq(0.975,2)),edof=2,width=3*N/((2*N^2+1)*delta_t))))
### }

