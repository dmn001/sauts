### Hanning data taper
###
### NOTE: requires is_even

hanning_taper <- function(N,parm=NULL)
{
    taper <- rep(0,N)
    range <- 1:floor((N+1)/2)
    taper[range] <- 1 - cos(pi*range/(floor(N/2)+1))
    taper[-range] <- rev(taper)[-range]
    return(taper * if(is_even(N)) sqrt(2/(3*N-2)) else sqrt(2/(3*N+3)))
}
