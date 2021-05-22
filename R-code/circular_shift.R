### circularly shifts elements to right (if shift>0) or left (if shift<0)

circular_shift <- function(x,shift=0)
{
    N <- length(x)
    K <- shift %% N
    if(K == 0) x else c(x[(N-K+1):N],x[1:(N-K)])
}
