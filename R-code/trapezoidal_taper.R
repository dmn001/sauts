### trapezoidal data taper (see Equation (8), Pukk79a)

trapezoidal_taper <- function(N,a=6.25/N)
{
    if(is.null(a)) a <- 6.25/N
    taper <- rep(1,N)
    M <- floor(a*N)
    if(M > 0)
    {
        range <- 1:M
        taper[range] <- range/(N*a)
        taper[(N+1-M):N] <- taper[M:1]
    }
    return(taper/sqrt(sum(taper^2)))
}
