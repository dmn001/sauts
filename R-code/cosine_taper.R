### cosine data taper

cosine_taper <- function(N,p=1.0)
{
    if(is.null(p)) p <- 1.0
    taper <- rep(1,N)
    M <- floor(p*N/2)
    if(M > 0)
    {
        range <- 1:M
        taper[range] <- 0.5 * (1 - cos(pi*range/(M + 1)))
        taper[(N+1-M):N] <- taper[M:1]
    }
    return(taper/sqrt(sum(taper^2)))
}
