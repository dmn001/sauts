### zeroth-order Slepian data taper as approximated by Equation (196b)

slepian_taper <- function(N,NW=4)
{
    if(is.null(NW)) NW <- 4
    W <- NW/N
    temp <- pi*W*(N-1)
    taper <- besselI(temp*sqrt(1-((2*(0:(N-1)) + 1 - N)/N)^2),0)
    C <- sqrt(sum(taper^2))
    return(taper/C)
}
