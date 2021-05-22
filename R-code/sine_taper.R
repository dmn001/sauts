### kth-order sine data 

sine_taper <- function(N,k=0)
{
    if(is.null(k)) k <- 0
    return(sqrt(2/(N+1))*sin((k+1)*pi*(1:N)/(N+1)))
}
