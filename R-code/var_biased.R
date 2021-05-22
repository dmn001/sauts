var_biased <- function(x)
{
    N <- length(x)
    return((N-1)*var(x)/N)
}
