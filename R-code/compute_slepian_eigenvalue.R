### 

compute_slepian_eigenvalue <- function(eigenvector,NW=4)
{
    N <- length(eigenvector)
    W <- NW/N
    taus <- 1:(N-1)
    acvs_bl_wn <- sin(2*pi*W*taus)/(pi*taus)
    acs_eigenvector <- acvs(eigenvector,center=FALSE,acs=TRUE)$acvs[-1]
    return(2*(W + sum(acvs_bl_wn*acs_eigenvector)))
}
