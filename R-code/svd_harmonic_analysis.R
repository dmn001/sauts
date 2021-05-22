###

svd_harmonic_analysis <- function(ts, p=2, p_prime=25)
{
    N <- length(ts)
    Z_vec <- c(ts[(p_prime+1):N],Conj(ts[1:(N-p_prime)]))
    A_mat <- rbind(embed(ts[-N],p_prime),Conj(embed(rev(ts[-1]),p_prime)[(N-p_prime):1,]))
    svd_A_mat <- svd(A_mat)    
    U_mat <- svd_A_mat$u
    V_mat <- svd_A_mat$v
    Lambda_mat <- diag(svd_A_mat$d)
    ## followin works for special case of p = 2
    ## varphi_hat_prime <- as.vector((V_mat[,1] %o% Conj(U_mat[,1])/svd_A_mat$d[1] + V_mat[,2] %o% Conj(U_mat[,2])/svd_A_mat$d[2]) %*% Z_vec)
    varphi_hat_prime <- as.vector(Reduce(`+`,lapply(1:p,function(k){V_mat[,k] %o% Conj(U_mat[,k])/svd_A_mat$d[k]})) %*% Z_vec)
    return(list(roots=polyroot(rev(c(1,-varphi_hat_prime))),
                eigenvalues=svd_A_mat$d))
}
