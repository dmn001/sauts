###

gen_cis_for_ar_sdf <- function(N_ts,ar_coeffs,innov_var,N_pad=N_ts,delta_t=1,alpha=0.05,zap_me=if(is_odd(N_pad)) 1 else c(1,N_pad/2+1))
{
    p <- length(ar_coeffs)
    sdf_unit_iv <- ar_coeffs_to_sdf(ar_coeffs, innov_var=1, delta_t=delta_t, N_pad=N_pad)
    G_hat_abs <- sqrt(delta_t/sdf_unit_iv$sdfe[-zap_me])
    G_hat_freqs <- sdf_unit_iv$freqs[-zap_me]
    A <- matrix(0,nrow=p,ncol=p)
    B <- matrix(0,nrow=p,ncol=p)
    diag(A) <- 1
    diag(B) <- ar_coeffs[p]
    for(j in 1:(p-1))
    {
        A_sub_diag <- -ar_coeffs[j]
        B_sub_diag <-  ar_coeffs[p-j]
        for(k in 1:(p-j))
        {
            A[j+k,k] <- A_sub_diag
            B[j+k,k] <- B_sub_diag
        }
    }
    sig_2_Gamma_inv <- t(A) %*% A - t(B) %*% B
    e_vecs <- matrix(nrow=p,ncol=length(G_hat_freqs))
    one_to_p <- 1:p
    for(k in 1:length(G_hat_freqs))
        e_vecs[,k] <- exp(complex(imag=2*pi*G_hat_freqs[k]*one_to_p*delta_t))
    beta <- 1 - sqrt(1-alpha)
    percent_point <- qchisq(p=1-beta,df=2)
    r_hat_2 <- sapply(1:length(G_hat_freqs),function(k) sqrt(Re(percent_point*t(Conj(e_vecs[,k])) %*% sig_2_Gamma_inv %*% e_vecs[,k]/N_ts)))
    lower <- 1/(G_hat_abs + r_hat_2)^2
    upper <- 1/(G_hat_abs - r_hat_2)^2
    return(list(freqs=G_hat_freqs,
                G_hat_abs=G_hat_abs,
                r_hat_2=r_hat_2,
                sdfe=innov_var*delta_t/G_hat_abs^2,
                lower=N_ts*innov_var*delta_t*lower/qchisq(p=1-beta/2,df=N_ts-p),
                upper=N_ts*innov_var*delta_t*upper/qchisq(p=beta/2,df=N_ts-p)))
}
