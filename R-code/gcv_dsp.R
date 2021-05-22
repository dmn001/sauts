##########################################################################
### implementation of bandwidth selection method described by
### Ombao et al. (2001) - see their Equation (9)
###
### NOTE: Ombao et al. used a rectangular smoother in their simulation
###       studies.  The default smoother here is Gaussian,
###       but setting Gaussian_p to FALSE yields rectangular
##########################################################################

gcv_dsp <- function(m,pgram_or_ts,pgram_p=TRUE,keep_0_p=TRUE,Gaussian_p=TRUE,swap_p=FALSE,delta_t=1)
{
    N <- length(pgram_or_ts)
    ## the_pgram should have values of the periodogram at Fourier
    ## frequencies ranging from zero to just below twice Nyquist ...
    the_pgram <- if(pgram_p) pgram_or_ts else pgram(ts,two_sided=TRUE,start=TRUE,delta_t=delta_t)$sdfe
    Nyquist_p <- is_even(N)
    ##
    M_chi_2 <- floor((N-1)/2)  # number of Fourier frequencies > 0 and < Nyquist
    M <- M_chi_2 + (if(keep_0_p) 1 else 0) + (if(Nyquist_p) 1 else 0)
    ##
    dft_pgram <- dft(the_pgram)
    temp <- if(Gaussian_p) dnorm(0:floor((N-1)/2),sd=1/(sqrt(2)*m))
            else rep(1,floor((N+1)/2))
    g_raw <- c(temp,if(Nyquist_p) 0 else NULL,rev(temp[-1]))
    g_m <- g_raw/sum(g_raw)
    edofs <- 2/sum(g_m^2)
    dft_g_m <- dft(g_m)
    B_U <- B_U_ds_pgram(g_m,delta_t=delta_t)
    dofs_model <- (1-g_m[1])^2
    ds_pgram <- Re(inverse_dft(dft_pgram*dft_g_m))
    temp <- if(swap_p) (-log(ds_pgram/the_pgram) + (ds_pgram-the_pgram)/the_pgram)/dofs_model else (-log(the_pgram/ds_pgram) + (the_pgram-ds_pgram)/ds_pgram)/dofs_model
    gcv <- (sum(temp[2:(M_chi_2+1)]) +
            if(keep_0_p) 0.5*temp[1] else 0 + if(Nyquist_p) 0.5*temp[N/2+1] else 0)/M
    return(list(N=N,
                Nyquist_p=Nyquist_p,
                M_chi_2=M_chi_2,
                M=M,
                the_pgram=the_pgram,
                g_m=g_m,
                edofs=edofs,
                B_U=B_U,
                dofs_model=dofs_model,
                cc=list(up = 10*log10(edofs/qchisq(0.025,edofs)),
                        down = -10*log10(edofs/qchisq(0.975,edofs)),
                        edof = edofs,
                        width = B_U),
                ds_pgram=ds_pgram,
                freqs=(0:(N-1))/(N*delta_t),
                sdfe=ds_pgram,
                gcv=gcv))
}
