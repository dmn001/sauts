### compute spectral window corresponding to direct spectral estimate

spec_window <- function(the_taper,pad_factor=4,delta_t=1,two_sided_p=TRUE,fix_nulls_p=FALSE,...)
{
    N <- length(the_taper)
    N_padded <- round(N * pad_factor)
    taper_padded <- c(the_taper,rep(0,N_padded-N))
    M <- floor(N_padded/2)  # number of nonzero frequencies
    freqs <- (0:M)/(N_padded*delta_t)
    sw <- delta_t*abs(fft(taper_padded)[1:(M+1)])^2
    sw <- if(two_sided_p) c(rev(sw[-1]),sw) else sw
    if(fix_nulls_p)
    {
        sw_null_fixer <- function(sw,dB_level=-120.0,first_p=TRUE,last_p=TRUE,tol=5)
        {
            eps <- 10^(dB_level/10)
            tol <- 10^(tol/10)  # not needed???
            N <- length(sw)
            r <- 2:(N-1)
            patch_me <- (sw[r-1] > sw[r]) & (sw[r+1] > sw[r]) # & (sw[r-1] > (sw[r]+tol))
            sw[c(first_p > eps, patch_me, last_p)] <- eps
            return(sw)
        }
        sw <- sw_null_fixer(sw,...)
    }
    return(list(freqs=if(two_sided_p) c(-rev(freqs[-1]),freqs) else freqs,
                sw=sw))
}

###

spec_window_central_lobe <- function(the_sw,freq_center=0,sw_mult=1,dB_cutoff=-20)
{
    freqs_all <- the_sw$freqs
    sw_all <- the_sw$sw
    i_freq_0 <- which(freqs_all == 0)  # should only be frequency equal to 0
    test_me <- (i_freq_0+1):length(freqs_all)
    gg <- c(TRUE,sw_all[test_me-1] > sw_all[test_me])
    i_adjust <- match(FALSE,gg)
    if(!is.na(i_adjust)) gg[i_adjust:length(gg)] <- FALSE
    if(i_freq_0 != 1)    gg <- c(rev(gg)[-1],gg)
    freqs_central_lobe_all <- freqs_all[gg]
    sw_central_lobe_all <- sw_all[gg]
    gg <- sw_central_lobe_all >= 10^(dB_cutoff/10)
    return(list(freqs=freqs_central_lobe_all[gg]+freq_center,
                sw=sw_central_lobe_all[gg] * sw_mult))
}
