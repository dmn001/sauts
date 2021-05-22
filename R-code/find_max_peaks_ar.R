###

find_max_peaks_ar <- function(tss,p,ar_method=ar_fbls,N_peaks=1,pad_factor=512,N_pad=NULL,center=TRUE,delta_t=1,J=10)
{
    g_prime <- function(omega) - sum(the_acvs*taus*sin(omega*taus))
    g_prime_prime <- function(omega) - sum(the_acvs*taus^2*cos(omega*taus))
    tss <- as.matrix(tss)
    N_time_series <- ncol(tss)
    N      <- nrow(tss)
    peaks_ar <- matrix(nrow=N_peaks,ncol=N_time_series)
    taus <- 0:p
    for(n in 1:N_time_series)
    {
        a_ts <- tss[,n]
        the_ar_coeffs <- ar_method(a_ts,p,center=center)$coeffs
        the_acvs <- acvs(c(-1,the_ar_coeffs),center=FALSE)$acvs
        temp <- ar_coeffs_to_sdf(the_ar_coeffs,pad_factor=pad_factor,N_pad=N_pad)
        the_freqs <- temp$freqs
        the_ar_sdfe <- temp$sdfe
        ## find one or more peaks
         if(N_peaks == 1) freqs_of_biggest_peaks <- 2*pi*the_freqs[which.max(the_ar_sdfe)]
         else
         {
             ar_sdfe_extend <- c(the_ar_sdfe[2],the_ar_sdfe,the_ar_sdfe[length(the_ar_sdfe)-1])
             freq_extend <- c(NA,the_freqs,NA)
             N_ar_sdfee <- length(ar_sdfe_extend)
             temp_1 <- ar_sdfe_extend[-c(1,N_ar_sdfee)] > ar_sdfe_extend[-c(N_ar_sdfee-1,N_ar_sdfee)]
             temp_2 <- ar_sdfe_extend[-c(1,N_ar_sdfee)] > ar_sdfe_extend[-c(1,2)]
             peaks <- ar_sdfe_extend[-c(1,N_ar_sdfee)] > ar_sdfe_extend[-c(N_ar_sdfee-1,N_ar_sdfee)] & ar_sdfe_extend[-c(1,N_ar_sdfee)] > ar_sdfe_extend[-c(1,2)]
             freqs_of_biggest_peaks <- 2*pi*(the_freqs[peaks])[sort(rev(order(the_ar_sdfe[peaks]))[1:N_peaks])]
         }
        for(k in 1:N_peaks)
        {
            omega_j_minus_1 <- freqs_of_biggest_peaks[k]
            for(j in 1:J)
            {
                omega_j <- omega_j_minus_1 - g_prime(omega_j_minus_1)/g_prime_prime(omega_j_minus_1)
                omega_j_minus_1 <- omega_j
            }
            peaks_ar[k,n] <- omega_j/(2*pi*delta_t)
        }
    }
    return(peaks_ar)
}

### one frequency version (deprecated)
### 
### find_max_peak_ar <- function(tss,p=16,J=10,ar_method=ar_fbls,pad_factor=512,N_pad=NULL,center=TRUE,delta_t=1)
### {
###     g_prime <- function(omega) - sum(the_acvs*taus*sin(omega*taus))
###     g_prime_prime <- function(omega) - sum(the_acvs*taus^2*cos(omega*taus))
###     tss <- as_matrix(tss)
###     N_reps <- ncol(tss)
###     N      <- nrow(tss)
###     peaks_ar <- rep(NA,N_reps)
###     taus <- 0:p
###     for(k in 1:N_reps)
###     {
###         a_ts <- tss[,k]
###         the_ar_coeffs <- ar_method(a_ts,p,center=center)$coeffs
###         temp <- ar_coeffs_to_sdf(the_ar_coeffs,pad_factor=pad_factor,N_pad=N_pad)
###         omega_j_minus_1 <- 2*pi*temp$freqs[which.max(temp$sdfe)]
###         the_acvs <- acvs(c(-1,the_ar_coeffs),center=FALSE)$acvs
###         for(j in 1:J)
###         {
###             omega_j <- omega_j_minus_1 - g_prime(omega_j_minus_1)/g_prime_prime(omega_j_minus_1)
###             omega_j_minus_1 <- omega_j
###         }
###         peaks_ar[k] <- omega_j/(2*pi*delta_t)
###     }
###     return(peaks_ar=peaks_ar)
### }
