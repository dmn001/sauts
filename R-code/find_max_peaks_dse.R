###

find_max_peaks_dse <- function(tss,the_taper=default_taper(nrow(tss)),N_peaks=1,pad_factor=4,N_pad=NULL,center=TRUE,delta_t=1,J=10)
{
    g_prime <- function(omega) - sum(the_acvs*taus*sin(omega*taus))
    g_prime_prime <- function(omega) - sum(the_acvs*taus^2*cos(omega*taus))
    tss <- as.matrix(tss)
    N_time_series <- ncol(tss)
    N <- nrow(tss)
    peaks_dse <- matrix(nrow=N_peaks,ncol=N_time_series)
    taus <- 0:(N-1)
    for(n in 1:N_time_series)
    {
        a_ts <- tss[,n]
        the_acvs <- acvs(a_ts,center=center,taper=the_taper)$acvs
        temp <- direct_sdf_est(a_ts,the_taper=the_taper,center=center,pad=pad_factor,N_pad=N_pad)
        the_freqs <- temp$freqs
        the_dse <- temp$sdfe
        ## find one or more peaks
         if(N_peaks == 1) freqs_of_biggest_peaks <- 2*pi*the_freqs[which.max(the_dse)]
         else
         {
             dse_extend <- c(the_dse[2],the_dse,the_dse[length(the_dse)-1])
             freq_extend <- c(NA,the_freqs,NA)
             N_dsee <- length(dse_extend)
             temp_1 <- dse_extend[-c(1,N_dsee)] > dse_extend[-c(N_dsee-1,N_dsee)]
             temp_2 <- dse_extend[-c(1,N_dsee)] > dse_extend[-c(1,2)]
             peaks <- dse_extend[-c(1,N_dsee)] > dse_extend[-c(N_dsee-1,N_dsee)] & dse_extend[-c(1,N_dsee)] > dse_extend[-c(1,2)]
             freqs_of_biggest_peaks <- 2*pi*(the_freqs[peaks])[sort(rev(order(the_dse[peaks]))[1:N_peaks])]
         }
        for(k in 1:N_peaks)
        {
            omega_j_minus_1 <- freqs_of_biggest_peaks[k]
            for(j in 1:J)
            {
                omega_j <- omega_j_minus_1 - g_prime(omega_j_minus_1)/g_prime_prime(omega_j_minus_1)
                omega_j_minus_1 <- omega_j
            }
            peaks_dse[k,n] <- omega_j/(2*pi*delta_t)
        }
    }
    return(peaks_dse)
}

### one frequency version (deprecated)
### 
### find_max_peak_dse <- function(tss,J=10,the_taper=default_taper(nrow(tss)),center=TRUE,delta_t=1)
### {
###     g_prime <- function(omega) - sum(the_acvs*taus*sin(omega*taus))
###     g_prime_prime <- function(omega) - sum(the_acvs*taus^2*cos(omega*taus))
###     tss <- as.matrix(tss)
###     N_time_series <- ncol(tss)
###     N <- nrow(tss)
###     peaks_dse <- rep(NA,N_time_series)
###     taus <- 0:(N-1)
###     for(k in 1:N_time_series)
###     {
###         a_ts <- tss[,k]
###         the_acvs <- acvs(a_ts,center=center,taper=the_taper)$acvs
###         temp <- direct_sdf_est(a_ts,the_taper=the_taper,center=center,pad=4)
###         omega_j_minus_1 <- 2*pi*temp$freqs[which.max(temp$sdfe)]
###         for(j in 1:J)
###         {
###             omega_j <- omega_j_minus_1 - g_prime(omega_j_minus_1)/g_prime_prime(omega_j_minus_1)
###             omega_j_minus_1 <- omega_j
###         }
###         peaks_dse[k] <- omega_j/(2*pi*delta_t)
###     }
###     return(peaks_dse)
### }
