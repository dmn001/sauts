###

find_ar_3_dB_width <- function(f_peak,innov_var,coeffs,delta_t=1,N_pad=65536)
{
    the_sdfe <- ar_coeffs_to_sdf(coeffs,innov_var,delta_t=delta_t,N_pad=N_pad)
    sdfe_peak <- ar_coeffs_to_sdf_single_freqs(f_peak,innov_var=innov_var,coeffs=coeffs,delta_t=delta_t)
    freqs <- the_sdfe$freqs
    sdfe  <- the_sdfe$sdfe
    j_peak <- which.min(abs(f_peak-freqs))
    for(j in j_peak:1)
    {
        if(dB(sdfe_peak/sdfe[j]) > 3)
        {
            j_left <- j
            break
        }
    }
    f_left <- optimize(function(x) abs(dB(ar_coeffs_to_sdf_single_freqs(x,innov_var=innov_var,coeffs=coeffs,delta_t=delta_t)/sdfe_peak)-3),c(freqs[j_left],freqs[j_left+1]))$min
    for(j in j_peak:length(freqs))
    {
        if(dB(sdfe_peak/sdfe[j]) > 3)
        {
            j_right <- j
            break
        }
    }
    f_right <- optimize(function(x) abs(dB(ar_coeffs_to_sdf_single_freqs(x,innov_var=innov_var,coeffs=coeffs,delta_t=delta_t)/sdfe_peak)-3),c(freqs[j_right-1],freqs[j_right]))$min
    return(list(j_left=j_left,
                j_right=j_right,
                f_left=f_left,
                f_right=f_right,
                width_3dB=f_right-f_left))
}
