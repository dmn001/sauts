###

find_pgram_3_dB_width <- function(ts,f_peak,delta_t=1,center=TRUE,N_pad=4096)
{
    N <- length(ts)
    temp <- pgram(ts,pad=N_pad/N,delta_t=delta_t,center=center)
    freqs <- temp$freqs
    sdfe  <- temp$sdfe
    sdfe_peak <- pgram_single_freqs(f_peak,ts,delta_t=delta_t,center=center)$sdfe
    j_peak <- which.min(abs(f_peak-freqs))
    for(j in j_peak:1)
    {
        if(dB(sdfe_peak/sdfe[j]) > 3)
        {
            j_left <- j
            break
        }
    }
    f_left <- optimize(function(x) abs(dB(pgram_single_freqs(x,ts,delta_t=delta_t,center=center)$sdfe/sdfe_peak)-3),c(freqs[j_left],freqs[j_left+1]))$min
    for(j in j_peak:length(freqs))
    {
        if(dB(sdfe_peak/sdfe[j]) > 3)
        {
            j_right <- j
            break
        }
    }
    f_right <- optimize(function(x) abs(dB(pgram_single_freqs(x,ts,delta_t=delta_t,center=center)$sdfe/sdfe_peak)-3),c(freqs[j_right-1],freqs[j_right]))$min
    return(list(j_left=j_left,
                j_right=j_right,
                f_left=f_left,
                f_right=f_right,
                width_3dB=f_right-f_left,
                temp=temp))
}
