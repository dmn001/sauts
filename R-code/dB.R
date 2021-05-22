### dB and related functions

dB <- function(x) 10*log10(x)

careful_dB <- function(x, zero_mapping=-100.0)
  {
    temp <- function(y,zero_mapping) {if(y > 0) dB(y) else zero_mapping}
    if(length(x) <= 1) temp(x,zero_mapping)
    else
      sapply(x,temp,zero_mapping)
  }

convert_from_dB <- function(x) 10^(x/10)
