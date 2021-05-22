### implementation of DFT and inverse DFT as defined in SAUTS
### using R's fft function (note need for renormalization in
### inverse_dft by division by length(x) in addition to setting
### inverse=TRUE)

dft <- function(x) fft(x)

inverse_dft <- function(x) fft(x,inverse=TRUE)/length(x)
