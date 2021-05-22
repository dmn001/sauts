### compute SDF for AR process over one or more selected frequencies

ar_coeffs_to_sdf_single_freqs <- function(f,innov_var=0.002,coeffs=c(2.7607, -3.8106, 2.6535, -0.9238),delta_t=1)
{
    p <- length(coeffs)
    return(sapply(f,function(f) innov_var*delta_t/abs( 1- sum(coeffs*exp(complex(imag=-2*pi*f*delta_t*(1:p)))))^2))
}

### deprecated version that allows only one frequency
### 
### ar_coeffs_to_sdf_single_freq <- function(f,innov_var=0.002,coeffs=c(2.7607, -3.8106, 2.6535, -0.9238),delta_t=1)
### {
###     p <- length(coeffs)
###     innov_var*delta_t/abs( 1- sum(coeffs*exp(complex(imag=-2*pi*f*delta_t*(1:p)))))^2
### }



