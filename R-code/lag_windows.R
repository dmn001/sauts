### lag windows

default_lag_window <- function(tau,m=NULL) 1

###

bartlett_lag_window <- function(tau,m=1)
{
    a_tau <- abs(tau)
    if(a_tau < m) 1 - a_tau/m else 0
}

###

daniell_lag_window <- function(tau,m=1)
{
    x <- pi*tau/m
    if(tau != 0) sin(x)/x else 1
}

###

parzen_lag_window <- function(tau,m=1)
{
    a_tau <- abs(tau)
    x <- a_tau/m
    if(a_tau <= m/2)
        1 - 6*(x^2 -x^3)
    else
    {
        if(a_tau <= m)
            2*(1-x)^3
        else
            0
    }
}

###

papoulis_lag_window <- function(tau,m=1)
{
    a_tau <- abs(tau)
    a_tau_o_m <- a_tau/m
    b <- pi*a_tau_o_m
    if(a_tau < m) sin(b)/pi + (1-a_tau_o_m)*cos(b) else 0
}

###

gaussian_lag_window <- function(tau,m=1)
{
    exp(-tau^2/m^2)
}

###

quadratic_lag_window <- function(tau,m=1)
{
    x <- pi*tau/m
    if(tau != 0) 3*(sin(x)/x - cos(x))/x^2 else 1
}

bartlett_priestley_lag_window <- quadratic_lag_window
epanechnikov_lag_window       <- quadratic_lag_window

### 

split_cosine_lag_window <- function(tau,parms=c(1,1))
{
    m <- parms[1]
    p <- if(length(parms) == 1) 1 else parms[2]
    a_tau <- abs(tau)
    M <- floor(p*(2*m-1)/2) + 1
    if(a_tau <= m-M)
        1
    else
    {
        if(a_tau <= m-1)
            0.5*(1-cos(pi*(m-a_tau)/M))
        else
            0
    }
}

### 

modified_daniell_lag_window <- function(tau,parms=c(3,64))
{
    M <- parms[1]
    N <- parms[2]
    g_j <- c(rep(1/(2*M),M),1/(4*M))
    for_dft <- rep(0,2*N)
    for_dft[1:(M+1)] <- g_j
    for_dft[((2*N)-M+1):(2*N)] <- rev(g_j[-1])
    return(Re(dft(for_dft))[tau+1])
}

modified_daniell_lag_window_one_sided <- function(M,N)
{
    g_j <- c(rep(1/(2*M),M),1/(4*M))
    for_dft <- rep(0,2*N)
    for_dft[1:(M+1)] <- g_j
    for_dft[((2*N)-M+1):(2*N)] <- rev(g_j[-1])
    return(Re(dft(for_dft))[1:N])
}

### 

reshaped_lag_window <- function(tau,parms=list(37,64,parzen_lag_window,slepian_taper))
{
    a_tau <- abs(tau)
    m <- parms[[1]]
    N <- parms[[2]]
    lag_window_func <- parms[[3]]
    taper_func <- parms[[4]]
    the_taper <- taper_func(N)
    return(lag_window_func(tau,m)*sum(the_taper[1:(N-a_tau)]*the_taper[(1+a_tau):N])/(1-a_tau/N))
}
