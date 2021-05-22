### functions for lag window-derived tapers

h_bartlett_lag_window <- function(t)
{
    abs_t <- abs(t)
    return(if(abs_t < 1) 1 - abs_t else 0)
}

### 

h_daniell_lag_window <- function(t,m=1)
{
    x <- pi*t/m
    return(if(t != 0) sin(x)/x else 1)
}

###

h_quadratic_lag_window <- function(t,m=1)
{
    x <- pi*t/m
    return(if(t != 0) 3*(sin(x)/x - cos(x))/x^2 else 1)
}

h_bartlett_priestley_lag_window <- h_quadratic_lag_window

###

h_parzen_lag_window <- function(t)
{
    abs_t <- abs(t)
    x <- abs_t
    return(if(abs_t <= 1/2)
               1 - 6*(abs_t^2 -abs_t^3)
           else
           {
               if(abs_t <= 1)
                   2*(1-abs_t)^3
               else
                   0
           })
}

### 

h_papoulis_lag_window <- function(t)
{
    abs_t <- abs(t)
    b <- pi*abs_t
    return(if(abs_t < 1) sin(b)/pi + (1-abs_t)*cos(b) else 0)
}

### 

h_gaussian_lag_window <- function(t,m=1)
{
    return(exp(-t^2/m^2))
}
