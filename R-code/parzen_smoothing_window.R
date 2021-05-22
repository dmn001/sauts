### NOTE: good for m even ...

parzen_smoothing_window <- function(f,m=2)
{
    if(f!=0)
    {
        x <- sin(pi*f)
        (4*(3-2*x^2)*(sin(m*pi*f/2))^4)/(m^3 * x^4)
    }
    else
    {
        3*m/4
    }
}


