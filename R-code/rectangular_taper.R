### rectangular data taper (also known as the default data taper)

rectangular_taper <- function(N,parm=NULL) rep(1/sqrt(N),N)

default_taper <- rectangular_taper
