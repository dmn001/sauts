### NOTE: WHEN TIME PERMITS, REDO THIS TO TAKE ADVANTAGE OF isTRUE(...,...))
###       SEE is.odd.R FOR AN EXAMPLE AND FOX & WEISBERG p. 103 DISCUSSION

is_power_of_2 <- function(x) if(!is.numeric(x) || (x <= 0) || (x < 1)) FALSE else {if(x==1) TRUE else is_power_of_2(x/2)}

