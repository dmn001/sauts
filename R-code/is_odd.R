### is n an odd integer?

is_odd <- function(n) isTRUE(all.equal(n/2 - floor(n/2),0.5))
