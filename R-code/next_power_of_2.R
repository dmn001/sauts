### NOTE: WHEN TIME PERMITS, REDO THIS TO TAKE ADVANTAGE OF isTRUE(...,...))
###       SEE is.odd.R FOR AN EXAMPLE AND FOX & WEISBERG p. 103 DISCUSSION

next_power_of_2 <- function(n) 2^(ceiling(log2(n)))

### NOTE (15 Nov 2017): the following function (part of ifultools evidently)
### does much the same as next_power_of_2, but note the `x[x < 1]' business
###
### nextDyadic <- function(x)
### {
###   # returns the smallest integer n such that
###   # (i)  n is a power of two (e.g., 1, 2, 4, 8, 16, ...) and
###   # (ii) n is greater than or equal to x
###   x[x < 1] <- 1
###   as.integer(round(2^(ceiling(logb(x, base=2)))))
### }

