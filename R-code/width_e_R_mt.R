### multitaper width measure of Equation (355a)

width_e_R_mt <- function(the_tapers,delta_t=1)
{
    the_sq_tapers <- as.matrix(the_tapers)^2
    K <- ncol(the_sq_tapers)
    return(sum(rowSums(the_sq_tapers)^2)/(K*delta_t))
}
