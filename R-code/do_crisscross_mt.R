### NOTE: multitapers are in the columns of the_tapers - each column should
###       have a sum of squares equal to unity

do_crisscross_mt <- function(the_tapers,delta_t=1,weights=NULL)
{
    K <- ncol(as.matrix(the_tapers))
    if(is.null(weights)) weights <- rep(1/K,K)
    edof <- 2/sum(weights^2)
    return(list(up = 10*log10(edof/qchisq(0.025,edof)),
                down = -10*log10(edof/qchisq(0.975,edof)),
                edof = edof,
                width = B_H_bar(the_tapers,delta_t,weights)))
}
