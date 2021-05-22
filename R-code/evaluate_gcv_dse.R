###

evaluate_gcv_dse <- function(the_dse,true_sdf,the_ms=some_ms_ff,the_measure=mse)
{
    N_dse <- length(the_dse)
    N_ms <- length(the_ms)
    the_dft_dse <- dft(the_dse)
    i_best_gridded_m <- which.min(sapply(the_ms,m_to_measure,the_dft_dse,true_sdf,the_measure))
    if(i_best_gridded_m-1 >= 1 && i_best_gridded_m+1 <= N_ms)
    {
        best_overall_m <- optimize(m_to_measure,the_ms[c(i_best_gridded_m-1,i_best_gridded_m+1)], the_dft_dse, true_sdf, the_measure)$min
        best_measure <- m_to_measure(best_overall_m,the_dft_dse,true_sdf,the_measure)
        i_gcv_best_gridded_m <- which.min(sapply(the_ms,function(an_m) gcv_dsp(an_m,the_dse)$gcv))
        if(i_gcv_best_gridded_m-1 >= 1 && i_gcv_best_gridded_m+1 <= N_ms)
        {
            gcv_best_overall_m <- optimize(function(an_m) gcv_dsp(an_m,the_dse)$gcv, the_ms[c(i_gcv_best_gridded_m-1,i_gcv_best_gridded_m+1)])$min
            gcv_measure <- m_to_measure(gcv_best_overall_m,the_dft_dse,true_sdf,the_measure)
            return(list(gcv_measure=gcv_measure,
                        best_measure=best_measure,
                        ratio=gcv_measure/best_measure,
                        gcv_m=gcv_best_overall_m,
                        best_m=best_overall_m))
        }
        else
            return(list(gcv_measure=NA,
                        best_measure=best_measure,
                        ratio=NA,
                        gcv_m=NA,
                        best_m=best_overall_m))
    }
    else 
        return(list(gcv_measure=NA,
                    best_measure=NA,
                    ratio=NA,
                    gcv_m=NA,
                    best_m=NA))
}


### function used in evaluation of GCV scheme ...

m_to_measure <- function(my_m, my_dft_dse, my_truth, my_measure)
{
    N_dse <- length(my_dft_dse)
    N_truth <- length(my_truth)
    temp <- dnorm(0:(N_dse/2-1),sd=1/(sqrt(2)*my_m))
    g_raw <- c(temp,0,rev(temp[-1]))
    dft_g_m <- dft(g_raw/sum(g_raw))
    return(my_measure(Re(inverse_dft(my_dft_dse*dft_g_m))[1:N_truth],my_truth))
}
