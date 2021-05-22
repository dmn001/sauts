### sin(t)/t taper (motivated by Equation (30) of Tsen81a,
### which references Geck78a);
### implementation is such that weight of zero would be assigned to values
### just before and just after actual time series;
### i.e., weights are all nonzero; this window maximimzes
### the area under the main lobe

sin_t_over_t_taper <- function(N)
{
    unnorm <- sapply(1:N,function(n) {temp <- -pi+n*2*pi/(N+1); if(isTRUE(all.equal(temp,0))) 1 else sin(temp)/temp})
    return(unnorm/sqrt(sum(unnorm^2)))
}
