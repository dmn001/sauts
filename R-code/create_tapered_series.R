### create a tapered time series
###
### following function takes a time series, centers it by subtracting
### its sample mean if center is TRUE, tapers the resulting series,
### recenters the tapered series is recenter is TRUE and normalizes
### the resulting series to have the same sample variance as the original
### time series if normalize is TRUE

create_tapered_series <- function(ts,the_taper,center=TRUE,recenter=FALSE,normalize=FALSE)
{
    ts_centered <- if(center) ts - mean(ts) else ts
    ts_tapered <- ts_centered*the_taper
    if(recenter) ts_tapered <- ts_tapered - mean(ts_tapered)
    if(normalize)
    {
        var_ts_centered <- sum(ts_tapered_centered^2)/length(ts)
        var_ts_tapered <- sum(ts_tapered^2)
        return(sqrt(var_ts_centered/var_ts_tapered)*ts_tapered)
    }
    else return(ts_tapered)
}
