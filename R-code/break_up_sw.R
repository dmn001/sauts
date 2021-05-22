### function for breaking up smoothing window with negative sidelobes
### into positive & negative components

break_up_sw <- function(sw,patch_with_this=1.0e-10)
{
    neg_guys <- which(sw$sw <= 0)
    pos_guys <- which(sw$sw > 0)
    neg_guys_rev <- which(rev(sw$sw) <= 0)
    pos_guys_rev <- which(rev(sw$sw) > 0)
    temp_pos <- c(neg_guys[which(diff(neg_guys)>1)]+1,pos_guys[which(diff(pos_guys)>1)])
    temp_neg <- c(neg_guys[which(diff(neg_guys)>1)],pos_guys[which(diff(pos_guys)>1)]+1)
    temp_pos_rev <- c(neg_guys_rev[which(diff(neg_guys_rev)>1)]+1,pos_guys_rev[which(diff(pos_guys_rev)>1)])
    temp_neg_rev <- c(neg_guys_rev[which(diff(neg_guys_rev)>1)],pos_guys_rev[which(diff(pos_guys_rev)>1)]+1)
    N_lw <- length(sw$sw)
    patch_these_pos <- unique(c(temp_pos,N_lw+1-temp_pos_rev))
    patch_these_neg <- unique(c(temp_neg,N_lw+1-temp_neg_rev))
    patch_me <- sw$sw
    patch_me[patch_these_pos] <-  patch_with_this
    patch_me[patch_these_neg] <- -patch_with_this
    return(list(sw_pos=list(freqs=sw$freqs[pos_guys],
                            sw=patch_me[pos_guys]),
                sw_neg=list(freqs=c(sw$freqs[neg_guys],0.5),
                            sw=c(patch_me[neg_guys],-patch_with_this))))
}
