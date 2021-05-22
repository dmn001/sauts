###

mse     <- function(est,truth) sum((est-truth)^2)/length(est)
nmse    <- function(est,truth) sum(((est-truth)/truth)^2)/length(est)
msle    <- function(est,truth) sum((log(est)-log(truth))^2)/length(est)
kl      <- function(est,truth) sum((truth/est)-log(truth/est)-1)/length(est)
lk      <- function(est,truth) kl(truth,est)
kllk    <- function(est,truth) kl(est,truth) + kl(truth,est)
max_abs_dB <- function(est,truth) max(abs(dB(est/truth)))
