### R CODE FOR REPRODUCING CONTENT OF FIGURES AND TABLES IN CHAPTER 9 ...

ar2_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_1.txt")
ar4_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_1.txt")
ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_wave.txt")
ac_time_differences <- scan("http://faculty.washington.edu/dbp/sauts/Data/maser_deglitched.txt")
ac_fractional_frequencies <- diff(ac_time_differences)*100/6

### functions used to compute content of figures in Chapter 9 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_least_squares.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_reduced_likelihood.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/burg_algorithm.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H_bar.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/circular_shift.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/direct_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_lwe.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/gen_cis_for_ar_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_odd.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_windows.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ma_acvs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/mt_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/rectangular_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ts_to_lag_window_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/yule_walker_algorithm_given_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/yule_walker_algorithm_given_data.R")

### NOTE: to install the sapa library, uncomment the following three statements
###       and execute them:
###
### install.packages("devtools")
### devtools::install_github("wconstan/ifultools")
### devtools::install_github("wconstan/sapa")

library(sapa)

###

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figures 458, 459, 461 and 462 ###

fig_458 <- function(the_acvs,tag_1,tag_2=" ",coeffs=ar4_coeffs,innov_var=ar4_innov_var,main="Figure 458")
{
    ar_yw <- yule_walker_algorithm_given_acvs(the_acvs)
    est_sdf <- ar_coeffs_to_sdf(ar_yw$coeffs,ar_yw$innov_var,N_pad=1024)
    plot(est_sdf$freqs,dB(est_sdf$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="dB",
         typ="l",lwd=0.5, col="gray40",axes=FALSE,
         main=main)
    true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(true_sdf$freqs,dB(true_sdf$sdf),lwd=1.0)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,40,20),las=2)
    axis(2,at=seq(-60,40,10),label=FALSE,tcl=-0.25)
    text(0.25,-56.0,tag_1,pos=2)
    text(0.25,15.5,tag_2,pos=3)
    box(bty="l")
}

### Figure 458, top row

fig_458(acvs(ar4_1[1:16],center=FALSE)$acvs[1:5],expression(italic(N == 16)),expression(italic(p == 4)))
fig_458(acvs(ar4_1[1:64],center=FALSE)$acvs[1:5],expression(italic(N == 64)))

### Figure 458, bottom row

fig_458(acvs(ar4_1[1:256],center=FALSE)$acvs[1:5],expression(italic(N == 256)))
fig_458(acvs(ar4_1,center=FALSE)$acvs[1:5],expression(italic(N == 1024)))

### Figure 459, top row

fig_458(acvs(ar4_1[1:16],center=FALSE)$acvs[1:9],expression(italic(N == 16)),expression(italic(p == 8)),main="Figure 459")
fig_458(acvs(ar4_1[1:64],center=FALSE)$acvs[1:9],expression(italic(N == 64)),main="Figure 459")

### Figure 459, bottom row

fig_458(acvs(ar4_1[1:256],center=FALSE)$acvs[1:9],expression(italic(N == 256)),main="Figure 459")
fig_458(acvs(ar4_1,center=FALSE)$acvs[1:9],expression(italic(N == 1024)),main="Figure 459")

### Figure 461, top row

fig_458(acvs(ar4_1[1:16],taper=slepian_taper(16,1),center=FALSE)$acvs[1:5],expression(italic(N == 16)),expression(paste(italic(p == 4),", ",italic(NW == 1))),main="Figure 461")
fig_458(acvs(ar4_1[1:64],taper=slepian_taper(64,1),center=FALSE)$acvs[1:5],expression(italic(N == 64)),main="Figure 461")

### Figure 461, bottom row

fig_458(acvs(ar4_1[1:256],taper=slepian_taper(256,1),center=FALSE)$acvs[1:5],expression(italic(N == 256)),main="Figure 461")
fig_458(acvs(ar4_1,taper=slepian_taper(1024,1),center=FALSE)$acvs[1:5],expression(italic(N == 1024)),main="Figure 461")

### Figure 462, top row

fig_458(acvs(ar4_1[1:16],taper=slepian_taper(16,2),center=FALSE)$acvs[1:5],expression(italic(N == 16)),expression(paste(italic(p == 4),", ",italic(NW == 2))),main="Figure 462")
fig_458(acvs(ar4_1[1:64],taper=slepian_taper(64,2),center=FALSE)$acvs[1:5],expression(italic(N == 64)),main="Figure 462")

### Figure 462, bottom row

fig_458(acvs(ar4_1[1:256],taper=slepian_taper(256,2),center=FALSE)$acvs[1:5],expression(italic(N == 256)),main="Figure 462")
fig_458(acvs(ar4_1,taper=slepian_taper(1024,2),center=FALSE)$acvs[1:5],expression(italic(N == 1024)),main="Figure 462")

### Table 459 ###

ar4_1_coeffs_4 <- yule_walker_algorithm_given_acvs(acvs(ar4_1,center=FALSE)$acvs[1:5])$coeffs
ar4_1_coeffs_8 <- yule_walker_algorithm_given_acvs(acvs(ar4_1,center=FALSE)$acvs[1:9])$coeffs

### Table 459, top to bottom rows

round(c(ar4_coeffs,rep(0,4)),4)
round(c(ar4_1_coeffs_4,rep(0,4)),4)
round(ar4_1_coeffs_8,4)
round(c(abs(ar4_1_coeffs_4-ar4_coeffs),rep(0,4)),4)
round(abs(ar4_1_coeffs_8-c(ar4_coeffs,rep(0,4))),4)

### Figure 469 ###

fig_469 <- function(ts,tag_1,tag_2=" ",coeffs=ar4_coeffs,innov_var=ar4_innov_var,p=4,upper_y_lim=20)
{
    ar_burg <- burg_algorithm(ts,p,center=FALSE)
    est_sdf <- ar_coeffs_to_sdf(ar_burg$coeffs,ar_burg$innov_var,N_pad=1024)
    plot(est_sdf$freqs,dB(est_sdf$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,upper_y_lim),yaxs="i",ylab="dB",
         typ="l",lwd=0.5, col="gray40",axes=FALSE,
         main="Figure 469")
    true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(true_sdf$freqs,dB(true_sdf$sdf),lwd=1.0)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,40,20),las=2)
    axis(2,at=seq(-60,40,10),label=FALSE,tcl=-0.25)
    text(0.25,-48.0,tag_1,pos=1)
    text(0.25,15.5,tag_2,pos=3)
    box(bty="l")
}

### Figure 469, top row

fig_469(ar4_1[1:16],expression(italic(N == 16)),expression(italic(p == 4)),upper_y_lim=30)
fig_469(ar4_1[1:64],expression(italic(N == 64)),upper_y_lim=30)

### Figure 469, bottom row

fig_469(ar4_1[1:256],expression(italic(N == 256)))
fig_469(ar4_1,expression(italic(N == 1024)))

### Figures 475 and 476 ###

fig_475_acvs <- function(the_acvs,tag,vline=NULL,main="Figure 475")
{
    N <- length(the_acvs)
    taus <- 0:(N-1)
    plot(taus,the_acvs,
         xlim=c(-0.75,32.75),xaxs="i",xlab=expression(tau),
         ylim=c(-0.6,1.05),yaxs="i",ylab="ACVS",
         typ="p",axes=FALSE,
         main=main)
    abline(h=0,lty="dashed")
    if(!is.null(vline)) abline(v=vline,lty="dashed")
    axis(1,at=seq(0,32,16))
    axis(1,at=seq(0,32,4),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-1,1,1),las=2)
    axis(2,at=seq(-1,1,0.5),label=FALSE,tcl=-0.25)
    text(33,0.75,tag,pos=2,cex=1.2)
    box(bty="l")
}

fig_475_sdf <- function(the_sdf,upper_y_lim=18,main="Figure 475")
{
    plot(the_sdf$freqs,the_sdf$sdf,
         xlim=c(0,0.1),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,upper_y_lim),yaxs="i",ylab="SDF",
         typ="l",axes=FALSE,
         main=main)
    axis(1,at=seq(0,0.1,0.1))
    axis(1,at=seq(0,0.01,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,10),las=2)
    axis(2,at=seq(0,20,5),label=FALSE,tcl=-0.25)
    box(bty="l")
}

acvs_at_tau <- function(tau,alpha=0)
{
    abs_tau <- abs(tau)
    if(abs_tau <= 8) 1 - abs_tau/8
    else
    {
        if(abs_tau <= 24) alpha*(1 - abs(abs_tau - 16)/8)
        else
            0
    }
}


### Figure 475, left-hand column

fig_475_acvs(sapply(0:32,acvs_at_tau),expression(alpha == 0),8)
fig_475_acvs(sapply(0:32,acvs_at_tau,0.5),expression(alpha == 0.5),8)
fig_475_acvs(sapply(0:32,acvs_at_tau,-0.5),expression(alpha == -0.5),8)
fig_475_acvs(ar_coeffs_to_acvs(yule_walker_algorithm_given_acvs(sapply(0:8,acvs_at_tau))$coeffs, max_lag=32),"maximum entropy",8)

### Figure 475, right-hand column

fig_475_sdf(ma_acvs_to_sdf(c(sapply(0:7,acvs_at_tau),rep(0,1024-8))))
fig_475_sdf(ma_acvs_to_sdf(c(sapply(0:23,acvs_at_tau,0.5),rep(0,1024-24))))
fig_475_sdf(ma_acvs_to_sdf(c(sapply(0:23,acvs_at_tau,-0.5),rep(0,1024-24))))
ar_8 <- yule_walker_algorithm_given_acvs(sapply(0:8,acvs_at_tau))
fig_475_sdf(ar_coeffs_to_sdf(ar_8$coeffs,ar_8$innov_var,N_pad=1024))

### Figure 476, left-hand column

fig_475_acvs(sapply(0:32,acvs_at_tau),"MA(7)",main="Figure 476")
fig_475_acvs(ar_coeffs_to_acvs(yule_walker_algorithm_given_acvs(sapply(0:8,acvs_at_tau))$coeffs, max_lag=32),"AR(8)",8,main="Figure 476")
fig_475_acvs(ar_coeffs_to_acvs(yule_walker_algorithm_given_acvs(sapply(0:9,acvs_at_tau))$coeffs, max_lag=32),"AR(9)",9,main="Figure 476")
fig_475_acvs(ar_coeffs_to_acvs(yule_walker_algorithm_given_acvs(sapply(0:16,acvs_at_tau))$coeffs, max_lag=32),"AR(16)",16,main="Figure 476")

### Figure 476, right-hand column

fig_475_sdf(ma_acvs_to_sdf(c(sapply(0:7,acvs_at_tau),rep(0,1024-8))),upper_y_lim=12,main="Figure 476")
fig_475_sdf(ar_coeffs_to_sdf(ar_8$coeffs,ar_8$innov_var,N_pad=1024),upper_y_lim=12,main="Figure 476")
ar_9 <- yule_walker_algorithm_given_acvs(sapply(0:9,acvs_at_tau))
fig_475_sdf(ar_coeffs_to_sdf(ar_9$coeffs,ar_9$innov_var,N_pad=1024),upper_y_lim=12,main="Figure 476")
ar_16 <- yule_walker_algorithm_given_acvs(sapply(0:16,acvs_at_tau))
fig_475_sdf(ar_coeffs_to_sdf(ar_16$coeffs,ar_16$innov_var,N_pad=1024),upper_y_lim=12,main="Figure 476")

### Figure 478 ###

fig_478 <- function(ts,tag_1,tag_2=" ",coeffs=ar4_coeffs,innov_var=ar4_innov_var,p=4,upper_y_lim=20)
{
    est_parms <- ar_fbls(ts,p,center=FALSE)
    est_sdf <- ar_coeffs_to_sdf(est_parms$coeffs,est_parms$innov_var,N_pad=1024)
    plot(est_sdf$freqs,dB(est_sdf$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,upper_y_lim),yaxs="i",ylab="dB",
         typ="l",lwd=0.5, col="gray40",axes=FALSE,
         main="Figure 478")
    true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(true_sdf$freqs,dB(true_sdf$sdf),lwd=1.0)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,40,20),las=2)
    axis(2,at=seq(-60,40,10),label=FALSE,tcl=-0.25)
    text(0.25,-48.0,tag_1,pos=1)
    text(0.25,15.5,tag_2,pos=3)
    box(bty="l")
}

### Figure 478, top row

fig_478(ar4_1[1:16],expression(italic(N == 16)),expression(italic(p == 4)),upper_y_lim=30)
fig_478(ar4_1[1:64],expression(italic(N == 64)),upper_y_lim=30)

### Figure 478, bottom row

fig_478(ar4_1[1:256],expression(italic(N == 256)))
fig_478(ar4_1,expression(italic(N == 1024)))

### Figure 483 ###

fig_483 <- function(ts,tag_1,tag_2=" ",coeffs=ar4_coeffs,innov_var=ar4_innov_var,p=4,upper_y_lim=20)
{
    est_parms <- ar.mle(ts,aic=FALSE,order=p,demean=FALSE)
    est_sdf <- ar_coeffs_to_sdf(est_parms$ar,est_parms$var.pred,N_pad=1024)
    plot(est_sdf$freqs,dB(est_sdf$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,upper_y_lim),yaxs="i",ylab="dB",
         typ="l",lwd=0.5, col="gray40",axes=FALSE,
         main="Figure 483")
    true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(true_sdf$freqs,dB(true_sdf$sdf),lwd=1.0)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,40,20),las=2)
    axis(2,at=seq(-60,40,10),label=FALSE,tcl=-0.25)
    text(0.25,-48.0,tag_1,pos=1)
    text(0.25,15.5,tag_2,pos=3)
    box(bty="l")
}

### Figure 483, top row

fig_483(ar4_1[1:16],expression(italic(N == 16)),expression(italic(p == 4)),upper_y_lim=30)
fig_483(ar4_1[1:64],expression(italic(N == 64)),upper_y_lim=30)

### Figure 483, bottom row

fig_483(ar4_1[1:256],expression(italic(N == 256)))
fig_483(ar4_1,expression(italic(N == 1024)))

### Table 483, Yule-Walker column

ar_reduced_likelihood(ar4_1[1:16],yule_walker_algorithm_given_data(ar4_1[1:16],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1[1:64],yule_walker_algorithm_given_data(ar4_1[1:64],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1[1:256],yule_walker_algorithm_given_data(ar4_1[1:256],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1,yule_walker_algorithm_given_data(ar4_1,4,center=FALSE)$coeffs,center=FALSE)

### Table 483, Burg column

ar_reduced_likelihood(ar4_1[1:16],burg_algorithm(ar4_1[1:16],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1[1:64],burg_algorithm(ar4_1[1:64],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1[1:256],burg_algorithm(ar4_1[1:256],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1,burg_algorithm(ar4_1,4,center=FALSE)$coeffs,center=FALSE)

### Table 483, FBLS column

ar_reduced_likelihood(ar4_1[1:16],ar_fbls(ar4_1[1:16],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1[1:64],ar_fbls(ar4_1[1:64],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1[1:256],ar_fbls(ar4_1[1:256],4,center=FALSE)$coeffs,center=FALSE)
ar_reduced_likelihood(ar4_1,ar_fbls(ar4_1,4,center=FALSE)$coeffs,center=FALSE)

### Table 483, ML column

ar_reduced_likelihood(ar4_1[1:16],ar.mle(ar4_1[1:16],aic=FALSE,order=4,demean=FALSE)$ar,center=FALSE)
ar_reduced_likelihood(ar4_1[1:64],ar.mle(ar4_1[1:64],aic=FALSE,order=4,demean=FALSE)$ar,center=FALSE)
ar_reduced_likelihood(ar4_1[1:256],ar.mle(ar4_1[1:256],aic=FALSE,order=4,demean=FALSE)$ar,center=FALSE)
ar_reduced_likelihood(ar4_1,ar.mle(ar4_1,aic=FALSE,order=4,demean=FALSE)$ar,center=FALSE)

### Figures 489 and 490 ###

fig_489 <- function(ts,coeffs,innov_var,tag_1,tag_2=" ",main="Figure 489")
{
    p <- length(coeffs)
    est_parms <- ar.mle(ts,aic=FALSE,order=p,demean=FALSE)
    est_sdf <- ar_coeffs_to_sdf(est_parms$ar,est_parms$var.pred,N_pad=1024)
    plot(est_sdf$freqs,dB(est_sdf$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(if(p==2) -20 else -60,40),yaxs="i",ylab="dB",
         typ="l",lwd=0.5, col="gray40",axes=FALSE,
         main=main)
    true_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(true_sdf$freqs,dB(true_sdf$sdf),lwd=1.0)
    cis <- gen_cis_for_ar_sdf(length(ts),est_parms$ar,est_parms$var.pred,N_pad=2048)
    lines(cis$freqs,dB(cis$lower),lty="dotted",lwd=1.0)
    lines(cis$freqs,dB(cis$upper),lty="dotted",lwd=1.0)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,40,20),las=2)
    axis(2,at=seq(-60,40,10),label=FALSE,tcl=-0.25)
    text(0.075,if(p==2) -17.0 else -55.0,tag_1,pos=4)
    text(0.25,if(p==2) 36.25 else 33.75,tag_2,pos=3)
    box(bty="l")
}

### Figure 489, top row

fig_489(ar2_1[1:16],ar2_coeffs,ar2_innov_var,expression(italic(N == 16)),expression(italic(p == 2)))
fig_489(ar2_1[1:64],ar2_coeffs,ar2_innov_var,expression(italic(N == 64)))

### Figure 489, bottom row

fig_489(ar2_1[1:256],ar2_coeffs,ar2_innov_var,expression(italic(N == 256)))
fig_489(ar2_1,ar2_coeffs,ar2_innov_var,expression(italic(N == 1024)))

### Figure 490, top row

fig_489(ar4_1[1:16],ar4_coeffs,ar4_innov_var,expression(italic(N == 16)),expression(italic(p == 4)),main="Figure 490")
fig_489(ar4_1[1:64],ar4_coeffs,ar4_innov_var,expression(italic(N == 64)),main="Figure 490")

### Figure 490, bottom row

fig_489(ar4_1[1:256],ar4_coeffs,ar4_innov_var,expression(italic(N == 256)),main="Figure 490")
fig_489(ar4_1,ar4_coeffs,ar4_innov_var,expression(italic(N == 1024)),main="Figure 490")

### Figure 496 ###

fig_496 <- function(ts,p=5,taper=slepian_taper(length(ts),2),delta_t=1/4)
{
    yw_parms <- yule_walker_algorithm_given_data(ts,p,center=TRUE)
    yw_sdf <- ar_coeffs_to_sdf(yw_parms$coeffs,yw_parms$innov_var,N_pad=1024,delta_t=delta_t)
    plot(yw_sdf$freqs,dB(yw_sdf$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 496")
    burg_parms <- burg_algorithm(ts,p,center=TRUE)
    burg_sdf <- ar_coeffs_to_sdf(burg_parms$coeffs,burg_parms$innov_var,N_pad=1024,delta_t=delta_t)
    lines(burg_sdf$freqs,dB(burg_sdf$sdfe),lty="dashed")
    dse <- direct_sdf_est(ts,taper,center=TRUE,delta_t=delta_t)
    lines(dse$freqs,dB(dse$sdfe),lwd=0.25,col="gray40")
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 496

fig_496(ocean_wave)

### Figure 497 ###

fig_497_a <- function(ts,p=5,tag="(a)",delta_t=1/4)
{
    pred_errors <- burg_algorithm(ts,p,center=TRUE)$forward[[p]]
    times <- (p:(length(ts)-1))*delta_t
    plot(times,pred_errors,
         xlim=c(-1,257),xaxs="i",xlab="time  (sec)",
         ylim=c(-80,80),ylab="prediction errors",
        ,type="l",lwd=0.25,col="gray40",axes=FALSE,
         main="Figure 497(a)")
    axis(1,at=seq(0,256,64))
    axis(2,at=seq(-80,80,40),las=2)
    axis(2,at=seq(-80,80,10),label=FALSE,tcl=-0.25)
    text(235,70,tag,pos=4)
    box(bty="l")
}

draw_cc_for_fig_497 <- function(cc,y_cc)
{
    x_cc <- 0.16
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
}

fig_497_b <- function(ts,p=5,m=55,tag="(b)",delta_t=1/4)
{
    pred_errors <- burg_algorithm(ts,p,center=TRUE)$forward[[p]]
    sdfe_raw <- pgram(pred_errors,delta_t=delta_t)
    plot(sdfe_raw$freqs,dB(sdfe_raw$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main="Figure 497(b)")
    sdfe_sm <- ts_to_lag_window_sdf_est(pred_errors,m=m,lag_window=parzen_lag_window,delta_t=delta_t)
    lines(sdfe_sm$freqs,dB(sdfe_sm$sdfe))
    draw_cc_for_fig_497(sdfe_sm$cc,-20)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.835,68.0,tag,pos=4)
    box(bty="l")
}

fig_497_c <- function(ts,taper=slepian_taper(length(ts),2),p=5,m=55,tag="(c)",delta_t=1/4)
{
    sdfe_raw <- direct_sdf_est(ts,taper,center=TRUE,delta_t=delta_t)
    plot(sdfe_raw$freqs,dB(sdfe_raw$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main="Figure 497(c)")
    temp <- burg_algorithm(ts,p,center=TRUE)
    burg_coeffs <- temp$coeffs
    pred_errors <- temp$forward[[p]]
    postcoloring <- ar_coeffs_to_sdf(burg_coeffs,N_pad=1019)
    sdfe_sm <- ts_to_lag_window_sdf_est(pred_errors,m=m,lag_window=parzen_lag_window,delta_t=delta_t)
    sdfe_postcolored <- sdfe_sm$sdfe[seq(1,1019,2)]*postcoloring$sdf
    lines(postcoloring$freqs/delta_t,dB(sdfe_postcolored))
    draw_cc_for_fig_497(sdfe_sm$cc,20)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.835,68.0,tag,pos=4)
    box(bty="l")
}

fig_497_d <- function(ts,taper=slepian_taper(length(ts),2),p=5,m_parzen=55,m_gaussian=23.666,tag="(d)",delta_t=1/4)
{
    sdfe_gaussian <- ts_to_lag_window_sdf_est(ts,m=m_gaussian,lag_window=gaussian_lag_window,taper=taper,delta_t=delta_t)
    plot(sdfe_gaussian$freqs,dB(sdfe_gaussian$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main="Figure 497(d)")
    draw_cc_for_fig_497(sdfe_gaussian$cc,0)
    temp <- burg_algorithm(ts,p,center=TRUE)
    burg_coeffs <- temp$coeffs
    pred_errors <- temp$forward[[p]]
    postcoloring <- ar_coeffs_to_sdf(burg_coeffs,N_pad=1019)
    sdfe_parzen <- ts_to_lag_window_sdf_est(pred_errors,m=m_parzen,lag_window=parzen_lag_window,delta_t=delta_t)
    sdfe_postcolored <- sdfe_parzen$sdfe[seq(1,1019,2)]*postcoloring$sdf
    lines(postcoloring$freqs/delta_t,dB(sdfe_postcolored))
    draw_cc_for_fig_497(sdfe_parzen$cc,20)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.835,68.0,tag,pos=4)
    box(bty="l")
}

### Figure 497, top to bottom plots

fig_497_a(ocean_wave)
fig_497_b(ocean_wave)
fig_497_c(ocean_wave)
fig_497_d(ocean_wave)

### Figure 499 ###

fig_499 <- function(ts,tag,K=40)
{
    N <- length(ts)
    which_ar <- if(tag=="Burg") burg else yule_walker
    ys <- which_ar(ts,K,center=TRUE)$pacs
    plot(1:K,ys,
         xlim=c(0,K+0.5),xaxs="i",xlab=expression(italic(k)),
         ylim=c(-1.05,1.05),yaxs="i",ylab="PACS",
         typ="h",axes=FALSE,
         main="Figure 499")
    abline(h=0,lwd=0.5)
    abline(h=c(-2,2)/sqrt(N),lwd=0.5,lty="dashed")
    axis(1,at=seq(0,K,5))
    axis(1,at=1:K,label=FALSE,tcl=-0.25)
    axis(2,at=seq(-1,1,1),las=2)
    axis(2,at=seq(-1,1,0.2),label=FALSE,tcl=-0.25)
    text(K/2+0.25,0.925,tag,pos=1)
    box(bty="l")
}

### Figure 499, top plot

fig_499(ocean_wave,"Yule-Walker")

### Figure 499, bottom plot

fig_499(ocean_wave,"Burg")

### Figure 500 ###

fig_500_a <- function(ts,p=25,delta_t=1/4)
{
    burg_est <- burg(ts,p,center=TRUE)
    burg_sdf <- ar_coeffs_to_sdf(burg_est$coeffs,burg_est$innov_var,N_pad=1024,delta_t=delta_t)
    plot(burg_sdf$freqs,dB(burg_sdf$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 500(a)")
    cis <- gen_cis_for_ar_sdf(length(ts),burg_est$coeffs,burg_est$innov_var,delta_t=delta_t)
    lines(cis$freqs,dB(cis$lower),lwd=0.25,col="gray40")
    lines(cis$freqs,dB(cis$upper),lwd=0.25,col="gray40")
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.97,65.6,"(a)",pos=2,cex=1.2)
    box(bty="l")
}

fig_500_b <- function(ts,p=25,p_prewhiten=5,m=55,delta_t=1/4)
{
    burg_est <- burg(ts,p,center=TRUE)
    burg_sdf <- ar_coeffs_to_sdf(burg_est$coeffs,burg_est$innov_var,N_pad=1024,delta_t=delta_t)
    plot(burg_sdf$freqs,dB(burg_sdf$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 500(b)")
    temp <- burg_algorithm(ts,p_prewhiten,center=TRUE)
    burg_coeffs <- temp$coeffs
    pred_errors <- temp$forward[[p_prewhiten]]
    postcoloring <- ar_coeffs_to_sdf(burg_coeffs,N_pad=1019)
    sdfe_parzen <- ts_to_lag_window_sdf_est(pred_errors,m=m,lag_window=parzen_lag_window,delta_t=delta_t)
    sdfe_postcolored <- sdfe_parzen$sdfe[seq(1,1019,2)]*postcoloring$sdf
    lines(postcoloring$freqs/delta_t,dB(sdfe_postcolored),lwd=0.25,col="gray40")
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.97,65.6,"(b)",pos=2,cex=1.2)
    box(bty="l")
}

### Figure 500, top plot

fig_500_a(ocean_wave)

### Figure 500, bottom plot

fig_500_b(ocean_wave)

### Table 500 ###

fd_N <- length(ac_fractional_frequencies) - 1
fd_p <- 41
fd_sig_sq_hat_yw <- yule_walker(diff(ac_fractional_frequencies),fd_p,center=TRUE)$innov_var
fd_sig_sq_hat_burg <- burg(diff(ac_fractional_frequencies),fd_p,center=TRUE)$innov_var

N <- length(ac_fractional_frequencies)
p <- 42
sig_sq_hat_yw <- yule_walker(ac_fractional_frequencies,p,center=TRUE)$innov_var
sig_sq_hat_burg <- burg(ac_fractional_frequencies,p,center=TRUE)$innov_var

### Table 500, top row, numbers from left to right

fd_p                                                     # 41
round(fd_sig_sq_hat_yw,5)                                # 0.02259
round(sqrt(2*fd_sig_sq_hat_yw^2/fd_N),5)                 # 0.00051
round(fd_sig_sq_hat_yw*(fd_N+fd_p+1)/(fd_N-fd_p-1),5)    # 0.02307

fd_p                                                     # 41
round(fd_sig_sq_hat_burg,5)                              # 0.0223
round(sqrt(2*fd_sig_sq_hat_burg^2/fd_N),5)               # 5e-04
round(fd_sig_sq_hat_burg*(fd_N+fd_p+1)/(fd_N-fd_p-1),5)  # 0.02277

### Table 500, bottom row, numbers from left to right

p                                                        # 42
round(sig_sq_hat_yw,5)                                   # 0.02222
round(sqrt(2*sig_sq_hat_yw^2/N),5)                       # 5e-04
round(sig_sq_hat_yw*(N+p+1)/(N-p-1),5)                   # 0.0227

p                                                        # 42
round(sig_sq_hat_burg,5)                                 # 0.02216
round(sqrt(2*sig_sq_hat_burg^2/N),5)                     # 5e-04
round(sig_sq_hat_burg*(N+p+1)/(N-p-1),5)                 # 0.02264

### Figure 501 ###

fig_501_a <- function(ff,p=42,m=80,tag="(a)",delta_t=1)
{
    burg_est <- burg(ff,p,center=TRUE)
    burg_sdf <- ar_coeffs_to_sdf(burg_est$coeffs,burg_est$innov_var,N_pad=1024,delta_t=delta_t)
    plot(burg_sdf$freqs,dB(burg_sdf$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/minute)")),
         ylim=c(-30,-5),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 501(a)")
    lwe <- ts_to_lag_window_sdf_est(ff,m,lag_window=gaussian_lag_window,center=TRUE)
    lines(lwe$freqs,dB(lwe$sdfe),lwd=0.5,col="gray")
    x_cc <- 0.3
    y_cc <- -25.0
    lines(c(x_cc,x_cc),y_cc+c(lwe$cc$up,-lwe$cc$down),lwd=0.5)
    lines(x_cc+c(-lwe$cc$width/2,lwe$cc$width/2),c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-30,0,10),las=2)
    axis(2,at=seq(-30,0,5),label=FALSE,tcl=-0.25)
    text(0.01,-6.5,tag,pos=4)
    box(bty="l")
}


fig_501_b <- function(ff,p=42,K=7,M=140,tag="(b)",x_upper_lim_log=log10(0.06))
{
    burg_est <- burg(ff,p,center=TRUE)
    burg_sdf <- ar_coeffs_to_sdf(burg_est$coeffs,burg_est$innov_var,N_pad=2*length(ff))
    plot(log10(burg_sdf$freqs[-1]),dB(burg_sdf$sdfe[-1]),
         xlim=c(log10(burg_sdf$freqs[2]),x_upper_lim_log),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/minute)")),
         ylim=c(-30,-5),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 501(b)")
    mtse <- mt_sdf_est(ff,t(as.matrix(taper("sine",length(ff),n.taper=K))),center=TRUE)
    log10_freqs <- log10(mtse$freqs[-1])         
    dB_mtse <- dB(mtse$sdfe[-1])        
    M_0 <- ceiling(K/2)
    reg_coeffs <- coef(lm(dB_mtse[M_0:M] ~ log10_freqs[M_0:M]))
    reg_xs <- log10(c(M_0,M)/length(ff))
    reg_ys <- reg_coeffs[1] + reg_coeffs[2]*reg_xs
    lines(reg_xs,reg_ys,lwd=1.5)
    ##
    dB_burg <- dB(burg_sdf$sdfe[seq(3,length(burg_sdf$freqs),2)])        
    reg_coeffs <- coef(lm(dB_burg[M_0:M] ~ log10_freqs[M_0:M]))
    reg_xs <- log10(c(M_0,M)/length(ff))
    reg_ys <- reg_coeffs[1] + reg_coeffs[2]*reg_xs
    lines(reg_xs,reg_ys,lwd=1.25,lty="dashed")
    ##
    axis(1, at=log10(c(1/100000,1/10000,1/1000,1/100,1/10,1)), labels=expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
    ttn <- 2:9
    axis(1, at=log10(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)), labels=FALSE, tcl=-0.25)
    axis(2,at=seq(-30,0,10),las=2)
    axis(2,at=seq(-30,0,5),label=FALSE,tcl=-0.25)
    text(-1.65,-6.5,tag,pos=4)
    box(bty="l")
}

### Figure 501, left-hand plot

fig_501_a(ac_fractional_frequencies)

### Figure 501, right-hand plot

fig_501_b(ac_fractional_frequencies)

### Figure 508 ###

fig_508 <- function(ts,tag)
{
    sdfe <- if(tag=="(a)") pgram(ts,center=FALSE)
            else
            {
                yw_est <- yule_walker(ts,length(ts)-1,center=FALSE)
                ar_coeffs_to_sdf(yw_est$coeffs,yw_est$innov_var,N_pad=1024)
            }
    plot(sdfe$freqs,dB(sdfe$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/minute)")),
         ylim=c(-40,20),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main=paste("Figure 508",tag,sep=""))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(0.45,15,tag,pos=1)
    box(bty="l")
}

### Figure 508, left-hand plot

fig_508(ar4_1,"(a)")

### Figure 508, right-hand plot

fig_508(ar4_1,"(b)")
