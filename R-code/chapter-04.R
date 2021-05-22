### R CODE FOR REPRODUCING CONTENT OF FIGURES IN CHAPTER 4 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 121 ###

fig_121_IS <- function(int_spec)
{
    plot(int_spec[[1]]$freqs,int_spec[[1]]$spec,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-0.02,1.005),yaxs="i",ylab="integrated spectrum",
         typ="l",axes=FALSE,
         main="Figure 121")
    if(length(int_spec) > 1) for(k in 2:length(int_spec)) lines(int_spec[[k]]$freqs,int_spec[[k]]$spec,lty=c("dotted","solid","dotted","solid")[k-1])
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(2,at=seq(0,1,1),las=2)
    box(bty="l")
}

phi <- 0.49
N_pad <- 2048
freqs_sdf <- (0:(N_pad/2))/N_pad
ar_sdf <- (1-phi^2)/abs(fft(c(1,-phi,rep(0,N_pad-2)))[1:((N_pad/2)+1)])^2
ar_sdf_two_sided <- c(rev(ar_sdf[-1]), ar_sdf)
freqs_sdf_two_sided <- c(-rev(freqs_sdf[-1]),freqs_sdf)
ar_ispec <- freqs_sdf[2]*cumsum(ar_sdf_two_sided)[-length(ar_sdf_two_sided)]
freqs_ispec_two_sided <- freqs_sdf_two_sided[-1]
int_spec_pure_cont <- list(list(freqs=freqs_ispec_two_sided,spec=ar_ispec))

freq <- 1/16
int_spec_pure_disc <- list(list(freqs=c(-1/2,-freq),
                                spec=rep(0,2)),
                           list(freqs=rep(-freq,2),
                                spec=c(0,1/2)),
                           list(freqs=c(-freq,freq),
                                spec=rep(1/2,2)),
                           list(freqs=rep(freq,2),
                                spec=c(1/2,1)),
                           list(freqs=c(freq,1/2),
                                spec=rep(1,2)))

int_spec_disc <- list(list(freqs=c(-1/2,-freq),
                           spec=c(0,7/32)),
                      list(freqs=rep(-freq,2),
                           spec=c(7/32,15/32)),
                      list(freqs=c(-freq,freq),
                           spec=c(15/32,17/32)),
                      list(freqs=rep(freq,2),
                           spec=c(17/32,25/32)),
                      list(freqs=c(freq,1/2),
                           spec=c(25/32,1)))

j_1 <- which(freqs_ispec_two_sided==-1/16)
j_2 <- which(freqs_ispec_two_sided==1/16)
int_spec_mixed <- list(list(freqs=freqs_ispec_two_sided[1:j_1],
                           spec=0.5*ar_ispec[1:j_1]),
                      list(freqs=rep(-freq,2),
                           spec=c(0.5*ar_ispec[j_1],0.25+0.5*ar_ispec[j_1])),
                      list(freqs=freqs_ispec_two_sided[j_1:j_2],
                           spec=0.25+0.5*ar_ispec[j_1:j_2]),
                      list(freqs=rep(freq,2),
                           spec=c(0.25+0.5*ar_ispec[j_2],0.5+0.5*ar_ispec[j_2])),
                      list(freqs=freqs_ispec_two_sided[j_2:length(ar_ispec)],
                           spec=0.5+0.5*ar_ispec[j_2:length(ar_ispec)]))

### Figure 121, top row

fig_121_IS(int_spec_pure_cont)

fig_121_IS(int_spec_pure_disc)

fig_121_IS(int_spec_disc)

fig_121_IS(int_spec_mixed)

###

fig_121_SDF <- function(the_sdf,jump=0,freq=1/16)
{
    plot(the_sdf$freqs,the_sdf$sdf,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,3),yaxs="i",ylab="SDF",
         typ="l",axes=FALSE,
         main="Figure 121")
    if(jump > 0)
    {
        lines(rep(-freq,2),c(0,jump))
        lines(rep(freq,2),c(0,jump))
        points(freq*c(-1,1),rep(jump,2),pch=16)
    }
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(2,at=seq(0,3,1),las=2)
    box(bty="l")
}

### Figure 121, 2nd row

fig_121_SDF(list(freqs=freqs_sdf_two_sided,
                 sdf=ar_sdf_two_sided))

fig_121_SDF(NULL,0.5)

fig_121_SDF(list(freqs=c(-0.5,0.5),
                 sdf=rep(0.5,2)),
            0.25)

fig_121_SDF(list(freqs=freqs_sdf_two_sided,
                 sdf=0.5*ar_sdf_two_sided),
            0.25)

###

fig_121_ACVS <- function(the_acvs)
{
    plot(0:32,the_acvs,
         xlim=c(0,32),xlab=expression(tau),
         ylim=c(-1,1),ylab="ACVS",
         typ="b",cex=0.5,lwd=0.25,axes=FALSE,
         main="Figure 121")
    axis(1,at=seq(0,32,16))
    axis(2,at=seq(-1,1,1),las=2)
    box(bty="l")
}

phi <- 0.49
freq <- 1/16

### Figure 121, 3rd row

fig_121_ACVS(phi^(0:32))

fig_121_ACVS(cos(2*pi*freq*(0:32)))

fig_121_ACVS(c(0.5,rep(0,32)) + 0.5*cos(2*pi*freq*(0:32)))

fig_121_ACVS(0.5*phi^(0:32)+ 0.5*cos(2*pi*freq*(0:32)))

###

fig_121_TS <- function(the_ts)
{
    plot(0:32,the_ts,
         xlim=c(0,32),xlab=expression(italic(t)),
         ylim=c(-3,3),ylab="time series",
         typ="l",lwd=0.25,axes=FALSE,
         main="Figure 121")
    axis(1,at=seq(0,32,16))
    axis(2,at=seq(-3,3,3),las=2)
    box(bty="l")
}

phi <- 0.49
freq <- 1/16
set.seed(5)
ts_1 <- sim_ar_process(33,phi)
ts_2 <- sqrt(2) * cos(2*pi*freq*(0:32) + runif(1,-pi,pi))
ts_3 <- sqrt(2) * cos(2*pi*freq*(0:32) + runif(1,-pi,pi)) + rnorm(33,sd=sqrt(0.5))
ts_4 <- sqrt(2) * cos(2*pi*freq*(0:32) + runif(1,-pi,pi)) + sim_ar_process(33,phi,0.5)

### Figure 121, 3rd row

fig_121_TS(ts_1)

fig_121_TS(ts_2)

fig_121_TS(ts_3)

fig_121_TS(ts_4)

### Figure 123 ###

fig_123 <- function(delta_t,tag,k_max=10,f_max=3)
{
    sdf_cont <- function(f) 1.9*exp(-2*f^2) + exp(-6*(abs(f) - 1.25)^2)
    freqs_disc <- seq(0,1/(2*delta_t),length=1000)
    ys_disc <- rep(0,length(freqs_disc))
    for(j in 1:length(freqs_disc))
    {
        freq <- freqs_disc[j]
        ys_disc[j] <- sdf_cont(freq)
        for(k in 1:k_max)
            ys_disc[j] <- ys_disc[j] + sdf_cont(freq + k/delta_t) + sdf_cont(freq - k/delta_t)
    }
    plot(freqs_disc,ys_disc,
         xlim=c(0,f_max),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,2),ylab="SDFs",
         typ="l",lwd=1.5,axes=FALSE,
         main="Figure 123")
    fN <- tail(freqs_disc,1)
    abline(v=fN,lty="dotted")
    ## periodic extension of discrete parameter SDF
    freqs_ext <- c(freqs_disc + fN,freqs_disc + 2*fN)
    ys_ext <- c(rev(ys_disc),ys_disc)
    lines(freqs_ext, ys_ext, lwd=1.5, col="gray")
    ## continuous parameter SDF
    freqs_cont <- seq(0,f_max,length=1000)
    ys_cont <- sapply(freqs_cont,sdf_cont)
    lines(freqs_cont,ys_cont,lwd=0.5)
    axis(1,at=seq(0,f_max,1))
    axis(2,at=seq(0,2,1),las=2)
    text(f_max,1.75,tag,pos=2)
    box(bty="l")
}

### Figure 123, upper plot

fig_123(1/2,expression(Delta[t] == 1/2))

### Figure 123, lower plot

fig_123(1/4,expression(Delta[t] == 1/4))

### Figure 125 ###

fig_125_SDF <- function(the_sdf,tag)
{
    freqs <- seq(0,0.5,length=length(the_sdf)) 
    plot(freqs,dB(the_sdf),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="SDF  (dB)",
         typ="l",axes=FALSE,
         main=paste("Figure 125",tag,sep=""))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(0.5,10,tag,pos=2)
    box(bty="l")
}

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)
ar4_sdf <- ar_coeffs_to_sdf(ar4_coeffs,ar4_innov_var,N_pad=1024)$sdf
    
ar2_innov_var <- 0.0004
ar2_coeffs    <- c(-1.15,-0.96)
ar2_sdf <- ar_coeffs_to_sdf(ar2_coeffs,ar2_innov_var,N_pad=1024)$sdf

### Figure 125, left-hand column

fig_125_SDF(ar4_sdf,"(a)")

fig_125_SDF(ar4_sdf+ar2_sdf,"(c)")

###

fig_125_ACVS <- function(the_acvs,tag)
{
    taus <- 0:(length(the_acvs)-1)
    N_acvs <- length(the_acvs)
    plot(taus,the_acvs,
         xlim=c(0,N_acvs),xlab=expression(tau),
         ylim=c(-2,2),ylab="ACVS",
         typ="b",axes=FALSE,
         main=paste("Figure 125",tag,sep=""))
    axis(1,at=seq(0,N_acvs,20))
    axis(1,at=seq(0,N_acvs,10),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-2,2,2),las=2)
    axis(2,at=seq(-2,2,1),label=FALSE,tcl=-0.25)
    text(N_acvs,1.75,tag,pos=2)
    box(bty="l")
}

ar4_acvs <- ar_coeffs_to_acvs(ar4_coeffs,64,ar4_innov_var,FALSE)
ar2_acvs <- ar_coeffs_to_acvs(ar2_coeffs,62,ar4_innov_var,FALSE)

### Figure 125, right-hand column

fig_125_ACVS(ar4_acvs,"(b)")

fig_125_ACVS(ar4_acvs+ar2_acvs,"(d)")

### Figure 129a ###

fig_129a <- function(the_sdf,tag)
{
    freqs <- seq(0,0.5,length=length(the_sdf)) 
    plot(freqs,the_sdf,
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,10),yaxs="i",ylab="SDF  (dB)",
         typ="l",axes=FALSE,
         main=paste("Figure 129",tag,sep=""))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,10,5),las=2)
    axis(2,at=seq(0,10,1),label=FALSE,tcl=-0.25)
    text(0.5,9,tag,pos=2)
    box(bty="l")
}

ar1_innov_var <- 0.36
ar1_coeffs    <- 0.8
ar1_sdf <- ar_coeffs_to_sdf(ar1_coeffs,ar1_innov_var,N_pad=1024)$sdf

ar2_innov_var <- 0.48
ar2_coeffs    <- c(-0.8,-0.6)
ar2_sdf <- ar_coeffs_to_sdf(ar2_coeffs,ar2_innov_var,N_pad=1024)$sdf

wn_1_sdf <- rep(3,2)

wn_2_sdf <- rep(8,2)

### Figure 129a, top row

fig_129a(ar1_sdf,"(a)")

fig_129a(ar2_sdf,"(b)")

### Figure 129a, bottom row

fig_129a(wn_1_sdf,"(c)")

fig_129a(wn_2_sdf,"(d)")

### NOTE: code to recreate Figure 129b is not provided - to do so
###       would reveal the solution to Exercise [4.10]!
