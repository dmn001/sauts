### R CODE FOR REPRODUCING CONTENT OF FIGURES IN CHAPTER 1 ...

wind_speed <- scan("http://faculty.washington.edu/dbp/sauts/Data/wind_speed_128.txt")
atomic_clock <- scan("http://faculty.washington.edu/dbp/sauts/Data/atomic_clock_128.txt")
Willamette_River <- scan("http://faculty.washington.edu/dbp/sauts/Data/Willamette_River_128.txt")
ocean_noise <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_noise_128.txt")
st_paul <- scan("http://faculty.washington.edu/dbp/sauts/Data/st_paul.txt")
sim_wind_speed <- scan("http://faculty.washington.edu/dbp/sauts/Data/sim_wind_speed.txt")
sim_atomic_clock <- scan("http://faculty.washington.edu/dbp/sauts/Data/sim_atomic_clock.txt")
sim_Willamette_River <- scan("http://faculty.washington.edu/dbp/sauts/Data/sim_Willamette_River.txt")
sim_ocean_noise <- scan("http://faculty.washington.edu/dbp/sauts/Data/sim_ocean_noise.txt")

### functions used to compute content of figures in Chapter 1 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")

###

compute_tspec <- function(alpha,beta=1,N=128) beta/(1 + alpha^2 - 2*alpha*cos(2*pi*(1:(N/2))/N))

###

compute_Sj <- function(j,xs)
{
    N <- length(xs)
    times <- 0:(N-1)
    f_j <- j/N
    temp <- (sum(xs*cos(2*pi*f_j*times)))^2
    if(j == N/2) temp/N^2
    else
        2*(temp + (sum(xs*sin(2*pi*f_j*times)))^2)/N^2
}

###

compute_tacs <- function(tau,tspec,N=2*length(tspec)) sum(tspec * cos(2*pi*(1:length(tspec))*tau/N))/sum(tspec)

###

gen_via_Equation_8a <- function(tspec,ts_mean,N=2*length(tspec))
{
    f_j <- (1:length(tspec))/N
    A_j <- rnorm(length(tspec))*sqrt(tspec)
    B_j <- rnorm(length(tspec))*sqrt(tspec)
    return(ts_mean + sapply(0:(N-1),function(t) sum(A_j*cos(2*pi*f_j*t) + B_j*sin(2*pi*f_j*t)))) 
}

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 2 ###

fig_2 <- function(ts,y_lim=c(-6,4),tag="(a)")
{
    xs <- 0:127
    plot(xs,ts,
         xlab=expression(paste(italic(t),"  (time index)")),
         ylim=y_lim,ylab=expression(paste(italic(x[t]),"  (time series)")),
         typ="l",axes=FALSE,
         main=paste("Figure 2",tag,sep=""))
    axis(1,at=seq(0,128,64))
    axis(1,at=seq(0,128,32),label=FALSE,tcl=-0.25)
    axis(2,at=c(y_lim[1],mean(y_lim),y_lim[2]),las=2)
    text(120,y_lim[2],tag,pos=4)
    box(bty="l")
}

### Figure 2(a)

fig_2(wind_speed)

### Figure 2(b)

fig_2(atomic_clock,c(-15,15),"(b)")

### Figure 2(c)

fig_2(Willamette_River,c(8,12),"(c)")

### Figure 2(d)

fig_2(ocean_noise,c(-7,7),"(d)")

### Figure 3 ###

fig_3 <- function(ts,xy_lim=c(-6,4),tag="(a)",left_p=FALSE,y_lab=expression(italic(x[t+1])),lag=1)
{
    N <- length(ts)
    plot(ts[-((N-lag+1):N)],ts[-(1:lag)],
         xlim=xy_lim,xlab=expression(italic(x[t])),
         ylim=xy_lim,ylab=y_lab,
         pch=16,cex=0.5,axes=FALSE,
         main=paste("Figure 3",tag,sep=""))
    ats_xy <- c(xy_lim[1],mean(xy_lim),xy_lim[2])
    axis(1,at=ats_xy)
    axis(2,at=ats_xy,las=2)
    text(if(left_p) xy_lim[2]-0.1*diff(range(xy_lim)) else xy_lim[1],xy_lim[2],tag,pos=4)
    box(bty="l")
}

### Figure 3(a)

fig_3(wind_speed)

### Figure 3(b)

fig_3(atomic_clock,c(-15,15),"(b)",TRUE)

### Figure 3(c)

fig_3(Willamette_River,c(8,12),"(c)")

### Figure 3(d)

fig_3(ocean_noise,c(-7,7),"(d)",TRUE)

### Figure 4 ###

fig_4 <- function(ts,tag="(a)")
{
    taus <- 0:32
    the_acs <- acvs(ts,acs=TRUE)$acvs[1:33]
    plot(taus,the_acs,
         xlab=expression(tau),
         ylim=c(-1,1),ylab="ACS",
         typ="b",axes=FALSE,
         main=paste("Figure 4",tag,sep=""))
    abline(h=0)
    axis(1,at=seq(0,32,8))
    axis(1,at=seq(0,32,1),label=FALSE,tcl=-0.25)
    axis(2,at=c(-1,0,1),las=2)
    text(30,1,tag,pos=4)
    box(bty="l")
}

### Figure 4(a)

fig_4(wind_speed)

### Figure 4(b)

fig_4(atomic_clock,"(b)")

### Figure 4(c)

fig_4(Willamette_River,"(c)")

### Figure 4(d)

fig_4(ocean_noise,"(d)")

### Figure 6 ###

fig_6_top <- function(ts)
{
    xs <- seq(1820,1983-1/12,1/12)
    plot(xs,ts,
         xlim=c(1820,1844),xaxs="i",xlab=expression(paste(italic(t),"  (years)")),
         ylim=c(-30,30),ylab=expression(italic(x[t])),
         typ="l",axes=FALSE,
         main="Figure 6(a)")
    axis(1,at=seq(1820,1844,8))
    axis(1,at=seq(1820,1844,1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-30,30,30),las=2)
    axis(2,at=seq(-30,30,10),label=FALSE,tcl=-0.25)
    text(1842.2,30,"(a)",pos=4)
    box(bty="l")
}

### Figure 6(a)

fig_6_top(st_paul)

###

fig_6_bot <- function(ts,lag=6,tag="(b)",y_lab=expression(italic(x[t+6])))
{
    xy_lim=c(-30,30)
    N <- length(ts)
    plot(ts[-((N-lag+1):N)],ts[-(1:lag)],
         xlim=xy_lim,xlab=expression(italic(x[t])),
         ylim=xy_lim,ylab=y_lab,
         pch=".",axes=FALSE,
         main=paste("Figure 6",tag,sep=""))
    axis(1,at=seq(-30,30,30))
    axis(1,at=seq(-30,30,10),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-30,30,30),las=2)
    axis(2,at=seq(-30,30,10),label=FALSE,tcl=-0.25)
    text(24,30,tag,pos=4)
    box(bty="l")
}

### Figure 6(b)

fig_6_bot(st_paul)

### Figure 6(c)

fig_6_bot(st_paul,9,"(c)",expression(italic(x[t+9])))

### Figure 7 ###

set.seed(2)
times <- 0:127
sinusoids <- matrix(nrow=128,ncol=60)
for(k in 1:10)
  {
    j <- 2*k - 1
    sinusoids[,((k-1)*6+1):(k*6)] <- cbind(cos(2*pi*times*j/128),
                                           sin(2*pi*times*j/128),
                                           rnorm(1)*cos(2*pi*times*j/128),
                                           rnorm(1)*sin(2*pi*times*j/128),
                                           rnorm(1)*cos(2*pi*times*j/128),
                                           rnorm(1)*sin(2*pi*times*j/128))
  }
sinusoids_sum <- matrix(nrow=128,ncol=3)
cols_to_sum <- c(seq(1, 55, 6),seq(2, 56, 6))
sinusoids_sum[,1] <- rowSums(sinusoids[,cols_to_sum])
cols_to_sum <- c(seq(3, 57, 6),seq(4, 58, 6))
sinusoids_sum[,2] <- rowSums(sinusoids[,cols_to_sum])
cols_to_sum <- c(seq(5, 59, 6),seq(6, 60, 6))
sinusoids_sum[,3] <- rowSums(sinusoids[,cols_to_sum])

fig_7 <- function(ts_1,ts_2,y_lab,which_col,y_lim=c(-2.4,2.4))
{
    xs <- 0:127
    plot(xs,ts_1,
         xlim=c(0,128),xaxs="i",xlab=expression(italic(t)),
         ylim=y_lim,yaxs="i",ylab=y_lab,
         typ="l",axes=FALSE,
         main=if(is.expression(y_lab))
                  paste("Figure 7, ",which_col," Column, Time Series",sep="")
              else 
                  paste("Figure 7, ",which_col," Column, ",y_lab,sep=""))
    if(!is.null(ts_2)) lines(xs,ts_2,col="darkgray")
    axis(1,at=c(0,128))
    box(bty="l")
}

### Figure 7, left-hand column

fig_7(sinusoids[,1],sinusoids[,2],"1/128","Left-hand")
fig_7(sinusoids[,7],sinusoids[,8],"3/128","Left-hand")
fig_7(sinusoids[,13],sinusoids[,14],"5/128","Left-hand")
fig_7(sinusoids[,19],sinusoids[,20],"7/128","Left-hand")
fig_7(sinusoids[,25],sinusoids[,26],"9/128","Left-hand")
fig_7(sinusoids[,31],sinusoids[,32],"11/128","Left-hand")
fig_7(sinusoids[,37],sinusoids[,38],"13/128","Left-hand")
fig_7(sinusoids[,43],sinusoids[,44],"15/128","Left-hand")
fig_7(sinusoids[,49],sinusoids[,50],"17/128","Left-hand")
fig_7(sinusoids[,55],sinusoids[,56],"19/128","Left-hand")
fig_7(sinusoids_sum[,1],NULL,expression(italic(x[t])),"Left-hand",y_lim=c(-14,14))

### Figure 7, middle column

fig_7(sinusoids[,3],sinusoids[,4],"1/128","Middle")
fig_7(sinusoids[,9],sinusoids[,10],"3/128","Middle")
fig_7(sinusoids[,15],sinusoids[,16],"5/128","Middle")
fig_7(sinusoids[,21],sinusoids[,22],"7/128","Middle")
fig_7(sinusoids[,27],sinusoids[,28],"9/128","Middle")
fig_7(sinusoids[,33],sinusoids[,34],"11/128","Middle")
fig_7(sinusoids[,39],sinusoids[,40],"13/128","Middle")
fig_7(sinusoids[,45],sinusoids[,46],"15/128","Middle")
fig_7(sinusoids[,51],sinusoids[,52],"17/128","Middle")
fig_7(sinusoids[,57],sinusoids[,58],"19/128","Middle")
fig_7(sinusoids_sum[,2],NULL,expression(italic(x[t])),"Middle",y_lim=c(-14,14))

### Figure 7, right-hand column

fig_7(sinusoids[,5],sinusoids[,6],"1/128","Right-hand")
fig_7(sinusoids[,11],sinusoids[,12],"3/128","Right-hand")
fig_7(sinusoids[,17],sinusoids[,18],"5/128","Right-hand")
fig_7(sinusoids[,23],sinusoids[,24],"7/128","Right-hand")
fig_7(sinusoids[,29],sinusoids[,30],"9/128","Right-hand")
fig_7(sinusoids[,35],sinusoids[,36],"11/128","Right-hand")
fig_7(sinusoids[,41],sinusoids[,42],"13/128","Right-hand")
fig_7(sinusoids[,47],sinusoids[,48],"15/128","Right-hand")
fig_7(sinusoids[,53],sinusoids[,54],"17/128","Right-hand")
fig_7(sinusoids[,59],sinusoids[,60],"19/128","Right-hand")
fig_7(sinusoids_sum[,3],NULL,expression(italic(x[t])),"Right-hand",y_lim=c(-14,14))

### Figure 10a ###

acvs_a <- acvs(wind_speed)$acvs
var_a <- acvs_a[1] 
alpha_a <- acvs_a[2]/var_a
beta_a  <- var_a/sum(compute_tspec(alpha_a))
tspec_a   <- compute_tspec(alpha_a,beta_a)
all.equal(sum(tspec_a),var_a)  # TRUE

acvs_b <- acvs(atomic_clock)$acvs
var_b <- acvs_b[1] 
alpha_b <- acvs_b[2]/var_b
beta_b  <- var_b/sum(compute_tspec(alpha_b))
tspec_b   <- compute_tspec(alpha_b,beta_b)
all.equal(sum(tspec_b),var_b)  # TRUE

acvs_c <- acvs(Willamette_River)$acvs
var_c <- acvs_c[1] 
espec_c <- sapply(1:64,compute_Sj,Willamette_River)
tspec_c      <- rep(espec_c[11],64)
tspec_c[-11] <- (var_c-espec_c[11])/63
all.equal(sum(tspec_c),var_c)  # TRUE

alpha_d <- 0
var_d <- acvs(ocean_noise)$acvs[1]
beta_d  <- var_d/sum(compute_tspec(alpha_d))
tspec_d   <- compute_tspec(alpha_d,beta_d)
all.equal(sum(tspec_d),var_d)  # TRUE

fig_10a <- function(tspec,y_lim,tag)
{
    N <- length(tspec)
    freqs <- (1:N)/(2*N)
    plot(freqs,tspec,
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f[j])),
         ylim=y_lim,yaxs="i",ylab="spectrum",
         typ="l",axes=FALSE,
         main=paste("Figure 10a",tag,sep=""))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=y_lim,las=2)
    text(0.46,y_lim[2]*0.95,tag,pos=4)
    box(bty="l")
}

### Figure 10a(a)

fig_10a(tspec_a,c(0,0.9),"(a)")

### Figure 10a(b)

fig_10a(tspec_b,c(0,2),"(b)")

### Figure 10a(c)

fig_10a(tspec_c,c(0,0.4),"(c)")

### Figure 10a(d)

fig_10a(tspec_d,c(0,0.6),"(d)")

### Figure 10b ###

tacs_a <- sapply(0:32, compute_tacs, tspec_a)
tacs_b <- sapply(0:32, compute_tacs, tspec_b)
tacs_c <- sapply(0:32, compute_tacs, tspec_c)
tacs_d <- sapply(0:32, compute_tacs, tspec_d)

fig_10b <- function(the_acs,tag="(a)")
{
    taus <- 0:32
    plot(taus,the_acs,
         xlab=expression(tau),
         ylim=c(-1,1),ylab="ACS",
         typ="b",lwd=0.25,cex=0.5,axes=FALSE,
         main=paste("Figure 10b",tag,sep=""))
    abline(h=0)
    axis(1,at=seq(0,32,8))
    axis(1,at=seq(0,32,1),label=FALSE,tcl=-0.25)
    axis(2,at=c(-1,0,1),las=2)
    text(30,1,tag,pos=4)
    box(bty="l")
}

### Figure 10b(a)

fig_10b(tacs_a)

### Figure 10b(b)

fig_10b(tacs_b,"(b)")

### Figure 10b(c)

fig_10b(tacs_c,"(c)")

### Figure 10b(d)

fig_10b(tacs_d,"(d)")

### Figure 11 ###

fig_11 <- function(ts,y_lim=c(-6,4),tag="(a)")
{
    xs <- 0:127
    plot(xs,ts,
         xlab=expression(italic(t)),
         ylim=y_lim,ylab="simulated series",
         typ="l",axes=FALSE,
         main=paste("Figure 11",tag,sep=""))
    axis(1,at=seq(0,128,64))
    axis(1,at=seq(0,128,32),label=FALSE,tcl=-0.25)
    axis(2,at=c(y_lim[1],mean(y_lim),y_lim[2]),las=2)
    text(120,y_lim[2],tag,pos=4)
    box(bty="l")
}

### NOTE: see Figure 10a for computation of tspec_a, tspec_b, tspec_c & tspec_d

set.seed(42)

### Figure 11(a)

fig_11(gen_via_Equation_8a(tspec_a,mean(wind_speed)))

### Figure 11(b)

fig_11(gen_via_Equation_8a(tspec_b,mean(atomic_clock)),c(-15,15),"(b)")

### Figure 11(c)

fig_11(gen_via_Equation_8a(tspec_c,mean(Willamette_River)),c(8,12),"(c)")

### Figure 11(d)

fig_11(gen_via_Equation_8a(tspec_d,mean(ocean_noise)),c(-7,7),"(d)")

### Figure 12 ###

espec_a <- sapply(1:64,compute_Sj,wind_speed)
espec_b <- sapply(1:64,compute_Sj,atomic_clock)
espec_c <- sapply(1:64,compute_Sj,Willamette_River)
espec_d <- sapply(1:64,compute_Sj,ocean_noise)

fig_12 <- function(tspec,espec,y_lim,tag,right_p=TRUE)
{
    N <- length(tspec)
    freqs <- (1:N)/(2*N)
    plot(freqs,espec,
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f[j])),
         ylim=y_lim,yaxs="i",ylab="spectra",
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 12",tag,sep=""))
    lines(freqs,tspec)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=y_lim,las=2)
    text(if(right_p) 0.46 else 0.04,y_lim[2]*0.95,tag,pos=4)
    box(bty="l")
}

### NOTE: see Figure 10a for computation of tspec_a, tspec_b, tspec_c & tspec_d

### Figure 12(a)

fig_12(tspec_a,espec_a,c(0,0.9),"(a)")

### Figure 12(b)

fig_12(tspec_b,espec_b,c(0,2),"(b)",FALSE)

### Figure 12(c)

fig_12(tspec_c,espec_c,c(0,0.4),"(c)")

### Figure 12(d)

fig_12(tspec_d,espec_d,c(0,0.6),"(d)")

### Figure 14 ###

fig_14 <- function(tspec,espec,y_lim,tag)
{
    N <- length(tspec)
    freqs <- (1:N)/(2*N)
    plot(freqs,dB(espec),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f[j])),
         ylim=y_lim,yaxs="i",ylab="spectra  (dB)",
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 14",tag,sep=""))
    lines(freqs,dB(tspec))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(y_lim[1],y_lim[2],20),las=2)
    axis(2,at=seq(y_lim[1],y_lim[2],10),label=FALSE,tcl=-0.25)
    text(0.46,dB(y_lim[2])-3,tag,pos=4)
    box(bty="l")
}

### NOTE: see Figure 10a for computation of tspec_a, tspec_b, tspec_c & tspec_d;
###       see Figure 12 for computation of espec_a, espec_b, espec_c & espec_d

### Figure 14(a)

fig_14(tspec_a,espec_a,c(-50,10),"(a)")

### Figure 14(b)

fig_14(tspec_b,espec_b,c(-50,10),"(b)")

### Figure 14(c)

fig_14(tspec_c,espec_c,c(-60,0),"(c)")

### Figure 14(d)

fig_14(tspec_d,espec_d,c(-60,0),"(d)")
