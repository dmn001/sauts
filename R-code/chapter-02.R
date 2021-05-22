### R CODE FOR REPRODUCING CONTENT OF FIGURES IN CHAPTER 2 ...

ma_a <- scan("http://faculty.washington.edu/dbp/sauts/Data/map1.txt")
ma_b <- scan("http://faculty.washington.edu/dbp/sauts/Data/map2.txt")
ma_c <- scan("http://faculty.washington.edu/dbp/sauts/Data/map3.txt")
ma_d <- scan("http://faculty.washington.edu/dbp/sauts/Data/map4.txt")
ma_e <- scan("http://faculty.washington.edu/dbp/sauts/Data/mam1.txt")
ma_f <- scan("http://faculty.washington.edu/dbp/sauts/Data/mam2.txt")
ma_g <- scan("http://faculty.washington.edu/dbp/sauts/Data/mam3.txt")
ma_h <- scan("http://faculty.washington.edu/dbp/sauts/Data/mam4.txt")
ar2_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_1.txt")
ar2_2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_2.txt")
ar2_3 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_3.txt")
ar2_4 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_4.txt")
ar4_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_1.txt")
ar4_2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_2.txt")
ar4_3 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_3.txt")
ar4_4 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_4.txt")
rotor <- scan("http://faculty.washington.edu/dbp/sauts/Data/rotor.txt")
temp <- scan("http://faculty.washington.edu/dbp/sauts/Data/resistor.txt")
resistor_times <- temp[seq(1,1999,2)]
resistor_ts    <- temp[seq(2,2000,2)]

### functions used to compute content of figures in Chapter 2 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_arma_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")

###

generate_double_exponential <- function() -log(1-runif(1))/(if(runif(1) <= 0.5) sqrt(2) else -sqrt(2))

generate_discrete_rv <- function()
{
    rn <- runif(1)
    if(rn <= 0.02) -5 else if(rn <= 0.98) 0 else 5
}

generate_mismash_deviate <- function()
{
    rn <- runif(1)
    if(rn <= 0.25) rnorm(1) else if(rn <= 0.5) runif(1,-sqrt(3),sqrt(3)) else if(rn <= 0.75) generate_double_exponential() else generate_discrete_rv()
}

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 31 ###

fig_31 <- function(ts,tag="(a)")
{
    xs <- 0:99
    plot(xs,ts,
         xlim=c(0,100),xlab=expression(italic(t)),
         ylim=c(-5,5),ylab="white noise",
         typ="l",axes=FALSE,
         main=paste("Figure 31",tag,sep=""))
    axis(1,at=seq(0,100,50))
    axis(1,at=seq(0,100,25),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-5,5,5),las=2)
    axis(2,at=seq(-5,5,1),label=FALSE,tcl=-0.25)
    text(90,5,tag,pos=4)
    box(bty="l")
}

set.seed(7)
wn_1 <- rnorm(100)
wn_2 <- runif(100,-sqrt(3),sqrt(3))
wn_3 <- replicate(100,generate_double_exponential())
wn_4 <- replicate(100,generate_discrete_rv())
wn_5 <- replicate(100,generate_mismash_deviate())
wn_6 <- c(wn_1[1:25],wn_2[1:25],wn_3[1:25],wn_4[1:25])

### Figure 31(a)

fig_31(wn_1)

### Figure 31(b)

fig_31(wn_2,"(b)")

### Figure 31(c)

fig_31(wn_3,"(c)")

### Figure 31(d)

fig_31(wn_4,"(d)")

### Figure 31(e)

fig_31(wn_5,"(e)")

### Figure 31(f)

fig_31(wn_6,"(f)")

### Figure 33 ###

fig_33 <- function(ts,tag="(a)")
{
    xs <- 0:127
    plot(xs,ts,
         xlim=c(0,128),xlab=expression(italic(t)),
         ylim=c(-6,6),ylab="MA(1) series",
         typ="l",axes=FALSE,
         main=paste("Figure 33",tag,sep=""))
    axis(1,at=seq(0,128,64))
    axis(1,at=seq(0,128,32),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-6,6,6),las=2)
    axis(2,at=seq(-6,6,1),label=FALSE,tcl=-0.25)
    text(120,6,tag,pos=4)
    box(bty="l")
}

### Figure 33(a)

fig_33(ma_a)

### Figure 33(b)

fig_33(ma_b,"(b)")

### Figure 33(c)

fig_33(ma_c,"(c)")

### Figure 33(d)

fig_33(ma_d,"(d)")

### Figure 33(e)

fig_33(ma_e,"(e)")

### Figure 33(f)

fig_33(ma_f,"(f)")

### Figure 33(g)

fig_33(ma_g,"(g)")

### Figure 33(h)

fig_33(ma_h,"(h)")

### aside: generation of additional MA(1) series ...

temp <- rnorm(129)
fig_33(temp[-1]-temp[-129],"(w)")  # theta =  1
fig_33(temp[-1]+temp[-129],"(x)")  # theta = -1

fig_33(sim_arma_process(128,theta=1),"(y)")
fig_33(sim_arma_process(128,theta=-1),"(z)")

### Figure 34 ###

fig_34 <- function(ts,tag="(a)",two_or_four="2")
{
    xs <- 0:1023
    plot(xs,ts,
         xlim=c(0,1024),xlab=expression(italic(t)),
         ylim=c(-5,5),ylab=paste("AR(",two_or_four,") series",sep=""),
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 34",tag,sep=""))
    axis(1,at=seq(0,1024,512))
    axis(1,at=seq(0,1024,256),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-5,5,5),las=2)
    axis(2,at=seq(-5,5,1),label=FALSE,tcl=-0.25)
    text(1000,5,tag,pos=4)
    box(bty="l")
}

### Figure 34(a)

fig_34(ar2_1)

### Figure 34(b)

fig_34(ar2_2,"(b)")

### Figure 34(c)

fig_34(ar2_3,"(c)")

### Figure 34(d)

fig_34(ar2_4,"(d)")

### Figure 34(e)

fig_34(ar4_1,"(e)","4")

### Figure 34(f)

fig_34(ar4_2,"(f)","4")

### Figure 34(g)

fig_34(ar4_3,"(g)","4")

### Figure 34(h)

fig_34(ar4_4,"(h)","4")

### aside: generation of additional AR(2) and AR(4) series ...

ar2_innov_var <- 1.0
ar2_coeffs <- c(0.75, -0.5)

fig_34(sim_ar_process(1024,ar2_coeffs,ar2_innov_var),"(y)")

ar4_innov_var <- 0.002
ar4_coeffs <- c(2.7607, -3.8106, 2.6535, -0.9238)

fig_34(sim_ar_process(1024,ar4_coeffs,ar4_innov_var),"(z)","4")

### Figure 36 ###

fig_36 <- function(ts,tag="(a)")
{
    xs <- 0:99
    plot(xs,ts,
         xlim=c(0,100),xaxs="i",xlab=expression(italic(t)),
         ylim=c(-3,3),ylab="harmonic series",
         typ="l",axes=FALSE,
         main=paste("Figure 36",tag,sep=""))
    axis(1,at=seq(0,100,50))
    axis(1,at=seq(0,100,25),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-3,3,3),las=2)
    axis(2,at=seq(-3,3,1),label=FALSE,tcl=-0.25)
    text(95,3,tag,pos=4)
    box(bty="l")
}

set.seed(42)
xs <- 0:99
f_1 <- 1/20
harm_1 <- rnorm(1)*cos(2*pi*f_1*xs) + rnorm(1)*sin(2*pi*f_1*xs) 
harm_3 <- rnorm(1)*cos(2*pi*f_1*xs) + rnorm(1)*sin(2*pi*f_1*xs) 
harm_5 <- rnorm(1)*cos(2*pi*f_1*xs) + rnorm(1)*sin(2*pi*f_1*xs) 
harm_2 <- sqrt(2)*cos(2*pi*f_1*xs + runif(1,-pi,pi)) 
harm_4 <- sqrt(2)*cos(2*pi*f_1*xs + runif(1,-pi,pi)) 
harm_6 <- sqrt(2)*cos(2*pi*f_1*xs + runif(1,-pi,pi)) 

### Figure 36(a)

fig_36(harm_1)

### Figure 36(b)

fig_36(harm_2,"(b)")

### Figure 36(c)

fig_36(harm_3,"(c)")

### Figure 36(d)

fig_36(harm_4,"(d)")

### Figure 36(e)

fig_36(harm_5,"(e)")

### Figure 36(f)

fig_36(harm_6,"(f)")

### Figure 39 ###

fig_39 <- function(ts,y_lims=c(200,1000),y_lab="time series",tag="(a)")
{
    xs <- (221-length(ts)):220
    plot(xs,ts,
         xlim=c(0,222),xlab="bin number",
         ylim=y_lims,ylab=y_lab,
         typ="l",axes=FALSE,
         main=paste("Figure 39",tag,sep=""))
    axis(1,at=seq(0,222,111))
    axis(2,at=c(y_lims[1],mean(y_lims),y_lims[2]),las=2)
    text(200,y_lims[2],tag,pos=4)
    box(bty="l")
}

### Figure 39(a)

fig_39(rotor)

### Figure 39(b)

fig_39(resid(lm(rotor ~ I(0:220))),c(-40,40),"residuals","(b)")

### Figure 39(c)

fig_39(diff(rotor),y_lims=c(-80,80),y_lab="first difference","(c)")

### Figure 40 ###

fig_40 <- function(times,ts)
{
    plot(times,ts,
         xlim=c(0,2190),xlab="days  (from 1 January 1980)",
         ylim=c(27.8,28.2),ylab="time series",
         typ="l",axes=FALSE,
         main="Figure 40")
    axis(1,at=seq(0,2190,365))
    axis(2,at=seq(27.8,28.2,0.2),las=2)
    box(bty="l")
}

###

fig_40(resistor_times,resistor_ts)
