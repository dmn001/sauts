### R CODE FOR REPRODUCING CONTENT OF FIGURES AND TABLES IN CHAPTER 6 ...

wind_speed <- scan("http://faculty.washington.edu/dbp/sauts/Data/wind_speed_128.txt")
ar2_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_1.txt")
ar2_2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_2.txt")
ar2_3 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_3.txt")
ar2_4 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_4.txt")
ar4_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_1.txt")
ar4_2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_2.txt")
ar4_3 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_3.txt")
ar4_4 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_4.txt")
earth_20 <- scan("http://faculty.washington.edu/dbp/sauts/Data/earth_20.txt")
ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_wave.txt")
chaotic_beam <- scan("http://faculty.washington.edu/dbp/sauts/Data/chaotic_beam.txt")
ocean_noise <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_noise_128.txt")

### functions used to compute content of figures in Chapter 6 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/circular_shift.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/cosine_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/direct_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_DCTII.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_lag_window_sdf_estimator.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_shp.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_shp_squared.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/fejer_kernel.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/hanning_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_even.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/next_power_of_2.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/rectangular_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/spec_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")

###

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 168 ###

fig_168_top_row <- function(the_acvs,tag)
{
    N <- length(the_acvs)
    taus <- 0:(N-1)
    plot(taus,the_acvs,
         xlim=c(0,N),xlab=expression(tau),
         ylim=c(-2,2),ylab="ACVS",
         typ="b",lwd=0.25,cex=0.5,axes=FALSE,
         main="Figure 168")
    axis(1,at=seq(0,60,20))
    axis(1,at=seq(0,60,10),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-2,2,2),las=2)
    axis(2,at=seq(-2,2,1),label=FALSE,tcl=-0.25)
    text(60,1.8,tag,pos=2)
    box(bty="l")
}

fig_168_bot_rows <- function(biased,unbiased,y_lab)
{
    max_lag <- length(biased)-1
    taus <- 0:max_lag
    plot(taus,dB(unbiased),
         xlim=c(0,max_lag+1),xlab=expression(tau),
         ylim=c(-80,20),ylab=y_lab,
         typ="l",lwd=0.5,col="gray40",axes=FALSE,
         main="Figure 168")
    lines(taus,dB(biased))
    axis(1,at=seq(0,60,20))
    axis(1,at=seq(0,60,10),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,40,20),las=2)
    axis(2,at=seq(-80,40,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

ar2_acvs <- ar_coeffs_to_acvs(ar2_coeffs,63,ar2_innov_var,FALSE)
ar4_acvs <- ar_coeffs_to_acvs(ar4_coeffs,63,ar4_innov_var,FALSE)

b_to_u <- 64/(64:1)

### NOTE: evaluation of the following R code is time consuming
###       (particularly the two lines involving ev_shp_squared,
###       each of which took 45 minutes to execute on a 2017-vintage
###       MacBook Pro):
###
###         ev_shp_ar2 <- sapply(0:63,ev_shp,64,ar2_acvs)
###         ev_shp_squared_ar2 <- sapply(0:63,ev_shp_squared,64,ar2_acvs)
###         ev_shp_ar4 <- sapply(0:63,ev_shp,64,ar4_acvs)
###         ev_shp_squared_ar4 <- sapply(0:63,ev_shp_squared,64,ar4_acvs)
###
###       Evaluation of the following four load forms alleviates
###       having to recreate ev_shp_ar2 etc.

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ev_shp_ar2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ev_shp_squared_ar2.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ev_shp_ar4.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ev_shp_squared_ar4.Rdata"))

### Figure 168, plots in left-hand column from top to bottom

fig_168_top_row(ar2_acvs,"AR(2)")
fig_168_bot_rows((ev_shp_ar2-ar2_acvs)^2,(b_to_u*ev_shp_ar2-ar2_acvs)^2,"squared bias  (dB)")
fig_168_bot_rows(ev_shp_squared_ar2-ev_shp_ar2^2,b_to_u^2*(ev_shp_squared_ar2-ev_shp_ar2^2),"variance  (dB)")
fig_168_bot_rows(ev_shp_squared_ar2-ev_shp_ar2^2+(ev_shp_ar2-ar2_acvs)^2,b_to_u^2*(ev_shp_squared_ar2-ev_shp_ar2^2)+(b_to_u*ev_shp_ar2-ar2_acvs)^2,"MSE  (dB)")

### Figure 168, plots in right-hand column from top to bottom

fig_168_top_row(ar4_acvs,"AR(2)")
fig_168_bot_rows((ev_shp_ar4-ar4_acvs)^2,(b_to_u*ev_shp_ar4-ar4_acvs)^2,"squared bias  (dB)")
fig_168_bot_rows(ev_shp_squared_ar4-ev_shp_ar4^2,b_to_u^2*(ev_shp_squared_ar4-ev_shp_ar4^2),"variance  (dB)")
fig_168_bot_rows(ev_shp_squared_ar4-ev_shp_ar4^2+(ev_shp_ar4-ar4_acvs)^2,b_to_u^2*(ev_shp_squared_ar4-ev_shp_ar4^2)+(b_to_u*ev_shp_ar4-ar4_acvs)^2,"MSE  (dB)")

### Figure 169 ###

fig_169 <- function(ts)
{
    temp <- acvs(ts)
    taus <- temp$lags
    acvs_biased <- temp$acvs 
    acvs_unbiased <- acvs(ts,unbiased=TRUE)$acvs
    plot(taus,acvs_biased,
         xlim=c(0,length(ts)),xlab=expression(paste(tau," (in 0.025 sec)")),
         ylim=c(-4,4),ylab="ACVS",
         typ="l",axes=FALSE,
         main="Figure 169")
    lines(taus,acvs_unbiased,lwd=0.5,col="gray40")
    abline(h=0,lty="dashed")
    axis(1,at=seq(0,128,32))
    axis(2,at=seq(-4,4,2),las=2)
    axis(2,at=seq(-4,4,1),label=FALSE,tcl=-0.25)
    axis(4,at=seq(-4,4,2),label=FALSE)
    axis(4,at=seq(-4,4,1),label=FALSE,tcl=-0.25)
    box(bty="u")
}

### Figure 169

fig_169(wind_speed)

### Figure 172 ###

fig_172 <- function(ts,coeffs,innov_var,y_ats,tag)
{
    the_pgram <- pgram(ts,center=FALSE)
    plot(the_pgram$freqs,the_pgram$sdfe,
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,y_ats[length(y_ats)]),yaxs="i",ylab=paste("AR(",length(coeffs),") spectra",sep=""),
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main=paste("Figure 172",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,the_ar_spec$sdf)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=y_ats,las=2)
    text(x=0.5,y=0.95*y_ats[length(y_ats)],tag,pos=2)
    box(bty="l")
}

### Figure 172, top row of plots

fig_172(ar2_1,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(a)")
fig_172(ar2_2,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(b)")

### Figure 172, 2nd row of plots

fig_172(ar2_3,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(c)")
fig_172(ar2_4,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(d)")

### Figure 172, 3rd row of plots

fig_172(ar4_1,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(e)")
fig_172(ar4_2,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(f)")

### Figure 172, bottom row of plots

fig_172(ar4_3,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(g)")
fig_172(ar4_4,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(h)")

### Figure 173 ###

fig_173 <- function(ts,coeffs,innov_var,y_ats,tag)
{
    the_pgram <- pgram(ts,center=FALSE)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab=paste("AR(",length(coeffs),") spectra  (dB)",sep=""),
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main=paste("Figure 173",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    if(length(coeffs) == 4)
    {
        N <- length(ts)
        temp <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,N-1,innov_var,FALSE))
        lines(temp$freqs, dB(temp$sdf_ev), lwd=0.5)
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=10,tag,pos=2)
    box(bty="l")
}

### Figure 173, top row of plots

fig_173(ar2_1,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(a)")
fig_173(ar2_2,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(b)")

### Figure 173, 2nd row of plots

fig_173(ar2_3,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(c)")
fig_173(ar2_4,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(d)")

### Figure 173, 3rd row of plots

fig_173(ar4_1,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(e)")
fig_173(ar4_2,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(f)")

### Figure 173, bottom row of plots

fig_173(ar4_3,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(g)")
fig_173(ar4_4,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(h)")

### Figure 176 ###

fig_176 <- function(N,right_p=FALSE,tag=NULL)
{
    the_kernel <- fejer_kernel(N)
    plot(the_kernel$freqs,if(!right_p) dB(the_kernel$kernel) else the_kernel$kernel,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=if(!right_p) c(-40,20) else c(0,N),yaxs="i",ylab="spectral window",
         typ="l",lwd=0.25,axes=FALSE,
         main="Figure 176")
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=if(!right_p) seq(-40,20,20) else seq(0,N,N/2),las=2)
    if(!right_p)
    {
        axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
        text(x=0.5,y=10,tag,pos=2)
    }
    box(bty="l")
}

### Figure 176, left-hand column of plots

fig_176(4,tag=expression(italic(N==4)))
fig_176(16,tag=expression(italic(N==16)))
fig_176(64,tag=expression(italic(N==64)))

### Figure 176, right-hand column of plots

fig_176(4,tag=expression(italic(N==4)),right_p=TRUE)
fig_176(16,tag=expression(italic(N==16)),right_p=TRUE)
fig_176(64,tag=expression(italic(N==64)),right_p=TRUE)

### Figure 177 ###

fig_177 <- function(N,coeffs,innov_var,tag_1,tag_2,tag_3=NULL)
{
    temp <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,N-1,innov_var,FALSE),N_pad=1024)
    plot(temp$freqs,dB(temp$sdf_ev),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-10,10),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 177",tag_1,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-10,10,10),las=2)
    axis(2,at=seq(-10,10,2),label=FALSE,tcl=-0.25)
    text(x=0.5,y=9,tag_1,pos=2)
    text(x=0.25,y=-8,tag_2,pos=1)
    text(x=0.25,y=9,tag_3,pos=1)
    box(bty="l")
}

### Figure 177, left-hand and right-hand plots

fig_177(16,ar2_coeffs,ar2_innov_var,"(a)",expression(italic(N==16)),"AR(2)")
fig_177(64,ar2_coeffs,ar2_innov_var,"(b)",expression(italic(N==64)))

### Figure 178 ###

fig_178 <- function(N,coeffs,innov_var,tag_1,tag_2,tag_3=NULL,vlines=NULL)
{
    temp <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,N-1,innov_var,FALSE),N_pad=2048)
    plot(temp$freqs,dB(temp$sdf_ev),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 178",tag_1,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    abline(v=vlines, lty="dotted")
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=15,tag_1,pos=2)
    text(x=0.25,y=-50,tag_2,pos=1)
    text(x=0.25,y=20,tag_3,pos=1)
    box(bty="l")
}

### Figure 178, top row of plots

fig_178(16,ar4_coeffs,ar4_innov_var,"(a)",expression(italic(N==16)),"AR(4)")
fig_178(64,ar4_coeffs,ar4_innov_var,"(b)",expression(italic(N==64)),vlines=c(1/8,0.4))

### Figure 178, bottom row of plots

fig_178(256,ar4_coeffs,ar4_innov_var,"(c)",expression(italic(N==256)))
fig_178(1024,ar4_coeffs,ar4_innov_var,"(d)",expression(italic(N==1024)))

### Figure 180 ###

fig_180 <- function(the_kernel,mult_p=FALSE,v_line=1/8,trans=function(x) x,big_y_ats=seq(-40,20,20),little_y_ats=seq(-50,30,10),tag="(a)",word="and",the_sdf=two_sided_ar4_sdf)
{
    N_freqs <- length(the_kernel)
    freqs <- seq(-0.5+1/N_freqs,0.5,length=N_freqs)
    ys <- trans(if(mult_p) the_kernel*the_sdf else the_kernel)
    plot(freqs,ys,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(little_y_ats[1],little_y_ats[length(little_y_ats)]),yaxs="i",ylab=paste("kernel",word,"AR(4) SDF"),
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 180",tag,sep=""))
    if(!mult_p) lines(freqs,trans(the_sdf))
    abline(v=v_line,lty="dotted")
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=big_y_ats,las=2)
    axis(2,at=little_y_ats,label=FALSE,tcl=-0.25)
    text(x=-0.35,y=0.88*diff(range(little_y_ats))+little_y_ats[1],tag,pos=2)
    box(bty="l")
}

temp <- ar_coeffs_to_sdf(ar4_coeffs,ar4_innov_var,N_pad=2048)$sdf
two_sided_ar4_sdf <- c(rev(temp[c(-1,-length(temp))]),temp)

temp <- fejer_kernel(64)$kernel
fejer_shift_1 <- circular_shift(temp,256)
fejer_shift_2 <- circular_shift(temp,820)

### Figure 180, top row of plots

fig_180(fejer_shift_1,trans=dB)
fig_180(fejer_shift_1,big=c(0,40,80),little=c(0,40,80),tag="(b)")

### Figure 180, 2nd row of plots

fig_180(fejer_shift_1,trans=dB,mult_p=TRUE,tag="(c)",word="times")
fig_180(fejer_shift_1,big=c(0,250,500),little=c(0,250,500),mult_p=TRUE,tag="(d)",word="times")

### Figure 180, 3rd row of plots

fig_180(fejer_shift_2,trans=dB,v_line=0.4,tag="(e)")
fig_180(fejer_shift_2,big=c(0,40,80),little=c(0,40,80),v_line=0.4,tag="(f)")

### Figure 180, bottom row of plots

fig_180(fejer_shift_2,trans=dB,mult_p=TRUE,v_line=0.4,tag="(g)",word="times")
fig_180(fejer_shift_2,big=c(0,1,2),little=c(0,1,2),mult_p=TRUE,v_line=0.4,tag="(h)",word="times")

### Figure 182 ###
###
### NOTE: fig_182 is virtually the same as fig_173, the only
###       difference being the addition of pad=2 in the call
##        to pgram (fig_173 uses the default pad=1)

fig_182 <- function(ts,coeffs,innov_var,y_ats,tag)
{
    the_pgram <- pgram(ts,center=FALSE,pad=2)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab=paste("AR(",length(coeffs),") spectra  (dB)",sep=""),
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main=paste("Figure 182",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    if(length(coeffs) == 4)
    {
        N <- length(ts)
        temp <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,N-1,innov_var,FALSE))
        lines(temp$freqs, dB(temp$sdf_ev), lwd=0.5)
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=10,tag,pos=2)
    box(bty="l")
}

### Figure 182, top row of plots

fig_182(ar2_1,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(a)")
fig_182(ar2_2,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(b)")

### Figure 182, 2nd row of plots

fig_182(ar2_3,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(c)")
fig_182(ar2_4,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(d)")

### Figure 182, 3rd row of plots

fig_182(ar4_1,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(e)")
fig_182(ar4_2,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(f)")

### Figure 182, bottom row of plots

fig_182(ar4_3,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(g)")
fig_182(ar4_4,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(h)")

### Figure 183a ###

fig_183a <- function(ts,tag,coeffs=ar4_coeffs)
{
    N <- length(ts)
    p <- length(coeffs)
    plot(0:3,ts[(N-3):N],
         xlim=c(1,6),xlab=expression(italic(t)),
         ylim=c(-5,5),ylab="AR(4) series",
         typ="o",axes=FALSE,
         main=paste("Figure 183a",tag,sep=""))
    pred <- as.vector(coeffs%*%ts[N:(N-p+1)])
    lines(3:4,c(ts[N],pred), type="b", pch=" ", lty="dotted")
    points(4,pred, pch=3)
    lines(4:7,ts[1:4], type="o")
    axis(1,at=1:6,labels=c(1021,NA,1023,0,1,2))
    axis(2,at=seq(-5,5,5),las=2)
    axis(2,at=seq(-5,5,1),label=FALSE,tcl=-0.25)
    text(x=5.5,y=4.5,tag,pos=2)
    box(bty="l")
}

### Figure 183a, top row

fig_183a(ar4_1,"(e)")
fig_183a(ar4_2,"(f)")

### Figure 183a, bottom row

fig_183a(ar4_3,"(g)")
fig_183a(ar4_4,"(h)")

### Figure 183b ###

fig_183b <- function(x,y)
{
    plot(x,y,
         xlim=c(-0.15,5.5),xaxs="i",xlab="absolute prediction error",
         ylim=c(-52,-25),yaxs="i",ylab="dB",
         typ="p",cex=0.625,axes=FALSE,
         main="Figure 183b")
    lines(lowess(x,y))
    abline(h=c(-47.20893,-30.3018),lty=c("dotted","dashed"))
    axis(1,at=0:5)
    axis(2,at=seq(-50,-30,10),las=2)
    box(bty="l")
}

set.seed(1)
N_rep <- 100
x_results <- rep(0,100)
y_results <- rep(0,100)
LD_ar4 <- step_down_LD_recursions(ar4_coeffs,ar4_innov_var,proc=FALSE)
for(n in 1:N_rep)
  {
    ar_ts <- sim_ar_process(1024,LD=LD_ar4)
    x_results[n] <- abs(as.numeric(ar_ts[1024:1021] %*% ar4_coeffs) - ar_ts[1])
    y_results[n] <- dB(mean(pgram(c(ar_ts,rep(0,1024)),center=FALSE)$sdfe[821:1025]))
  }

### Figure 183b

fig_183b(x_results,y_results)

### Figure 185 ###

fig_185 <- function(ys,big_y_ats=seq(-5,5,5),little_y_ats=NULL,y_lab="AR(4) series")
{
    N <- length(ys)
    plot(0:(N-1),ys,
         xlim=c(0,N),xlab=expression(italic(t)),
         ylim=c(big_y_ats[1],big_y_ats[length(big_y_ats)]),ylab=y_lab,
         typ="l",lwd=0.25,axes=FALSE,
         main="Figure 185")
    axis(1,at=seq(0,1024,512))
    axis(1,at=seq(0,1024,256),label=FALSE,tcl=-0.25)
    axis(2,at=big_y_ats,las=2)
    axis(2,at=little_y_ats,label=FALSE,tcl=-0.25)
    box(bty="l")
}

the_taper <- hanning_taper(1024)

### Figure 185, top to bottom plots

fig_185(ar4_1,little=seq(-5,5,1))
fig_185(the_taper,big=seq(0,0.06,0.02),y_lab="Hanning taper")
fig_185(the_taper*ar4_1,big=seq(-0.2,0.2,0.1),y_lab="tapered series")

### Figure 187 ###

fig_187 <- function(ts,coeffs,innov_var,y_ats,tag)
{
    the_dse <- direct_sdf_est(ts,hanning_taper(length(ts)),center=FALSE,pad=2)
    plot(the_dse$freqs,dB(the_dse$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab=paste("AR(",length(coeffs),") spectra  (dB)",sep=""),
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main=paste("Figure 187",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=10,tag,pos=2)
    box(bty="l")
}

### Figure 187, top row of plots

fig_187(ar2_1,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(a)")
fig_187(ar2_2,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(b)")

### Figure 187, 2nd row of plots

fig_187(ar2_3,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(c)")
fig_187(ar2_4,ar2_coeffs,ar2_innov_var,seq(0,25,5),"(d)")

### Figure 187, 3rd row of plots

fig_187(ar4_1,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(e)")
fig_187(ar4_2,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(f)")

### Figure 187, bottom row of plots

fig_187(ar4_3,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(g)")
fig_187(ar4_4,ar4_coeffs,ar4_innov_var,seq(0,150,50),"(h)")

### Figure 190 ###

fig_190 <- function(the_taper,left_tag,right_tag)
{
    N <- length(the_taper)
    plot(0:(N-1),the_taper,
         xlim=c(0,N),xlab=expression(italic(t)),
         ylim=c(0,0.3),ylab="data taper",
         typ="p",pch=20,cex=0.2,axes=FALSE,
         main=paste("Figure 190",left_tag,sep=""))
    axis(1,at=seq(0,64,32))
    axis(1,at=seq(0,64,16),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0.0,0.3,0.1),las=2)
    text(x=0,y=0.29,left_tag,pos=4)
    text(x=64,y=0.29,right_tag,pos=2)
    box(bty="l")
}

### Figure 190, left-hand column

fig_190(rectangular_taper(64),"(a)",expression(paste("rectangular (",italic(p==0),")",sep="")))
fig_190(cosine_taper(64,0.2),"(b)",expression(italic(p==0.2)))
fig_190(cosine_taper(64,0.5),"(c)",expression(italic(p==0.5)))
fig_190(hanning_taper(64),"(d)",expression(paste("Hanning (",italic(p==1),")",sep="")))

### Figure 190, right-hand column

fig_190(slepian_taper(64,1),"(e)",expression(italic(NW==1)))
fig_190(slepian_taper(64,2),"(f)",expression(italic(NW==2)))
fig_190(slepian_taper(64,4),"(g)",expression(italic(NW==4)))
fig_190(slepian_taper(64,8),"(h)",expression(italic(NW==8)))

### Figure 191 ###

fig_191 <- function(the_taper,left_tag,right_tag,v_line=NULL)
{
    temp <- spec_window(the_taper,pad_factor=16,fix_nulls_p=TRUE,first_p=FALSE)
    freqs <- temp$freqs
    ys <- dB(temp$sw)
    plot(freqs,ys,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-100,20),yaxs="i",ylab="spectral window  (dB)",
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 191",left_tag,sep=""))
    abline(v=v_line*c(-1,1),lty="dotted")
    ## add 3 dB down width
    i_max <- which.max(ys)
    three_dB_down <- ys[i_max] - 3
    i <- which(ys[i_max:length(ys)] <= three_dB_down)[1] + i_max - 1
    lines(freqs[c(2*i_max-i,i)],c(three_dB_down,three_dB_down))
    ## add variance width
    bw_v <- function(taper)
    {
        N <- length(taper)
        Nm1 <- N - 1
        autocor <- Re(fft(abs(fft(c(taper,rep(0,N)))^2)))/(2*N)
        return(sqrt(1 + sum(((-1)^(1:Nm1))*autocor[2:N]/(1:Nm1)^2)*12/pi^2))
    }
    lines(bw_v(the_taper)*c(-0.5,0.5),c(three_dB_down-5,three_dB_down-5))
    ## add autocorrelation width
    lines(B_H(the_taper)*c(-0.5,0.5),c(three_dB_down-10,three_dB_down-10))
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-100,20,20),las=2)
    axis(2,at=seq(-100,20,10),label=FALSE,tcl=-0.25)
    text(x=-0.5,y=10,left_tag,pos=4)
    text(x=0.5,y=10,right_tag,pos=2)
    box(bty="l")
}

### Figure 191, left-hand column

fig_191(rectangular_taper(64),"(a)","rectangular")
fig_191(cosine_taper(64,0.2),"(b)",expression(italic(p==0.2)))
fig_191(cosine_taper(64,0.5),"(c)",expression(italic(p==0.5)))
fig_191(hanning_taper(64),"(d)","Hanning")

### Figure 191, right-hand column

fig_191(slepian_taper(64,1),"(e)",expression(italic(NW==1)),v_line=1/64)
fig_191(slepian_taper(64,2),"(f)",expression(italic(NW==2)),v_line=1/32)
fig_191(slepian_taper(64,4),"(g)",expression(italic(NW==4)),v_line=1/16)
fig_191(slepian_taper(64,8),"(h)",expression(italic(NW==8)),v_line=1/8)

### Figure 193 ###

fig_193 <- function(the_taper,tag,coeffs=ar4_coeffs,innov_var=ar4_innov_var)
{
    ev_dse <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,length(the_taper)-1,innov_var,FALSE),the_taper,N_pad=1024)
    plot(ev_dse$freqs, dB(ev_dse$sdf_ev),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 193")
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=15,tag,pos=2)
    box(bty="l")
}

### Figure 193, left-hand column

fig_193(rectangular_taper(64),expression(paste("rectangular (",italic(p==0),")",sep="")))
fig_193(cosine_taper(64,0.2),expression(italic(p==0.2)))
fig_193(cosine_taper(64,0.5),expression(italic(p==0.5)))
fig_193(hanning_taper(64),expression(paste("Hanning (",italic(p==1),")",sep="")))

### Figure 193, right-hand column

fig_193(slepian_taper(64,1),expression(italic(NW==1)))
fig_193(slepian_taper(64,2),expression(italic(NW==2)))
fig_193(slepian_taper(64,4),expression(italic(NW==4)))
fig_193(slepian_taper(64,8),expression(italic(NW==8)))

### Figure 199 ###

fig_199 <- function(pw_filter,tag,right_p=FALSE,extra_p=FALSE,ts=ar4_2,coeffs=ar4_coeffs,innov_var=ar4_innov_var)
{
    pw_ts <- convolve(ts,pw_filter,type="filter")
    N_pad <- 2048
    pgram_pw_ts <- pgram(pw_ts,center=FALSE,pad=N_pad/length(pw_ts))
    freqs <- pgram_pw_ts$freqs
    squared_gain <- abs(fft(c(pw_filter,rep(0,N_pad-length(pw_filter))))[1:((N_pad/2)+1)])^2
    plot(freqs,dB(if(right_p) pgram_pw_ts$sdfe/squared_gain else pgram_pw_ts$sdfe), 
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="spectra  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 199",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=N_pad)$sdf
    lines(freqs,dB(if(right_p) the_ar_spec else the_ar_spec * squared_gain))
    if(extra_p)
    {
        N <- length(ts)
        L <- length(pw_filter)
        ar_acvs <- ar_coeffs_to_acvs(coeffs,N+2*L,innov_var,FALSE)
        pre_acvs <- rep(0,N-L+1)
        for(tau in 0:(N-L))
            for(k in 1:L)
                for(l in 1:L)
                {
                    pre_acvs[tau+1] <- pre_acvs[tau+1] + pw_filter[k]*pw_filter[l]*ar_acvs[abs(tau+k-l)+1]
                }
        temp <- ev_lag_window_sdf_estimator(pre_acvs,rep(1/sqrt(N-L+1),N-L+1),N_pad=N_pad)
        pc <- abs(fft(c(pw_filter,rep(0,N_pad-L)))[1:((N_pad/2)+1)])^2
        lines(0.25+temp$freqs[1:410],dB(temp$sdf_ev[1:410]/pc[1:410]),lwd=0.25)
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=15,tag,pos=2)
    box(bty="l")
}

LD_ar4 <- step_down_LD_recursions(ar4_coeffs,ar4_innov_var,FALSE)

### Figure 199, top row

fig_199(c(1,-ar4_coeffs),"(a)")
fig_199(c(1,-ar4_coeffs),"(b)",right_p=TRUE)

### Figure 199, 2nd row

fig_199(c(1,-0.99),"(c)")
fig_199(c(1,-0.99),"(d)",right_p=TRUE)

### Figure 199, 3rd row

fig_199(c(1,-LD_ar4$coeffs[[2]]),"(e)")
fig_199(c(1,-LD_ar4$coeffs[[2]]),"(f)",right_p=TRUE,extra_p=TRUE)

### Figure 199, bottom row

fig_199(c(1,-1.3,0.8),"(g)")
fig_199(c(1,-1.3,0.8),"(h)",right_p=TRUE,extra_p=TRUE)

### Figure 200 ###

fig_200 <- function(pwf_1,pwf_3,pwf_4)
{
    N_pad <- 2048
    squared_gain <- function(filter) abs(fft(c(filter,rep(0,N_pad-length(filter))))[1:((N_pad/2)+1)])^2
    freqs <- seq(0.0,0.5,1/N_pad)
    plot(freqs,dB(squared_gain(pwf_1)), 
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-50,30),yaxs="i",ylab="squared gain function  (dB)",
         typ="l",axes=FALSE,
         main="Figure 200")
    lines(freqs,dB(squared_gain(pwf_3)),lwd=0.25)
    lines(freqs,dB(squared_gain(pwf_4)),lty="dotted")
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,30,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

LD_ar4 <- step_down_LD_recursions(ar4_coeffs,ar4_innov_var,FALSE)

### Figure 200

fig_200(c(1,-ar4_coeffs),c(1,-LD_ar4$coeffs[[2]]),c(1,-1.3,0.8))

### Figure 206 ###

fig_206 <- function(the_taper,B_H_multiplier,tag,N_pad=8192)
{
    N_pad_half <- N_pad/2
    freqs <- (-(N_pad_half-1):N_pad_half)/N_pad
    N <- length(the_taper)
    temp <- abs(fft(c(the_taper,rep(0,N_pad-N))))
    H_abs <- c(temp[(N_pad_half+2):N_pad],temp[1:(N_pad_half+1)])
    B_H_taper <- B_H(the_taper)
    i <- round(N_pad*(1-B_H_taper*B_H_multiplier))
    H_abs_shifted <- c(H_abs[i:N_pad],H_abs[1:(i-1)])
    for_xlim <- 1/8 + 1/64
    plot(freqs,H_abs_shifted, 
         xlim=(1/8 + 1/64)*c(-1,1),xlab=expression(italic(v)),
         ylim=c(0,30),ylab=" ",
         typ="l",axes=FALSE,
         main=paste("Figure 206",tag,sep=""))
    lines(freqs,H_abs,lwd=0.5)
    lines(freqs,H_abs*H_abs_shifted,col="gray",lwd=2)
    abline(v=0,lty="dotted")
    abline(v=B_H_taper,lty="dotted")
    axis(1,at=seq(-1/8,1/8,1/8),labels=c("-1/8","0","1/8"))
    axis(1,at=seq(-1/2,1/2,1/64),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(0,30,10),las=2)
    text(1/8,28,tag,pos=2)
    box(bty="l")
}

### Figure 206, first row, left to right

fig_206(slepian_taper(64,2),0.5,"(a)")
fig_206(slepian_taper(64,2),1,"(b)")
fig_206(slepian_taper(64,2),2,"(c)")

### Figure 206, second row, left to right

fig_206(slepian_taper(64,4),0.5,"(d)")
fig_206(slepian_taper(64,4),1,"(e)")
fig_206(slepian_taper(64,4),2,"(f)")

### Figure 207 ###

fig_207 <- function(ts,tag_1,tag_2)
{
    N <- length(ts)
    the_pgram <- pgram(ts,center=FALSE,pad=2^(11-round(log2(N))))
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-40,20),yaxs="i",ylab="periodogram  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 207",tag_1,sep=""))
    abline(h=0)
    x_cc <- 7/16
    y_cc <- -30
    lines(c(x_cc,x_cc),y_cc+c(dB(2/qchisq(0.975,2)),dB(2/qchisq(0.025,2))),lwd=0.5)
    lines(x_cc+c(-0.5,0.5)*the_pgram$cc$width,c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,1/N),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=18,tag_1,pos=2)
    text(x=0.25,y=-35,tag_2,pos=1)
    box(bty="l")
}

set.seed(42)
ts_128 <- rnorm(128)

### Figure 207, first row, left to right

fig_207(ts_128[1:16],"(a)",expression(N==16))
fig_207(ts_128[1:32],"(b)",expression(N==32))

### Figure 207, second row, left to right

fig_207(ts_128[1:64],"(c)",expression(N==64))
fig_207(ts_128,"(d)",expression(N==128))

### Figure 208 ###

fig_208 <- function(taper,tag_1,tag_2,ts=ar2_1[1:128],coeffs=ar2_coeffs,innov_var=ar2_innov_var)
{
    N <- length(ts)
    the_dse <- direct_sdf_est(ts,taper,center=FALSE,pad=2^(11-round(log2(N))))
    plot(the_dse$freqs,dB(the_dse$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-40,20),yaxs="i",ylab="AR(2) spectra  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 208",tag_1,sep=""))
    ar_sdf <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=2048)
    lines(ar_sdf$freqs,dB(ar_sdf$sdf))
    x_cc <- 0.4
    y_cc <- -30
    lines(c(x_cc,x_cc),y_cc+c(dB(2/qchisq(0.975,2)),dB(2/qchisq(0.025,2))),lwd=0.5)
    lines(x_cc+c(-0.5,0.5)*the_dse$cc$width,c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=18,tag_1,pos=2)
    text(x=0.04,y=-36,tag_2,pos=4)
    box(bty="l")
}

### Figure 208, first row, left to right

fig_208(default_taper(128),"(a)","periodogram")
fig_208(slepian_taper(128,2),"(b)",expression(italic(NW==2)))

### Figure 208, second row, left to right

fig_208(slepian_taper(128,4),"(c)",expression(italic(NW==4)))
fig_208(slepian_taper(128,8),"(d)",expression(italic(NW==8)))

### Figure 210a ###

fig_210a <- function(right_p=FALSE)
{
    if(!right_p)
    {
        xs <- seq(0,10,0.01)
        ys <- exp(-xs/2)/2
        plot(xs,ys,
             xlim=c(0,10),xaxs="i",xlab=expression(italic(u)),
             ylim=c(0,0.5),yaxs="i",ylab="PDF",
             typ="l",lwd=0.5,axes=FALSE,
             main="Figure 210a(a)")
        xs_inner <- seq(5.9915,10,0.01)
        ys_inner <- exp(-xs_inner/2)/2
        polygon(c(5.9915,xs_inner,10),c(0,ys_inner,0),col="gray",border=NA)
        abline(v=2,lty="dotted")
        axis(1,at=seq(0,10,5))
        axis(1,at=seq(0,10,1),label=FALSE,tcl=-0.25)
        axis(2,at=seq(0,0.5,0.5),las=2)
        axis(2,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
        text(9.85,0.44,"(a)",pos=2)
    }
    else
    {
        xs <- seq(-20,20,0.04)
        pdf_log_chi2 <- function(x)
        {
            temp <- 10^(x/10)
            return(log(10) * temp * exp(-temp/2)/20)
        }
        ys <- pdf_log_chi2(xs)
        plot(xs,ys,
             xlim=c(-20,20),xaxs="i",xlab=expression(paste(italic(v),"  (dB)")),
             ylim=c(0,0.1),yaxs="i",ylab="PDF",
             typ="l",lwd=0.5,axes=FALSE,
             main="Figure 210a(b)")
        xs_inner <- seq(-20,-9.8891,0.04)
        ys_inner <- pdf_log_chi2(xs_inner)
        polygon(c(-20,xs_inner,-9.8891),c(0,ys_inner,0),col="gray",border=NA)
        abline(v=dB(2/exp(-digamma(1))),lty="dotted")
        axis(1,at=seq(-20,20,10))
        axis(2,at=seq(0,0.1,0.1),las=2)
        axis(2,at=seq(0,0.1,0.01),label=FALSE,tcl=-0.25)
        text(19.4,0.088,"(b)",pos=2)
    }
    box(bty="l")
}

### Figure 210a, left-hand plot

fig_210a()

### Figure 210a, right-hand plot

fig_210a(right_p=TRUE)

### Figure 210b ###

fig_210b <- function(ts,right_p=FALSE)
{
    trans <- if(right_p) dB else function(x) x
    N <- length(ts)
    the_pgram <- pgram(ts,center=FALSE)
    plot(the_pgram$freqs[-c(1,65)],trans(the_pgram$sdfe[-c(1,65)]),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=if(right_p) c(-40,20) else c(0,10),yaxs="i",ylab=paste("periodogram",if(right_p) "  (dB)" else NULL,sep=""),
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 210b",if(right_p) "(b)" else "(a)",sep=""))
    abline(h=trans(2))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=if(right_p) seq(-40,20,20) else seq(0,10,5),las=2)
    axis(2,at=if(right_p) seq(-40,20,10) else seq(0,10,1),label=FALSE,tcl=-0.25)
    text(x=0.5,y=if(right_p) 14 else 9,if(right_p) "(b)" else "(a)",pos=2)
    box(bty="l")
}

set.seed(4)
ts_128 <- rnorm(128)*sqrt(2)

### Figure 210b, left-hand plot

fig_210b(ts_128)

### Figure 210b, right-hand plot

fig_210b(ts_128,right_p=TRUE)

### Figure 212 ###

fig_212 <- function(N=64,N_pad=2048)
{
    taper_1 <- slepian_taper(N,1)
    taper_2 <- slepian_taper(N,2)
    taper_4 <- slepian_taper(N,4)
    taper_8 <- slepian_taper(N,8)
    R_1 <- abs(fft(c(taper_1^2,rep(0,N_pad-N)))[1:((N_pad/8)+1)])^2
    R_2 <- abs(fft(c(taper_2^2,rep(0,N_pad-N)))[1:((N_pad/8)+1)])^2
    R_4 <- abs(fft(c(taper_4^2,rep(0,N_pad-N)))[1:((N_pad/8)+1)])^2
    R_8 <- abs(fft(c(taper_8^2,rep(0,N_pad-N)))[1:((N_pad/8)+1)])^2
    freqs <- (0:(N_pad/8))/N_pad
    plot(freqs,R_1,
             xlim=c(0,0.13),xaxs="i",xlab=expression(paste(eta,"  (frequency lag)")),
             ylim=c(0,1),yaxs="i",ylab="correlation",
         typ="l",axes=FALSE,
         main="Figure 212")
    lines(freqs,R_2,lty="longdash")
    lines(freqs,R_4,lty="dashed")
    lines(freqs,R_8,lty="dotted")
    abline(v=c(1/64,1/32,1/16,1/8),
           lty=c("solid","longdash","dashed","dotted"))
    lines(c(sum(taper_1^4),0),c(0.6,0.6))
    lines(c(sum(taper_2^4),0),c(0.5,0.5),lty="longdash")
    lines(c(sum(taper_4^4),0),c(0.4,0.4),lty="dashed")
    lines(c(sum(taper_8^4),0),c(0.3,0.3),lty="dotted")
    axis(1,at=c(0,1/64,1/32,1/16,1/8),labels=c("0","1/64","1/32","1/16","1/8"))
    axis(1, at=seq(0,1/8,1/64), labels=FALSE, tcl=-0.25)
    axis(1, at=c(5/64), labels=c("5/64"), tcl=-0.25)
    axis(2, at=seq(0,1,0.5), las=2)
    axis(2, at=seq(0,1,0.1), labels=FALSE, tcl=-0.25, tcl=-0.25)
    box(bty="l")
}

### Figure 212

fig_212()

### Table 214 ###

N <- 64
the_tapers <- list(cosine_taper(N,0),
                   cosine_taper(N,0.2),
                   cosine_taper(N,0.5),
                   cosine_taper(N,1),
                   slepian_taper(N,1),
                   slepian_taper(N,2),
                   slepian_taper(N,4),
                   slepian_taper(N,8))
delta_f <- 1/N

### Table 214, first row (1.50 1.56 1.72 2.06 1.59 2.07 2.86 4.01)

round(unlist(lapply(the_tapers,B_H))/delta_f,2)

### Table 214, second row (1.00 1.11 1.35 1.93 1.40 1.99 2.81 3.97)

round(unlist(lapply(the_tapers,function(x) sum(x^4)))/delta_f,2)

### Table 214, third row (1.50 1.41 1.27 1.06 1.13 1.04 1.02 1.01)

round(unlist(lapply(the_tapers,B_H))/unlist(lapply(the_tapers,function(x) sum(x^4))),2)

### Figure 216 ###

fig_216 <- function(ts,big_y_ats=c(0,4,8),delta_y=1,inc,right_p=FALSE)
{
    y_upper_lim <- big_y_ats[length(big_y_ats)]
    temp <- pgram(ts,center=FALSE)
    N <- length(ts)
    zap_me <- c(1,if(is_even(N)) N/2+1 else NULL)
    freqs <- temp$freqs[-zap_me]
    the_pgram <- temp$sdfe[-zap_me]
    the_cumsum <- cumsum(the_pgram)
    xs <- if(right_p) freqs[-length(freqs)] else freqs 
    ys <- if(right_p) the_cumsum[-length(the_cumsum)]/the_cumsum[length(the_cumsum)] else the_pgram
    plot(xs,ys,
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=if(right_p) c(0,1) else c(0,y_upper_lim),yaxs="i",ylab=paste(if(right_p) "cumulative" else NULL,"periodogram"),
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 216",if(right_p) "(b)" else "(a)",sep=""))
    points(xs,ys)
    if(right_p)
    {
        M <- length(the_cumsum)
        D_0p05 <- 1.358/(sqrt(M-1) + 0.12 + 0.11/sqrt(M-1))
        L_u <- function(f)  D_0p05 - 1/(M-1) + N*f/(M-1)
        L_l <- function(f) -D_0p05 + N*f/(M-1)
        lines(c(0,0.5),c(L_u(0),L_u(0.5)))
        lines(c(0,0.5),c(L_l(0),L_l(0.5)))
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=if(right_p) c(0,1) else big_y_ats,las=2)
    if(!right_p) axis(2,at=seq(0,y_upper_lim,delta_y),label=FALSE,tcl=-0.25)
    text(x=0.5,y=if(right_p) 0.1 else 0.9*y_upper_lim,if(right_p) "(b)" else "(a)",pos=2)
    box(bty="l")
}

### Figure 216, left-hand plot

fig_216(ar2_1[1:32])

### Figure 216, right-hand plot

fig_216(ar2_1[1:32],right_p=TRUE)

### Figure 217 ###

fig_217 <- function(ts,coeffs,innov_var,tag)
{
    the_dct_pgram <- pgram(c(ts,rev(ts)),center=FALSE)
    the_dct_pgram$sdfe[1] <- the_dct_pgram$sdfe[1]/2
    plot(the_dct_pgram$freqs,dB(the_dct_pgram$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main=paste("Figure 217",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs, innov_var, N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=10,tag,pos=2)
    box(bty="l")
}

### Figure 217, top row of plots

fig_217(ar4_1,ar4_coeffs,ar4_innov_var,"(e)")
fig_217(ar4_2,ar4_coeffs,ar4_innov_var,"(f)")

### Figure 217, bottom row of plots

fig_217(ar4_3,ar4_coeffs,ar4_innov_var,"(g)")
fig_217(ar4_4,ar4_coeffs,ar4_innov_var,"(h)")

### Figure 218 ###

fig_218 <- function(N,coeffs,innov_var,tag_1,tag_2,N_pad_ar=1024)
{
    ar_acvs <- ar_coeffs_to_acvs(coeffs,N-1,ar4_innov_var,FALSE)
    ev_pgram <- ev_lag_window_sdf_estimator(ar_acvs) #,N_pad=N_pad)
    ev_DCT <- ev_DCTII(ar_acvs)
    plot(ev_pgram$freqs,dB(ev_pgram$sdf),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 218",tag_2,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=N_pad_ar)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf))
    lines(ev_DCT$freqs,dB(ev_DCT$sdf),lwd=2.0)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.25,y=-50,tag_1,pos=1)
    text(x=0.5,y=10,tag_2,pos=2)
    box(bty="l")
}

### Figure 218, top row of plots

fig_218(16,ar4_coeffs,ar4_innov_var,expression(italic(N==16)),"(a)")
fig_218(64,ar4_coeffs,ar4_innov_var,expression(italic(N==64)),"(b)")

### Figure 218, bottom row of plots

fig_218(256,ar4_coeffs,ar4_innov_var,expression(italic(N==256)),"(c)")
fig_218(1024,ar4_coeffs,ar4_innov_var,expression(italic(N==1024)),"(d)")

### Figures 223 and 224b ###

fig_223 <- function(ys,y_ats,x_lab,main="Figure 223")
{
    N <- length(ys)
    plot(0:(N-1),Re(ys),
         xlim=c(0,N),xlab=x_lab,
         ylim=c(y_ats[1],tail(y_ats,1)),ylab=" ",
         typ="n",axes=FALSE,
         main=main)
    if(!(sum(Im(ys)) == 0))
    {
        lines(0:(N-1), Im(ys), lwd=0.5, col="gray40")
        points(0:(N-1), Im(ys), pch=16, cex=0.5, col="gray40")
    }
    lines(0:(N-1), Re(ys), col="black")
    points(0:(N-1), Re(ys), pch=16, cex=0.5)
    axis(1,at=c(0,N/2,N))
    axis(2,at=y_ats,las=2)
    box(bty="l")
}

### Figure 223, top row, left-hand plot

(N <- length(earth_20))        # 20
(M <- next_power_of_2(2*N-1))  # 64
M-N  # 44
tXt <- c(earth_20-mean(earth_20),rep(0,M-N))

fig_223(tXt,seq(-6,6,6),expression(italic(t)))

### Figure 223, top row, right-hand plot

tXt_dft <- dft(tXt)

fig_223(tXt_dft,seq(-32.5,32.5,32.5),expression(italic(k)))

### Figure 223, 2nd row, left-hand plot

tht <- c(hanning_taper(N),rep(0,M-N))

fig_223(tht,c(0,1),expression(italic(t)))

### Figure 223, 2nd row, right-hand plot

tht_dft <- dft(tht)

fig_223(tht_dft,seq(-4,4,4),expression(italic(k)))

### Figure 223, 3rd row, left-hand plot

thttXt <- tht*tXt

fig_223(thttXt,seq(-2,2,2),expression(italic(t)))

### Figure 223, 3rd row, right-hand plot

thttXt_dft <- dft(thttXt)

fig_223(thttXt_dft,seq(-6,6,6),expression(italic(k)))

### Figure 223, 4th row, left-hand plot

tSdk <- abs(thttXt_dft)^2
tsdtau <- Re(inverse_dft(tSdk))

fig_223(tsdtau,seq(-8,8,8),expression(tau))

### Figure 223, 4th row, right-hand plot

fig_223(tSdk,seq(0,60,30),expression(italic(k)))

### Figure 224b, top row, left-hand plot

fig_223(tXt,seq(-6,6,6),expression(italic(t)),main="Figure 224b")

### Figure 224b, top row, right-hand plot

fig_223(tXt_dft,seq(-32.5,32.5,32.5),expression(italic(k)),main="Figure 224b")

### Figure 224b, 2nd row, left-hand plot

tSpk <- abs(tXt_dft)^2/N
tsptau <- Re(inverse_dft(tSpk))

fig_223(tsptau,seq(-8,8,8),expression(tau),main="Figure 224b")

### Figure 224b, 2nd row, right-hand plot

fig_223(tSpk,seq(0,60,30),expression(italic(k)),main="Figure 224b")

### Figure 225 ###

fig_225 <- function(ts,delta_t=1/4)
{
    N <- length(ts)
    plot((0:(N-1))*delta_t,ts,
         xlim=c(0,N*delta_t),xlab="time  (sec)",
         ylim=c(-1200,1700),ylab="relative height",
         typ="l",axes=FALSE,
         main="Figure 225")
    axis(1,at=seq(0,256,64))
    axis(2,at=seq(-1000,1000,1000),las=2)
    axis(2,at=seq(-1000,1500,500),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 225

fig_225(ocean_wave)

### Figures 226 and 227 ###

fig_226 <- function(ts,taper,tag_1,tag_2,pad=1,h_line=0,v_line_p=FALSE,delta_t=1/4,main="Figure 226")
{
    dse <- direct_sdf_est(ts,taper,center=TRUE,delta_t=delta_t,pad=pad)
    plot(dse$freqs,dB(dse$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main=paste(main,tag_2,sep=""))
    cc <- dse$cc
    x_cc <- 0.16
    y_cc <- 15
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
    if(v_line_p) lines(c(0.16,0.16),c(68.5,80),lwd=0.5)
    abline(h=h_line,lty="dashed",lwd=0.5)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.0,80,tag_1,pos=1)
    text(1.9,80,tag_2,pos=1)
    box(bty="l")
}

### Figure 226, top plot

fig_226(ocean_wave,default_taper(1024),"periodogram","(a)",v=TRUE)

### Figure 226, 2nd plot

fig_226(ocean_wave,slepian_taper(1024,1),expression(paste("Slepian, ", italic(NW==1/Delta[t]))),"(b)")

### Figure 226, 3rd plot

fig_226(ocean_wave,slepian_taper(1024,2),expression(paste("Slepian, ", italic(NW==2/Delta[t]))),"(c)")

### Figure 226, bottom plot

fig_226(ocean_wave,slepian_taper(1024,4),expression(paste("Slepian, ", italic(NW==4/Delta[t]))),"(d)")

### Figure 227

fig_226(ocean_wave,default_taper(1024),"periodogram",NULL,pad=2,h_line=NULL,main="Figure 227")

### Figure 228 ###

fig_228 <- function(ts,delta_t=0.001)
{
    N <- length(ts)
    plot((0:(N-1))*delta_t,ts,
         ,xlab="time  (sec)",
         ylab="speed",
         typ="l",axes=FALSE,
         main="Figure 228")
    axis(1,at=seq(0,2.0,0.5))
    axis(2,at=seq(5,15,5),las=2)
    axis(2,at=seq(0,20,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 228

fig_228(chaotic_beam)

### Figure 229 ###

fig_229 <- function(ts,taper,tag_1,tag_2,delta_t=0.001)
{
    dse <- direct_sdf_est(ts,taper,center=TRUE,delta_t=delta_t)
    plot(dse$freqs,dB(dse$sdfe),
         xlim=c(0,500),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-120,0),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main=paste("Figure 229",tag_2,sep=""))
    cc <- dse$cc
    x_cc <- 425
    y_cc <- -30
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
    abline(h=-86,lty="dashed",lwd=0.5)
    axis(1,at=seq(0,500,100))
    axis(1,at=seq(0,500,10),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-120,0,20),las=2)
    axis(2,at=seq(-120,0,10),label=FALSE,tcl=-0.25)
    text(250,0,tag_1,pos=1)
    text(475,0,tag_2,pos=1)
    box(bty="l")
}

### Figure 229, top plot

fig_229(chaotic_beam,default_taper(2048),"periodogram","(a)")

### Figure 229, 2nd plot

fig_229(chaotic_beam,slepian_taper(2048,1),expression(paste("Slepian, ", italic(NW==1/Delta[t]))),"(b)")

### Figure 229, 3rd plot

fig_229(chaotic_beam,slepian_taper(2048,2),expression(paste("Slepian, ", italic(NW==2/Delta[t]))),"(c)")

### Figure 229, bottom plot

fig_229(chaotic_beam,hanning_taper(2048),"Hanning (100% cosine)","(d)")

### Figure 231 ###

fig_231 <- function(ts,right_p=FALSE)
{
    temp <- pgram(ts,center=TRUE)
    N <- length(ts)
    zap_me <- c(1,N/2+2)
    freqs <- temp$freqs[-zap_me]
    the_pgram <- temp$sdfe[-zap_me]
    the_cumsum <- cumsum(the_pgram)
    xs <- if(right_p) freqs[-length(freqs)] else freqs 
    ys <- if(right_p) the_cumsum[-length(the_cumsum)]/the_cumsum[length(the_cumsum)] else the_pgram
    plot(xs,ys,
         xlim=c(0,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=if(right_p) c(0,1) else c(0,36),yaxs="i",ylab=paste(if(right_p) "cumulative" else NULL,"periodogram"),
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 231",if(right_p) "(b)" else "(a)",sep=""))
    if(right_p)
    {
        M <- length(the_cumsum)
        D_0p05 <- 1.358/(sqrt(M-1) + 0.12 + 0.11/sqrt(M-1))
        D_0p1  <- 1.224/(sqrt(M-1) + 0.12 + 0.11/sqrt(M-1))
        L_u <- function(f,D)  D - 1/(M-1) + N*f/(M-1)
        L_l <- function(f,D) -D + N*f/(M-1)
        lines(c(0,0.5),c(L_u(0,D_0p05),L_u(0.5,D_0p05)))
        lines(c(0,0.5),c(L_l(0,D_0p05),L_l(0.5,D_0p05)))
        lines(c(0,0.5),c(L_u(0,D_0p1),L_u(0.5,D_0p1)),lty="dashed")
        lines(c(0,0.5),c(L_l(0,D_0p1),L_l(0.5,D_0p1)),lty="dashed")
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=if(right_p) c(0,1) else c(0,18,36),las=2)
    if(!right_p) axis(2,at=seq(0,36,3),label=FALSE,tcl=-0.25)
    text(x=0.5,y=if(right_p) 0.1 else 32.4,if(right_p) "(b)" else "(a)",pos=2)
    box(bty="l")
}

### Figure 231, left-hand plot

fig_231(ocean_noise)

### Figure 231, right-hand plot

fig_231(ocean_noise,right_p=TRUE)

### NOTE: code to recreate Figure 239 is not provided - to do so
###       would reveal the solution to Exercise [6.11]!
