### R CODE FOR REPRODUCING CONTENT OF FIGURES AND TABLES IN CHAPTER 7 ...

ar2_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar2_1.txt")
ar4_1 <- scan("http://faculty.washington.edu/dbp/sauts/Data/ar4_1.txt")
earth_20 <- scan("http://faculty.washington.edu/dbp/sauts/Data/earth_20.txt")
ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_wave.txt")
rough <- scan("http://faculty.washington.edu/dbp/sauts/Data/rough.txt")
smooth <- scan("http://faculty.washington.edu/dbp/sauts/Data/smooth.txt")
ac_time_differences <- scan("http://faculty.washington.edu/dbp/sauts/Data/maser_deglitched.txt")
ac_fractional_frequencies <- diff(ac_time_differences)*100/6
ship_altitude <- scan("http://faculty.washington.edu/dbp/sauts/Data/ship_altitude_2048.txt")
co2 <- scan("http://faculty.washington.edu/dbp/sauts/Data/co2_1958_2017.txt")

### functions used to compute content of figures in Chapter 7 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/acvs_to_lag_window_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/bartlett_priestley_design_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/beta_W_swb_exact.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/break_up_sw.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U_ds_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U_ds_pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_W_swb_exact.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_W_P_swb_exact.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/circular_shift.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/compute_cD_and_Im.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/cosine_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/daniell_design_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dirichlet_kernel.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/discrepancy_measures.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_lwe.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/direct_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/evaluate_gcv_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ev_lag_window_sdf_estimator.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/gcv_dsp.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/hanning_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_even.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_windows.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_window_derived_tapers.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_window_to_smoothing_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/m_to_dsdse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/next_power_of_2.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/parzen_smoothing_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/rectangular_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/spec_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ts_to_lag_window_sdf_est.R")

### NOTE: to install the sapa library, uncomment the following three statements
###       and execute them:
###
### install.packages("devtools")
### devtools::install_github("wconstan/ifultools")
### devtools::install_github("wconstan/sapa")

library(sapa)

###

fig_lag_window <- function(lw,main=" ")
{
    taus <- 0:(length(lw)-1)
    plot(taus,lw,
         xlim=c(0,64),xaxs="i",xlab=expression(tau),
         ylim=c(-0.25,1.02),yaxs="i",ylab="lag window",
         typ="p",cex=0.125,axes=FALSE,
         main=main)
    abline(h=0,lwd=0.5)
    axis(1,at=seq(0,64,32))
    axis(1,at=seq(0,64,8),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(-0.2,1,0.2),label=FALSE,tcl=-0.25)
    text(x=63,y=0.85,"(a)",pos=2)
    box(bty="l")
}

###

fig_smoothing_window <- function(lw,neg_sidelobes_p=FALSE,sidelobe_fix_p=TRUE,nyquist_fix_p=TRUE,fix_value=convert_from_dB(-100),main=" ")
{
    sw <- lag_window_to_smoothing_window(lw,8*2048)
    if(neg_sidelobes_p)
    {
        temp <- break_up_sw(sw)
        plot(temp$sw_pos$freqs,dB(temp$sw_pos$sw),
             xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
             ylim=c(-80,20),yaxs="i",ylab="smoothing window  (dB)",
             typ="l",lwd=0.25,axes=FALSE,
             main=main)
        polygon(c(temp$sw_neg$freqs,head(temp$sw_neg$freqs,1)),
                careful_dB(-c(temp$sw_neg$sw,tail(temp$sw_neg$sw,1))),col="grey")
    }
    else
    {
        if(sidelobe_fix_p) for(j in 2:(length(sw$freqs)-1)) if(sw$sw[j-1] > sw$sw[j] && sw$sw[j+1] > sw$sw[j]) sw$sw[j] <- fix_value
        if(nyquist_fix_p) sw$sw[length(sw$sw)] <- fix_value
        plot(sw$freqs,careful_dB(sw$sw),
             xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
             ylim=c(-80,20),yaxs="i",ylab="smoothing window  (dB)",
             typ="l",lwd=0.25,axes=FALSE,
             main=main)
    }
    lines(c(0,sw$beta_W/2),rep(dB(sw$sw[1])-3,2),lwd=1)
    lines(c(0,sw$B_W/2),rep(dB(sw$sw[1])-6,2),lwd=1.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    text(x=0.4925,y=12.8,"(b)",pos=2)
    box(bty="l")
}

###

fig_spectral_window <- function(lw,taper,tag,neg_sidelobes_p=FALSE,main=" ")
{
    sw <- ev_lag_window_sdf_estimator(rep(1,64),taper,lw,N_pad=1024)
    if(neg_sidelobes_p)
    {
        temp <- break_up_sw(list(freqs=sw$freqs,
                                 sw=sw$sdf_ev))
        plot(temp$sw_pos$freqs,dB(temp$sw_pos$sw),
             xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
             ylim=c(-80,20),yaxs="i",ylab="spectral window  (dB)",
             typ="l",lwd=0.25,axes=FALSE,
             main=paste(main,tag,sep=""))
        polygon(c(temp$sw_neg$freqs,head(temp$sw_neg$freqs,1)),
                careful_dB(-c(temp$sw_neg$sw,tail(temp$sw_neg$sw,1))),col="grey")
    }
    else
    {
        plot(sw$freqs,dB(sw$sdf_ev),
             xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
             ylim=c(-80,20),yaxs="i",ylab="spectral window  (dB)",
             typ="l",lwd=0.25,axes=FALSE,
             main=paste(main,tag,sep=""))
    }
    lines(c(0,sw$B_U/2),rep(dB(sw$sdf_ev[1])-3,2),lwd=1)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    text(x=0.4925,y=12.8,tag,pos=2)
    box(bty="l")
}

###

fig_design_window <- function(dw,main=" ")
{
    plot(dw$freqs,careful_dB(dw$dw),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-80,20),yaxs="i",ylab="design window  (dB)",
         typ="l",lwd=0.25,axes=FALSE,
         main=main)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

###

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

ar4_innov_var <- 0.002
ar4_coeffs    <- c(2.7607, -3.8106, 2.6535, -0.9238)

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 246 ###

fig_246 <- function(ts,coeffs,innov_var,pad_factor,tag)
{
    the_pgram <- pgram(ts,pad=pad_factor,center=FALSE)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-40,20),yaxs="i",ylab="AR(2) spectra  (dB)",
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main=paste("Figure 246",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf),lwd=0.5)
    N_freqs <- length(the_pgram$freqs)
    if(pad_factor==1)
    {
        M <- 15
        ds <- filter(c(the_pgram$sdfe,rev(the_pgram$sdfe[c(-1,-N_freqs)])),rep(1/(2*M+1),2*M+1),circular=TRUE)[1:N_freqs]
        y_filter <- c(rep(-36,20),rep(-28,2*M+1),rep(-36,20))
        x_filter <- 0.04 + (1:(40+2*M+1))/1024
    }
    else
    {
        irs <- dnorm(seq(-35,35,0.5),sd=5*sqrt(3))
        irs <- irs/sum(irs)
        ds <- filter(c(the_pgram$sdfe,rev(the_pgram$sdfe[c(-1,-N_freqs)])),irs,circular=TRUE)[1:N_freqs]
        x_filter <- 0.041 + seq(0,70,0.5)/1024
        factor <- 8/dnorm(0,sd=30/sqrt(12))
        y_filter <- -36 + factor * dnorm(seq(-35,35,0.5),sd=30/sqrt(12))
    }
    lines(the_pgram$freqs,dB(ds))
    lines(x_filter,y_filter)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=10,tag,pos=2)
    box(bty="l")
}

### Figure 246, left-hand plot

fig_246(ar2_1,ar2_coeffs,ar2_innov_var,1,"(a)")

### Figure 246, right-hand plot

fig_246(ar2_1,ar2_coeffs,ar2_innov_var,2,"(b)")

### Figure 248 ###

fig_248 <- function(ys,y_lim,y_lab,tag)
{
    N <- length(ys)
    xs <- (-(N-1)):(N-1)
    plot(xs,c(rev(ys[-1]),ys),
         xlim=c(-128,128),xlab=expression(tau),
         ylim=y_lim,ylab=y_lab,
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 248",tag,sep=""))
    axis(1,at=seq(-128,128,64))
    axis(1,at=seq(-160,160,16),label=FALSE,tcl=-0.25)
    axis(2,at=c(y_lim[1],0,y_lim[2]),las=2)
    axis(2,at=seq(y_lim[1],y_lim[2],1),label=FALSE,tcl=-0.25)
    text(x=128,y=(y_lim[2]-y_lim[1])*0.8+y_lim[1],tag,pos=2)
    box(bty="l")
}

ar2_1_acvs <- acvs(ar2_1,center=FALSE)$acvs
parzen_lw <- sapply(0:1023,parzen_lag_window,64)

### Figure 248, plots from top to bottom

fig_248(ar2_1_acvs,c(-2,2),"AR(2) ACVS","(a)")
fig_248(parzen_lw,c(0,1),"lag window","(b)")
fig_248(parzen_lw*ar2_1_acvs,c(-2,2),"windowed ACVS","(c)")

### Figure 249 ###

fig_249 <- function(ts,coeffs,innov_var)
{
    the_pgram <- pgram(ts,pad=2,center=FALSE)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-40,20),yaxs="i",ylab="AR(2) spectra  (dB)",
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main="Figure 249")
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(the_ar_spec$freqs,dB(the_ar_spec$sdf),lwd=0.5)
    the_lwe <- ts_to_lag_window_sdf_est(ts,m=64,lag_window=parzen_lag_window,center=FALSE)
    lines(the_lwe$freqs,dB(the_lwe$sdfe))
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 249

fig_249(ar2_1,ar2_coeffs,ar2_innov_var)
    
### Figure 252 ###

fig_252 <- function()
{
    freqs <- seq(-0.5,0.5,0.001)
    smoo_wind <- sapply(freqs,parzen_smoothing_window,64)
    plot(freqs,smoo_wind,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,50),ylab="smoothing window",
         typ="l",lwd=0.25,axes=FALSE,
         main="Figure 252")
    m <- 64
    beta_W <- beta_W_swb_exact(sapply(0:m,parzen_lag_window,m))
    B_W    <- B_W_swb_exact(sapply(0:m,parzen_lag_window,m))
    B_W_P  <- B_W_P_swb_exact(sapply(0:m,parzen_lag_window,m))
    lines(c(-beta_W/2,beta_W/2),c(24,24),lwd=1)
    lines(c(-B_W/2,B_W/2),c(22,22),lwd=1.5)
    lines(c(-B_W_P/2,B_W_P/2),c(20,20),lwd=2)
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,50,10),las=2)
    box(bty="l")
}

### Figure 252

fig_252()

### Figure 254 ###

fig_254 <- function(bottom_p=FALSE)
{
    taus <- (-1023):1023
    N_dft <- if(bottom_p) 2048 else 1024
    M <- if(bottom_p) 30 else 15
    plot(taus,sapply(taus/N_dft,dirichlet_kernel,N=2*M+1),
         xlim=c(-1024,1024),xlab=expression(tau),
         ylim=c(-0.25,1),ylab="lag windows",
         typ="l",lwd=2,col="gray",axes=FALSE,
         main=paste("Figure 254",if(bottom_p) "(b)" else "(a)",sep=""))
    temp     <- dnorm(seq(-35,35,if(bottom_p) 0.5 else 1),sd=5*sqrt(3))
    irs_test <- temp/sum(temp)
    for_dft <- rep(0,N_dft)
    j_max <- if(bottom_p) 70 else 35
    for_dft[1:(j_max+1)] <- irs_test[(j_max+1):(2*j_max+1)]
    for_dft[(N_dft-j_max+1):N_dft] <- irs_test[1:j_max]
    lw_1_sided <- Re(fft(for_dft))[1:1024]
    lines(taus,c(rev(lw_1_sided[-1]),lw_1_sided),lwd=0.5)
    axis(1,at=seq(-1024,1024,512))
    axis(1,at=seq(-1024,1024,128),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    text(x=800,y=0.9,if(bottom_p) "(b)" else "(a)",pos=2)
    box(bty="l")
}

### Figure 254, top plot

fig_254()

### Figure 254, bottom plot

fig_254(bottom_p=TRUE)

### Table 255 ###

beta_g_swb <- function(g,N)
{
    g_less_0 <- g[-1]
    sqrt(24 * sum(((1:(length(g_less_0)))^2)*g_less_0))/N
}

B_g_swb <- function(g,N)  1/(sum(c(g,g[-1])^2) * N)

taus <- 0:1023

### Table 255, first row

lw_1 <- sapply(taus/1024,dirichlet_kernel,N=31)
g_1 <- rep(1/31,16)

round(beta_W_swb_exact(lw_1),5)  # 0.03025
round(B_W_swb_exact(lw_1),5)     # 0.01537
round(B_W_P_swb_exact(lw_1),5)   # 0.01537
round(beta_g_swb(g_1,1024),5)    # 0.03026
round(B_g_swb(g_1,1024),5)       # 0.03027
round(1/(g_1[1]*1024),5)         # 0.03027
round(31/1024,5)                 # 0.03027

### Table 255, second row

temp     <- dnorm(seq(-35,35,1),sd=5*sqrt(3))
irs_test <- temp/sum(temp)
irs_test[35:37]
g_2 <- irs_test[36:71]
for_dft <- rep(0,1024)
for_dft[1:36] <- irs_test[36:71]
for_dft[990:1024] <- irs_test[1:35]
lw_2 <- Re(fft(for_dft))

round(beta_W_swb_exact(lw_2),5)  # 0.02928
round(B_W_swb_exact(lw_2),5)     # 0.01522
round(B_W_P_swb_exact(lw_2),5)   # 0.01071
round(beta_g_swb(g_2,1024),5)    # 0.02929
round(B_g_swb(g_2,1024),5)       # 0.02998
round(1/(g_2[1]*1024),5)         # 0.0212 

### Table 255, third row

lw_3 <- sapply(taus/2048,dirichlet_kernel,N=61)
g_3 <- rep(1/61,31)

round(beta_W_swb_exact(lw_3),5)  # 0.02978
round(B_W_swb_exact(lw_3),5)     # 0.02979
round(B_W_P_swb_exact(lw_3),5)   # 0.0298 
round(beta_g_swb(g_3,2048),5)    # 0.02978
round(B_g_swb(g_3,2048),5)       # 0.02979
round(1/(g_3[1]*2048),5)         # 0.02979
round(61/2048,5)                 # 0.02979

### Table 255, fourth row

temp     <- dnorm(seq(-35,35,0.5),sd=5*sqrt(3))
irs_test <- temp/sum(temp)
irs_test[70:72]
g_4 <- irs_test[71:141]
for_dft <- rep(0,2048)
for_dft[1:71] <- irs_test[71:141]
for_dft[1979:2048] <- irs_test[1:70]
lw_4 <- Re(fft(for_dft))[1:1024]

round(beta_W_swb_exact(lw_4),5)  # 0.02928
round(B_W_swb_exact(lw_4),5)     # 0.02998
round(B_W_P_swb_exact(lw_4),5)   # 0.0212 
round(beta_g_swb(g_4,2048),5)    # 0.02928
round(B_g_swb(g_4,2048),5)       # 0.02998
round(1/(g_4[1]*2048),5)         # 0.0212  

### Figure 259 ###

fig_259_top <- function(ts,k=38,J=5)
{
    pgram_ts <- pgram(ts,center=FALSE)
    padded_pgram_ts <- pgram(ts,center=FALSE,pad=16)
    x <- padded_pgram_ts$freqs
    y <- padded_pgram_ts$sdfe
    plot(x,y,
         xlim=c(0.1,0.2),xlab=" ",
         ylim=c(0,max(y[ (x >= 0.1) & (x <= 0.2) ])),ylab="SDF estimates",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 259(a)")
    k_p_1 <- k+1
    good_guys <- (k_p_1-J):(k_p_1+J)
    points(pgram_ts$freqs[good_guys],pgram_ts$sdfe[good_guys])
    abline(v=k/256,lty="dotted")
    the_lwe <- ts_to_lag_window_sdf_est(ts,m=64,lag_window=parzen_lag_window,center=FALSE)
    lines(the_lwe$freqs,the_lwe$sdfe,lwd=1.5)
    sw_grid_freqs <- ((-J):J)/256
    sw_grid <- sapply(sw_grid_freqs,parzen_smoothing_window,m=64)
    points(k/256,sum(pgram_ts$sdfe[good_guys]*sw_grid)/256,pch=4)
    axis(1,at=pgram_ts$freqs[k_p_1]+c(-J,0,J)/256,label=c(expression(italic(f[k] + f[-J])),expression(italic(f[k])),expression(italic(f[k] + f[J]))))
    axis(1,at=pgram_ts$freqs,label=FALSE,tcl=-0.25)
    axis(2,las=2)
    text(0.2,11,"(a)")
    text(0.177,10.25,expression(paste(hat(S)^(D),"(.)")))
    text(0.206,2.2,expression(paste(hat(S)[m]^(LW),"(.)")),pos=2)
    box(bty="l")
}

fig_259_bottom <- function(k=38,J=5)
{
    x <- seq(-0.6,0.6,length=1025)
    y <- sapply(x,parzen_smoothing_window,m=64)
    plot(x,y,
         xlim=c(-0.05,0.05),xlab=" ",
         ylim=c(0,50),ylab="smoothing window",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 259(b)")
    sw_grid_freqs <- ((-J):J)/256
    sw_grid <- sapply(sw_grid_freqs,parzen_smoothing_window,m=64)
    points(sw_grid_freqs,sw_grid)
    abline(v=0,lty="dotted")
    lw <- sapply(0:255,parzen_lag_window,m=64)
    lw_2_sided <- c(rev(lw[-1]),lw)
    B_W <- 1/(sum(lw_2_sided))
    lines(B_W*c(-0.5,0.5), rep(22,2),lwd=1.5)
    axis(1,at=c(-J/256,0,J/256),label=c(expression(italic(f[-J])),expression(0),expression(italic(f[J]))))
    axis(1,at=seq(-4*J/256,4*J/256,1/256),label=FALSE,tcl=-0.25)
    axis(2,las=2)
    text(0.05,45,"(b)")
    text(0.012,40,expression(paste(italic(W[m]),"(.)")))
    text(0,24,expression(italic(B[W])))
    box(bty="l")
}

### Figure 259, top plot

fig_259_top(ar2_1[1:256])

### Figure 259, bottom plot

fig_259_bottom()

### Table 260 ###

tab_260 <- function(N,display_N=FALSE)
{
    if(display_N) cat(switch(5-floor(log10(N)),""," ","  ","   ","    ","     "),N,"  ")
    cat(round(N*sum(cosine_taper(N,0)^4),2),
        round(N*sum(cosine_taper(N,0.2)^4),2),
        round(N*sum(cosine_taper(N,0.5)^4),2),
        round(N*sum(cosine_taper(N,1)^4),2),
        "   ",
        round(N*sum(t(as.matrix(taper("dpss",N,n.taper=2,bandwidth=1)))[,1]^4),2),
        round(N*sum(t(as.matrix(taper("dpss",N,n.taper=2,bandwidth=2)))[,1]^4),2),
        round(N*sum(t(as.matrix(taper("dpss",N,n.taper=2,bandwidth=4)))[,1]^4),2),
        round(N*sum(t(as.matrix(taper("dpss",N,n.taper=2,bandwidth=8)))[,1]^4),2),
        "\n")
}


tab_260_via_Equation_196b <- function(N,display_N=FALSE)
{
    if(display_N) cat(switch(5-floor(log10(N)),""," ","  ","   ","    ","     "),N,"  ")
    cat(round(N*sum(cosine_taper(N,0)^4),2),
        round(N*sum(cosine_taper(N,0.2)^4),2),
        round(N*sum(cosine_taper(N,0.5)^4),2),
        round(N*sum(cosine_taper(N,1)^4),2),
        "   ",
        round(N*sum(slepian_taper(N,1)^4),2),
        round(N*sum(slepian_taper(N,2)^4),2),
        round(N*sum(slepian_taper(N,4)^4),2),
        round(N*sum(slepian_taper(N,8)^4),2),
            "\n")
}

### Table 260

tab_260(512)  # 1 1.12 1.35 1.94     1.34 1.96 2.8 3.98 

### following shows convergence of C_h approximation as N increases 

for(N in 2^(5:14)) tab_260(N,TRUE)

### same as above, but using Slepian approximation of Equation (196b)

for(N in 2^(5:14)) tab_260_via_Equation_196b(N,TRUE)

### Figure 266 ###

fig_266 <- function()
{
    x <- seq(1,1000,1)
    y_upper <- dB(x/qchisq(0.025,x))
    plot(log10(x),y_upper,
         xlim=c(log10(0.9005),3.0),xaxs="i",xlab=expression(nu),
         ylim=c(-10,30),ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 266")
    y_lower <- 10*log10(x/qchisq(0.975,x))
    lines(log10(x),y_lower)
    abline(h=0,lwd=0.5)
    abline(v=log10(c(2,66,146,581)),lwd=0.5)
    axis(1,at=seq(0,3,1), labels=c("1","10","100","1000"))
    axis(1,at=log10(c(2:9,seq(20,90,10),seq(200,900,100))), labels=FALSE, tcl=-0.25)
    axis(2,at=seq(-10,30,10),las=2)
    box(bty="l")
}

### Figure 266

fig_266()

### Figure 267 ###

fig_267 <- function(ts,coeffs,innov_var,tag,right_p=FALSE)
{
    slw64 <- ts_to_lag_window_sdf_est(ts,m=64,lag_window=parzen_lag_window, center=FALSE)
    transform <- if(right_p) dB else function(x) x
    plot(slw64$freqs,transform(slw64$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(if(right_p) -11 else 0,10),yaxs="i",ylab=paste("AR(2) spectra",if(right_p) "  (dB)"),
         typ="l",lwd=1.5,axes=FALSE,
         main=paste("Figure 267",tag,sep=""))
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(the_ar_spec$freqs,transform(the_ar_spec$sdf),lwd=0.5)
    edof <- slw64$cc$edof
    CI_lower <- edof*slw64$sdfe/qchisq(0.975,edof)
    CI_upper <- edof*slw64$sdfe/qchisq(0.025,edof)
    lines(slw64$freqs,transform(CI_lower),lty="dashed")
    lines(slw64$freqs,transform(CI_upper),lty="dashed")
    if(right_p)
    {
        ## add criss-cross
        x_cc <-  0.1
        y_cc <- -6.0
        cc_width <- slw64$cc$width
        cc_up    <- slw64$cc$up
        cc_down <-  slw64$cc$down
        lines(c(x_cc-cc_width/2,x_cc+cc_width/2), c(y_cc,y_cc), lwd=0.5)
        lines(c(x_cc,x_cc), c(y_cc-cc_down,y_cc+cc_up), lwd=0.5)
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=if(right_p) seq(-10,10,10) else seq(0,10,5),las=2)
    axis(2,at=if(right_p) seq(-15,10,1) else seq(0,10,1),label=FALSE,tcl=-0.25)
    text(x=0.5,y=9,tag,pos=2)
    box(bty="l")
}

### Figure 267, left-hand plot

fig_267(ar2_1,ar2_coeffs,ar2_innov_var,"(a)")

### Figure 267, right-hand plot

fig_267(ar2_1,ar2_coeffs,ar2_innov_var,"(b)",right_p=TRUE)

### Figure 270, Bartlett lag window ###

bartlett_lw <- sapply(0:63,bartlett_lag_window,30)

### Figure 270(a)
fig_lag_window(bartlett_lw,main="Figure 270(a)")
### Figure 270(b)
fig_smoothing_window(bartlett_lw,main="Figure 270(b)")
### Figure 270(c)
fig_spectral_window(bartlett_lw,rectangular_taper(64),"(c)",main="Figure 270")
### Figure 270(d)
fig_spectral_window(bartlett_lw,slepian_taper(64),"(d)",main="Figure 270")

### Figure 272a, Daniell design window ###

fig_design_window(daniell_design_window(20),main="Figure 272a")

### Figure 272b, Daniell lag window ###

daniell_lw <- sapply(0:63,daniell_lag_window,20)

### Figure 272b(a)
fig_lag_window(daniell_lw,main="Figure 272(a)")
### Figure 272b(b)
fig_smoothing_window(daniell_lw,TRUE,main="Figure 272(b)")
### Figure 272b(c)
fig_spectral_window(daniell_lw,rectangular_taper(64),"(c)",main="Figure 272")
### Figure 272b(d)
fig_spectral_window(daniell_lw,slepian_taper(64),"(d)",main="Figure 272")

### Figure 274a, Bartlett-Priestley design window ###

fig_design_window(bartlett_priestley_design_window(16.67),main="Figure 274a")
    
### Figure 274b(a), Bartlett-Priestley lag window ###

bp_lw <- sapply(0:63,bartlett_priestley_lag_window,16.67)

### Figure 274b(a)
fig_lag_window(bp_lw,main="Figure 274(a)")
### Figure 274b(b)
fig_smoothing_window(bp_lw,TRUE,main="Figure 274(b)")
### Figure 274b(c)
fig_spectral_window(bp_lw,rectangular_taper(64),"(c)",main="Figure 274")
### Figure 274b(d)
fig_spectral_window(bp_lw,slepian_taper(64),"(d)",main="Figure 274")

### Figure 276, Parzen lag window ###

parzen_lw <- sapply(0:63,parzen_lag_window,37)

### Figure 276(a)
fig_lag_window(parzen_lw,main="Figure 276(a)")
### Figure 276(b)
fig_smoothing_window(parzen_lw,side=FALSE,nyq=FALSE,main="Figure 276(b)")
### Figure 276(c)
fig_spectral_window(parzen_lw,rectangular_taper(64),"(c)",main="Figure 276")
### Figure 276(d)
fig_spectral_window(parzen_lw,slepian_taper(64),"(d)",main="Figure 276")

### Figure 278, Gaussian lag window ###

gaussian_lw <- sapply(0:63,gaussian_lag_window,16)

### Figure 278b(a)
fig_lag_window(gaussian_lw,main="Figure 278(a)")
### Figure 278b(b)
fig_smoothing_window(gaussian_lw,TRUE,main="Figure 278(b)")
### Figure 278b(c)
fig_spectral_window(gaussian_lw,rectangular_taper(64),"(c)",main="Figure 278")
### Figure 278b(d)
fig_spectral_window(gaussian_lw,slepian_taper(64),"(d)",main="Figure 278")

### Figure 279, Papoulis lag window ###

papoulis_lw <- sapply(0:63,papoulis_lag_window,34)

### Figure 279(a)
fig_lag_window(papoulis_lw,main="Figure 279(a)")
### Figure 279(b)
fig_smoothing_window(papoulis_lw,nyq=FALSE,main="Figure 279(b)")
### Figure 279(c)
fig_spectral_window(papoulis_lw,rectangular_taper(64),"(c)",main="Figure 279")
### Figure 279(d)
fig_spectral_window(papoulis_lw,slepian_taper(64),"(d)",main="Figure 279")

### Table 279 ###

### Table 279, asymptotic variance multipliers via Equation (260)

eqn_260 <- function(lw,Nm2,m) sum(sapply((-Nm1):Nm1,lw,m)^2)
Nm1 <- 999
m   <- 100

round(eqn_260(bartlett_lag_window,m,m)/m,2)              # 0.67
round(eqn_260(daniell_lag_window,Nm1,m)/m,1)             # 1   
round(eqn_260(bartlett_priestley_lag_window,Nm1,m)/m,2)  # 1.2 
round(eqn_260(parzen_lag_window,m,m)/m,2)                # 0.54
round(eqn_260(gaussian_lag_window,Nm1,m)/m,2)            # 1.25
round(eqn_260(papoulis_lag_window,m,m)/m,2)              # 0.59

### Table 279, nu multipliers via Equation (264b)

eqn_264b <- function(lw,Nm2,m) 2/sum(sapply((-Nm1):Nm1,lw,m)^2)
Nm1 <- 999
m   <- 100

round(m*eqn_264b(bartlett_lag_window,m,m),2)              # 3   
round(m*eqn_264b(daniell_lag_window,Nm1,m),1)             # 2   
round(m*eqn_264b(bartlett_priestley_lag_window,Nm1,m),2)  # 1.67
round(m*eqn_264b(parzen_lag_window,m,m),2)                # 3.71
round(m*eqn_264b(gaussian_lag_window,Nm1,m),2)            # 1.6 
round(m*eqn_264b(papoulis_lag_window,m,m),2)              # 3.41

### Table 279, B_W multipliers via Equation (251e)

eqn_251e <- function(lw,Nm1,m) 1/sum(sapply((-Nm1):Nm1,lw,m)^2)
Nm1 <- 999
m   <- 100

round(m*eqn_251e(bartlett_lag_window,m,m),2)              # 1.5
round(m*eqn_251e(daniell_lag_window,Nm1,m),1)             # 1
round(m*eqn_251e(bartlett_priestley_lag_window,Nm1,m),2)  # 0.83
round(m*eqn_251e(parzen_lag_window,m,m),2)                # 1.85
round(m*eqn_251e(gaussian_lag_window,Nm1,m),2)            # 0.8
round(m*eqn_251e(papoulis_lag_window,m,m),2)              # 1.7

### Table 279, beta_W multipliers via Equation (251b)

eqn_251b <- function(lw,Nm1,m) sqrt(1+12*sum(sapply(1:Nm1,lw,m)*(-1)^(1:Nm1)/(1:Nm1)^2)/pi^2)
Nm1 <- 999
m   <- 100

round(sqrt(m)*eqn_251b(bartlett_lag_window,m,m),2)        # 0.92
round(m*eqn_251b(daniell_lag_window,Nm1,m),2)             # 1
round(m*eqn_251b(bartlett_priestley_lag_window,Nm1,m),2)  # 0.77
round(m*eqn_251b(parzen_lag_window,m,m),2)                # 1.91
round(m*eqn_251b(gaussian_lag_window,Nm1,m),2)            # 0.78
round(m*eqn_251b(papoulis_lag_window,m,m),2)              # 1.73

### Figure 280, modified Daniell lag window ###

md_lw <- modified_daniell_lag_window_one_sided(3,64)

### Figure 280(a)
fig_lag_window(md_lw,main="Figure 280(a)")
### Figure 280(b)
fig_smoothing_window(md_lw,TRUE,main="Figure 280(b)")
### Figure 280(c)
fig_spectral_window(md_lw,rectangular_taper(64),"(c)",main="Figure 280")
### Figure 280(d)
fig_spectral_window(md_lw,slepian_taper(64),"(d)",main="Figure 280")

### Figure 282, reshaped lag window ###

reshaped_lw <- sapply(0:63,reshaped_lag_window)

### Figure 282(a)
fig_lag_window(reshaped_lw,main="Figure 282(a)")
### Figure 282(b)
fig_smoothing_window(reshaped_lw,TRUE,main="Figure 282(b)")

### Figure 282(c) ###

fig_282c <- function(ts,pad_factor=32)
{
    reshaped_lw <- sapply(0:63,reshaped_lag_window)
    lwe_reshaped <- ts_to_lag_window_sdf_est(ts,lag_window=reshaped_lw,center=FALSE,pad=pad_factor)
    temp <- break_up_sw(list(freqs=lwe_reshaped$freqs,
                             sw=lwe_reshaped$sdfe))
    plot(temp$sw_pos$freqs,dB(temp$sw_pos$sw),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-80,20),yaxs="i",ylab="smoothing window  (dB)",
         typ="l",axes=FALSE,
         main="Figure 282(c)")
    polygon(c(temp$sw_neg$freqs,head(temp$sw_neg$freqs,1)),
            careful_dB(-c(temp$sw_neg$sw,tail(temp$sw_neg$sw,1))),col="grey")
    lwe_parzen <- ts_to_lag_window_sdf_est(ts,m=37,taper_parameter=4,taper=slepian_taper,lag_window=parzen_lag_window,center=FALSE)
    lines(lwe_parzen$freqs,dB(lwe_parzen$sdfe),lwd=0.25)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    text(x=0.4925,y=12.8,"(c)",pos=2)
    box(bty="l")
}

### Figure 282(c)

fig_282c(ar4_1[1:64])

### Figure 282(d) ###

fig_282d <- function(coeffs,innov_var)
{
    the_ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    plot(the_ar_spec$freqs,dB(the_ar_spec$sdf),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-80,20),yaxs="i",ylab="smoothing window  (dB)",
         typ="l",axes=FALSE,
         main="Figure 282(d)")
    ev_lwe <- ev_lag_window_sdf_estimator(ar_coeffs_to_acvs(coeffs,63,innov_var,FALSE),lag_window=sapply(0:63,reshaped_lag_window))
    lines(ev_lwe$freqs,dB(ev_lwe$sdf_ev),lwd=0.25)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    text(x=0.4925,y=12.8,"(d)",pos=2)
    box(bty="l")
}

### Figure 282(d)

fig_282d(ar4_coeffs,ar4_innov_var)

### Figure 283 ###

fig_283 <- function(lw,taper)
{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    ## par(mar=c(5.1,4.1,4.1,2.1)) is default (bottom, left, top, right)
    par(mar=c(4.1,5.1,4.1,4.1))
    N <- length(taper)
    plot(1:N,lw,
         xlim=c(1,N),xlab=expression(paste(italic(t)," or ",tau)),
         ylim=c(0,1),ylab="data taper",
         typ="l",lwd=4,col="gray",axes=FALSE,
         main="Figure 283")
    lines(1:N,taper/taper[(N+1)/2],lwd=0.5)
    axis(1,at=seq(1,127,63),label=c("0 or -63","63 or 0","126 or 63"))
    axis(2,at=c(0,1),label=c("0","0.17"),las=2)
    axis(4,at=c(0,1),las=2)
    mtext("lag window",4,2)
    box(bty="u")
}

### Figure 283

fig_283(sapply((-63):63,parzen_lag_window,64),slepian_taper(127,NW=3.5))

### Figure 284, 100% cosine lag window ###

cos_lw <- sapply(0:63,split_cosine_lag_window,c(27,1))

### Figure 284(a)
fig_lag_window(cos_lw,main="Figure 284(a)")
### Figure 284(b)
fig_smoothing_window(cos_lw,TRUE,main="Figure 284(b)")
### Figure 284(c)
fig_spectral_window(cos_lw,rectangular_taper(64),"(c)",main="Figure 284")
### Figure 284(d)
fig_spectral_window(cos_lw,slepian_taper(64),"(d)",TRUE,main="Figure 284")

### Figure 285, autocorrelated 100% cosine lag window ###

auto_cos_lw <- Re(inverse_dft(abs(dft(c(cosine_taper(41,1),rep(0,215))))^2))[1:64]

### Figure 285(a)
fig_lag_window(auto_cos_lw,main="Figure 285(a)")
### Figure 285(b)
fig_smoothing_window(auto_cos_lw,main="Figure 285(b)")
### Figure 285(c)
fig_spectral_window(auto_cos_lw,rectangular_taper(64),"(c)",main="Figure 285")
### Figure 285(d)
fig_spectral_window(auto_cos_lw,slepian_taper(64),"(d)",main="Figure 285")

### Figure 286, Slepian lag window ###

temp <- slepian_taper(127,10.2)
slepian_lw <- temp[64:127]/temp[64]

### Figure 286(a)
fig_lag_window(slepian_lw,main="Figure 286(a)")
### Figure 286(b)
fig_smoothing_window(slepian_lw,main="Figure 286(b)")
### Figure 286(c)
fig_spectral_window(slepian_lw,rectangular_taper(64),"(c)",main="Figure 286")
### Figure 286(d)
fig_spectral_window(slepian_lw,slepian_taper(64),"(d)",main="Figure 286")

### Figure 287 ###

fig_287 <- function(the_taper,tag)
{
    N <- length(the_taper)
    plot(0:(N-1),the_taper,
         xlim=c(0,N),xlab=expression(italic(t)),
         ylim=c(-0.04,0.3),ylab="data taper",
         typ="p",pch=20,cex=0.2,axes=FALSE,
         main=paste("Figure 287",tag,sep=""))
    axis(1,at=seq(0,64,32))
    axis(1,at=seq(0,64,16),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,0.3,0.1),las=2)
    text(12,0.25,tag)
    box(bty="l")
}

xs <- 2*(1:64)/65 - 1

temp <- sapply(xs,h_daniell_lag_window,0.754)
taper_a <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_bartlett_priestley_lag_window,0.528)
taper_b <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_parzen_lag_window)
taper_c <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_gaussian_lag_window,0.435)
taper_d <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_papoulis_lag_window)
taper_e <- temp/sqrt(sum(temp^2))

taper_f <- slepian_taper(64,3.17)

### Figure 287, top row

fig_287(taper_a,"(a)")
fig_287(taper_b,"(b)")

### Figure 287, middle row

fig_287(taper_c,"(c)")
fig_287(taper_d,"(d)")

### Figure 287, bottom row

fig_287(taper_e,"(e)")
fig_287(taper_f,"(f)")

### Figure 288 ###

fig_288 <- function(the_taper,tag)
{
    temp <- spec_window(the_taper,pad_factor=16,fix_nulls_p=TRUE,first_p=FALSE)
    freqs <- temp$freqs
    ys <- dB(temp$sw)
    plot(freqs,ys,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-100,20),yaxs="i",ylab="spectral window  (dB)",
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 288",tag,sep=""))
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
    text(x=-0.5,y=10,tag,pos=4)
    box(bty="l")
}

xs <- 2*(1:64)/65 - 1

temp <- sapply(xs,h_daniell_lag_window,0.754)
taper_a <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_bartlett_priestley_lag_window,0.528)
taper_b <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_parzen_lag_window)
taper_c <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_gaussian_lag_window,0.435)
taper_d <- temp/sqrt(sum(temp^2))

temp <- sapply(xs,h_papoulis_lag_window)
taper_e <- temp/sqrt(sum(temp^2))

taper_f <- slepian_taper(64,3.17)

### Figure 288, top row

fig_288(taper_a,"(a)")
fig_288(taper_b,"(b)")

### Figure 288, middle row

fig_288(taper_c,"(c)")
fig_288(taper_d,"(d)")

### Figure 288, bottom row

fig_288(taper_e,"(e)")
fig_288(taper_f,"(f)")

### Figure 289 ###

fig_289 <- function(ts,lag_window,m,tag,the_taper=hanning_taper(length(ts)))
{
    dse <- direct_sdf_est(ts,the_taper,center=FALSE,pad=16)
    plot(dse$freqs,dB(dse$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="AR(4) spectra  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main=paste("Figure 289",tag,sep=""))
    lwe <- ts_to_lag_window_sdf_est(ts,m=m,taper=the_taper,lag_window=lag_window,center=FALSE,pad=2)
    lines(lwe$freqs,dB(lwe$sdfe),lwd=1.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(0.49,10,tag,pos=2)
    box(bty="l")
}

### Figure 289, top row

fig_289(ar4_1[1:64],bartlett_lag_window,30,"(a)")
fig_289(ar4_1[1:64],daniell_lag_window,20,"(b)")

### Figure 289, bottom row

fig_289(ar4_1[1:64],parzen_lag_window,37,"(c)")
fig_289(ar4_1[1:64],gaussian_lag_window,16,"(d)")

### Figure 291 ###

fig_291 <- function(ts,pad_factor=32,the_taper=hanning_taper(length(ts)))
{
    dse <- direct_sdf_est(ts,the_taper,center=FALSE,pad=16)
    plot(dse$freqs,dB(dse$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="AR(4) spectra  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 291")
    lwe <- ts_to_lag_window_sdf_est(ts,27,taper=the_taper,lag_window=split_cosine_lag_window,center=FALSE,pad=pad_factor)
    temp <- break_up_sw(list(freqs=lwe$freqs,
                             sw=lwe$sdfe))
    lines(temp$sw_pos$freqs,dB(temp$sw_pos$sw),lwd=1.5)
    polygon(c(temp$sw_neg$freqs,head(temp$sw_neg$freqs,1)),
            careful_dB(-c(temp$sw_neg$sw,tail(temp$sw_neg$sw,1))),lwd=1.5,col="grey")
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 291

fig_291(ar4_1[1:64])

### Figure 292 ###

fig_292 <- function(m,coeffs,innov_var,tag,N_pad=2048)
{
    N_pad_over_2 <- N_pad/2  # assumes N_pad is even
    ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=N_pad)
    plot(ar_spec$freqs,dB(ar_spec$sdf),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 292")
    ar_acvs <- ar_coeffs_to_acvs(coeffs, N_pad_over_2, innov_var, FALSE)
    ar_acvs_wrapped <- c(ar_acvs[-(N_pad_over_2+1)],rev(ar_acvs[-1]))
    temp <- sapply(0:1024,parzen_lag_window,m)
    parzen_wrapped <- c(temp[-(N_pad_over_2+1)],rev(temp[-1]))
    sw_and_ar_spec_convolved <- Re(dft(ar_acvs_wrapped*parzen_wrapped)[1:(N_pad_over_2+1)])
    lines(ar_spec$freqs,dB(sw_and_ar_spec_convolved),lwd=0.5)
    lw <- sapply(0:m,parzen_lag_window,m)
    sw <- lag_window_to_smoothing_window(lw,8*N_pad)
    freqs_two_sided <- c(rev(-sw$freqs[-1]),sw$freqs)
    sw_two_sided <- c(rev(sw$sw[-1]),sw$sw)
    lines(freqs_two_sided+0.1,dB(sw_two_sided)-31.8,lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=13.3,tag,pos=2)
    box(bty="l")
}

### Figure 292, top row

fig_292(32,ar4_coeffs,ar4_innov_var,expression(italic(m == 32)))
fig_292(64,ar4_coeffs,ar4_innov_var,expression(italic(m == 64)))

### Figure 292, bottom row

fig_292(128,ar4_coeffs,ar4_innov_var,expression(italic(m == 128)))
fig_292(256,ar4_coeffs,ar4_innov_var,expression(italic(m == 256)))

### Figure 294 ###

fig_294 <- function(ts,m,coeffs,innov_var,tag)
{
    the_pgram <- pgram(ts,center=FALSE)
    plot(the_pgram$freqs,dB(the_pgram$sdf),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-40,20),yaxs="i",ylab="AR(2) spectra  (dB)",
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main="Figure 294")
    ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(ar_spec$freqs,dB(ar_spec$sdf),lwd=0.5)
    lw <- ts_to_lag_window_sdf_est(ts,m,lag_window=parzen_lag_window,center=FALSE)
    lines(lw$freqs, dB(lw$sdfe), lwd=1)
    x_left <- 0.025
    x_mid <- x_left + lw$cc$width/2
    y_cc <- -30.0
    lines(c(x_left,x_left+lw$cc$width), c(y_cc,y_cc), lwd=0.5)
    lines(c(x_mid,x_mid), c(y_cc-lw$cc$down,y_cc+lw$cc$up), lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=15,tag,pos=2)
    box(bty="l")
}

### Figure 294, top row

fig_294(ar2_1,4,ar2_coeffs,ar2_innov_var,expression(italic(m == 4)))
fig_294(ar2_1,16,ar2_coeffs,ar2_innov_var,expression(italic(m == 16)))

### Figure 294, bottom row

fig_294(ar2_1,29,ar2_coeffs,ar2_innov_var,expression(italic(m == 29)))
fig_294(ar2_1,64,ar2_coeffs,ar2_innov_var,expression(italic(m == 64)))

### Figures 295 and 301 ###

fig_295 <- function(ts,m,coeffs,innov_var,tag,main="Figure 295")
{
    the_dse <- direct_sdf_est(ts,hanning_taper(length(ts)),center=FALSE)
    plot(the_dse$freqs,dB(the_dse$sdf),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="AR(4) spectra  (dB)",
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main=main)
    ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(ar_spec$freqs,dB(ar_spec$sdf),lwd=0.5)
    lw <- ts_to_lag_window_sdf_est(ts,m,lag_window=parzen_lag_window,taper=hanning_taper,center=FALSE)
    lines(lw$freqs, dB(lw$sdfe), lwd=1)
    x_left <- 0.025
    x_mid <- x_left + lw$cc$width/2
    y_cc <- -50.0
    lines(c(x_left,x_left+lw$cc$width), c(y_cc,y_cc), lwd=0.5)
    lines(c(x_mid,x_mid), c(y_cc-lw$cc$down,y_cc+lw$cc$up), lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=13.3,tag,pos=2)
    box(bty="l")
}

### Figure 295, top row

fig_295(ar4_1,64,ar4_coeffs,ar4_innov_var,expression(italic(m == 64)))
fig_295(ar4_1,128,ar4_coeffs,ar4_innov_var,expression(italic(m == 128)))

### Figure 295, bottom row

fig_295(ar4_1,179,ar4_coeffs,ar4_innov_var,expression(italic(m == 179)))
fig_295(ar4_1,565,ar4_coeffs,ar4_innov_var,expression(italic(m == 565)))

### Figure 301

fig_295(ar4_1,102,ar4_coeffs,ar4_innov_var,expression(italic(m == 102)),main="Figure 301")

### Figure 298 ###

fig_298 <- function()
{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    ## par(mar=c(5.1,4.1,4.1,2.1)) is default (bottom, left, top, right)
    par(mar=c(5.1,5.1,6.0,4.1))
    plot(c(0.21,0.21,0.77,0.77),c(0,0.7,0.7,0),
         xlim=c(0,1),xaxs="i",xlab=" ",
         ylim=c(0,1),yaxs="i",ylab="SDF",
         typ="l",lwd=1.5,axes=FALSE,
         main="Figure 298")
    axis(1,at=c(0.22,0.32,0.72), labels=expression(italic(f[l]),italic(f[l+1]),italic(f[l+K-1])))
    axis(1,at=seq(0.02,0.92,0.1), labels=FALSE, tcl=-0.25)
    axis(2, at=c(0,0.7), labels=expression(0,italic(C)),las=2)
    axis(3, at=c(0.21,mean(c(0.21,0.77)),0.77), labels=expression(italic(f-W),italic(f),italic(f+W)))
    text(0.5,0.35,expression(italic(2*W)))
    ## text(0.23,0.35,"<")
    text(0.23,0.35,expression(paste(" "%<-%" ")))
    text(0.75,0.35,expression(paste(" "%->%" ")))
    ## text(0.77-0.02,0.35,">")
    ## lines(c(0.23,0.4),c(0.35,0.35),lwd=1.5)
    box(bty="c")
}

### Figure 298

fig_298()

### Figure 304 ###

fig_304 <- function(ts,m_max=60,x_max=50)
{
    C_D <- dB(pgram(ts,center=FALSE)$sdf)
    Im <- compute_cD_and_Im(C_D,length(ts),m_max)$I_m
    plot(1:length(Im),Im,
         xlim=c(0,x_max),xlab=expression(italic(m)),
         ylim=c(0,max(Im)),ylab="estimated MISE",
         typ="p",axes=FALSE,
         main="Figure 304")
    m_min <- which(min(Im) == Im)
    points(m_min,Im[m_min],pch=16)
    axis(1,at=seq(0,x_max,10))
    axis(2,at=seq(0,25,5),las=2)  # HACK
    box(bty="l")
}

### Figure 304

fig_304(ar2_1)

### Figure 305 ###

fig_305 <- function(ts,m,coeffs,innov_var,tag,taper_p=FALSE,right_p=FALSE,m_max=100)
{
    N <- length(ts)
    ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    plot(ar_spec$freqs,dB(ar_spec$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab=paste("AR(",length(coeffs),") spectra  (dB)",sep=""),
         typ="n",axes=FALSE,
         main=paste("Figure 305",tag,sep=""))
    if(taper_p)
    {
        dse <- direct_sdf_est(ts,hanning_taper(N),center=FALSE)
        ugrid <- seq(1,N/2+1,2)
        i_max <- N/2
        N_for_cD_calc <- N/2
    }
    else
    {
        dse <- pgram(ts,center=FALSE)
        ugrid <- seq(1,N/2+1,1)
        i_max <- N/4
        N_for_cD_calc <- N
    }
    C_D <- dB(dse$sdf[ugrid])
    cD_Im <- compute_cD_and_Im(C_D,N_for_cD_calc,m_max)
    lwe <- ts_to_lag_window_sdf_est(ts, m=m, lag_window=parzen_lag_window, center=FALSE, taper=if(taper_p) hanning_taper(N) else rectangular_taper(N))
    lwe_log <- acvs_to_lag_window_sdf_est(cD_Im$c_D[1:i_max], m=m, lag_window=parzen_lag_window, two_sided_p=FALSE)
    if(right_p)
    {
        lines(ar_spec$freqs,dB(ar_spec$sdfe),lwd=0.5)
        lines(lwe_log$freqs,lwe_log$sdfe-digamma(1)*dB(exp(1)),lwd=1.5)
    }
    else
    {
        lines(dse$freqs,dB(dse$sdfe),lwd=0.25,col="gray")
        lines(ar_spec$freqs, dB(ar_spec$sdf), lwd=0.5)
        lines(lwe$freqs, dB(lwe$sdfe), lwd=1)
        lines(lwe_log$freqs, lwe_log$sdfe, lwd=1.5)
    }
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,20),label=FALSE,tcl=-0.25)
    text(0.49,10,tag,pos=2)
    box(bty="l")
}

### Figure 305, top row

fig_305(ar2_1,17,ar2_coeffs,ar2_innov_var,"(a)")
fig_305(ar2_1,17,ar2_coeffs,ar2_innov_var,"(b)",right_p=TRUE)

### Figure 305, bottom row

fig_305(ar4_1,45,ar4_coeffs,ar4_innov_var,"(c)",taper_p=TRUE)
fig_305(ar4_1,45,ar4_coeffs,ar4_innov_var,"(d)",taper_p=TRUE,right_p=TRUE)

### Figure 306 ###

fig_306 <- function(ts,m,coeffs,innov_var,m_max=100)
{
    N <- length(ts)
    ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    plot(ar_spec$freqs,ar_spec$sdfe,
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,50),yaxs="i",ylab="AR(4) spectra  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 306")
    dse <- direct_sdf_est(ts,hanning_taper(N),center=FALSE)
    C_D <- dB(dse$sdf[seq(1,N/2+1,2)])
    cD_Im <- compute_cD_and_Im(C_D,N/2,m_max)
    lwe <- ts_to_lag_window_sdf_est(ts, m=m, lag_window=parzen_lag_window, center=FALSE, taper=hanning_taper(N))
    lwe_log <- acvs_to_lag_window_sdf_est(cD_Im$c_D[1:(N/2)], m=m, lag_window=parzen_lag_window, two_sided_p=FALSE)
    lines(lwe$freqs,lwe$sdfe, lwd=1)
    lines(lwe_log$freqs,convert_from_dB(lwe_log$sdfe-digamma(1)*dB(exp(1))), lwd=1.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,50,10),las=2)
    box(bty="l")
}

### Figure 306

fig_306(ar4_1,45,ar4_coeffs,ar4_innov_var)

### Figure 309 ###

fig_309 <- function(ts,coeffs,innov_var,the_ms=seq(0.01,0.1,0.002),diddle_y=c(0,0))
{
    ts_pgram <- pgram(ts,two_sided=TRUE,start=TRUE,center=FALSE)
    the_gcvs <- sapply(the_ms,function(...)  gcv_dsp(...)$gcv,ts_pgram$sdfe)
    the_best_m <- evaluate_gcv_dse(ts_pgram$sdfe,ar_coeffs_to_sdf(coeffs,innov_var,N_pad=length(ts))$sdf,the_ms,the_measure=mse)$gcv_m
    plot(the_ms,the_gcvs,
         xlab=expression(italic(m)),
         ylim=range(the_gcvs)+diddle_y,ylab=expression(paste("GCV(",italic(m),")")),
         typ="p",axes=FALSE,
         main="Figure 309")
    abline(v=the_best_m,lty="dashed")
    m_min <- which.min(the_gcvs)
    points(the_ms[m_min],the_gcvs[m_min],pch=16)
    axis(1,at=seq(0,0.1,0.02))
    axis(1,at=seq(0,0.1,0.01),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0.59,0.64,0.01),las=2)
    box(bty="l")
}

### Figure 309

fig_309(ar2_1,ar2_coeffs,ar2_innov_var,diddle=c(-0.006,0))

### Figures 310 and 311 ###

fig_310 <- function(ts,m,coeffs,innov_var,tag,hanning_taper_p=FALSE,right_p=FALSE,main="Figure 310")
{
    N <- length(ts)
    taper <- if(hanning_taper_p) hanning_taper(N) else rectangular_taper(N)
    p <- length(coeffs)
    dse <- direct_sdf_est(ts,taper,center=FALSE)
    inc <- if(right_p && p == 4) 2 else 1
    plot(dse$freqs[seq(1,N/2+1,inc)],dB(dse$sdfe[seq(1,N/2+1,inc)]),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-40-10*(p-2),20),yaxs="i",ylab=paste("AR(",p,") spectra  (dB)",sep=""),
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main=main)
    ar_spec <- ar_coeffs_to_sdf(coeffs,innov_var,N_pad=1024)
    lines(ar_spec$freqs,dB(ar_spec$sdf),lwd=0.5)
    smoothed <- if(right_p)
                {
                    if(p==2) ts_to_lag_window_sdf_est(ts,m,lag_window=parzen_lag_window,center=FALSE)
                    else
                    m_to_dsdse(m,dft(direct_sdf_est(ts,taper,two_sided_p=TRUE,start_with_0_p=TRUE,center=FALSE)$sdfe[seq(1,N,2)]),hanning_taper_p)
                }
                else
                    m_to_dsdse(m,dft(direct_sdf_est(ts,taper,two_sided_p=TRUE,start_with_0_p=TRUE,center=FALSE)$sdfe),hanning_taper_p)
    lines(smoothed$freqs,dB(smoothed$sdf))
    cc <- smoothed$cc
    x_left <- 0.025
    x_mid <- x_left + cc$width/2
    y_cc <- -30.0-10*(p-2)
    lines(c(x_left,x_left+cc$width), c(y_cc,y_cc), lwd=0.5)
    lines(c(x_mid,x_mid), c(y_cc-cc$down,y_cc+cc$up), lwd=0.5)
    axis(1,at=seq(0,0.5,0.5))
    axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),label=FALSE,tcl=-0.25)
    text(0.5,10-(p-4)*1.25,tag,pos=2)
    box(bty="l")
}

### Figure 310, left-hand plot

fig_310(ar2_1,0.0286,ar2_coeffs,ar2_innov_var,expression(italic(m == 0.0286)))

### Figure 310, right-hand plot

fig_310(ar2_1,22,ar2_coeffs,ar2_innov_var,expression(italic(m == 22)),right_p=TRUE)

### Figure 311, left-hand plot

fig_310(ar4_1,0.2313,ar4_coeffs,ar4_innov_var,expression(italic(m == 0.2313)),hanning_taper_p=TRUE,main="Figure 311")

### Figure 311, right-hand plot

fig_310(ar4_1,0.2261,ar4_coeffs,ar4_innov_var,expression(italic(m == 0.2261)),hanning_taper_p=TRUE,right_p=TRUE,main="Figure 311")

### Figure 317 ###

fig_317 <- function(ys,y_ats,x_lab)
{
    N <- length(ys)
    plot(0:(N-1),Re(ys),
         xlim=c(0,N),xlab=x_lab,
         ylim=c(y_ats[1],tail(y_ats,1)),ylab=" ",
         typ="n",axes=FALSE,
         main="Figure 317")
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

### Figure 317, top row, left-hand plot

(N <- length(earth_20))        # 20
(M <- next_power_of_2(2*N-1))  # 44
M-N  # 44
tXt <- c(earth_20-mean(earth_20),rep(0,M-N))

fig_317(tXt,seq(-6,6,6),expression(italic(t)))

### Figure 317, top row, right-hand plot

tXt_dft <- dft(tXt)

fig_317(tXt_dft,seq(-32.5,32.5,32.5),expression(italic(k)))

### Figure 317, 2nd row, left-hand plot

tht <- c(hanning_taper(N),rep(0,M-N))

fig_317(tht,c(0,1),expression(italic(t)))

### Figure 317, 2nd row, right-hand plot

tht_dft <- dft(tht)

fig_317(tht_dft,seq(-4,4,4),expression(italic(k)))

### Figure 317, 3rd row, left-hand plot

thttXt <- tht*tXt

fig_317(thttXt,seq(-2,2,2),expression(italic(t)))

### Figure 317, 3rd row, right-hand plot

thttXt_dft <- dft(thttXt)

fig_317(thttXt_dft,seq(-6,6,6),expression(italic(k)))

### Figure 317, 4th row, left-hand plot

tSdk <- abs(thttXt_dft)^2
tsdtau <- Re(inverse_dft(tSdk))

fig_317(tsdtau,seq(-8,8,8),expression(tau))

### Figure 317, 4th row, right-hand plot

fig_317(tSdk,seq(0,60,30),expression(italic(k)))

### Figure 317, 5th row, left-hand plot

temp <- sapply(0:9,parzen_lag_window,10)
twtau <- c(temp,rep(0,M-19),rev(temp[-1]))

fig_317(twtau,c(0,1),expression(tau))

### Figure 317, 5th row, right-hand plot

twtau_dft <- Re(dft(twtau))

fig_317(twtau_dft,seq(0,8,4),expression(italic(k)))

### Figure 317, 6th row, left-hand plot

tslwtau <- twtau*tsdtau

fig_317(tslwtau,seq(-8,8,8),expression(tau))

### Figure 317, 6th row, right-hand plot

tSlwk <- Re(dft(tslwtau))

fig_317(tSlwk,seq(0,60,30),expression(italic(k)))

### Figures 318 and 319 ###

fig_318 <- function(ts,m,lw,tag_1,tag_2,taper=slepian_taper(length(ts),2),dw=NULL,delta_t=1/4,main="Figure 319")
{
    dse <- direct_sdf_est(ts,taper,center=TRUE,delta_t=delta_t)
    plot(dse$freqs,dB(dse$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main=paste(main,tag_2,sep=""))
    lwe <- ts_to_lag_window_sdf_est(ts,m,taper=taper,lag_window=lw,center=TRUE,delta_t=delta_t)
    lines(lwe$freqs, dB(lwe$sdfe))
    cc <- lwe$cc
    x_cc <- 0.16
    y_cc <- 20
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
    sw <- lag_window_to_smoothing_window(sapply(0:(length(ts)-1),lw,m),N_pad=8192,delta_t=delta_t)
    lines(c(0.25-sw$B_W/2,0.25+sw$B_W/2),c(-9,-9),lwd=0.5)
    if(is.null(dw))
    {
        freqs_two_sided <- c(rev(-sw$freqs[-1]),sw$freqs)
        dw_or_sw_two_sided <- careful_dB(c(rev(sw$sw[-1]),sw$sw))
    }
    else
    {
        temp <- dw(m,delta_t=delta_t)
        freqs_two_sided <- c(rev(-temp$freqs[-1]),temp$freqs)
        dw_or_sw_two_sided <- careful_dB(c(rev(temp$dw[-1]),temp$dw))
    }
    lines(freqs_two_sided+0.25,dw_or_sw_two_sided-max(dw_or_sw_two_sided)-14,lwd=0.5)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.0,80,tag_1,pos=1)
    text(1.9,80,tag_2,pos=1)
    box(bty="l")
}

### Figure 318

fig_318(ocean_wave,150,parzen_lag_window,expression(paste("Parzen, ", italic(m==150)))," ",main="Figure 318")

### Figure 319, upper plot

fig_318(ocean_wave,55,parzen_lag_window,expression(paste("Parzen, ", italic(m==55))),"(a)")

### Figure 319, 2nd plot

fig_318(ocean_wave,29.749,daniell_lag_window,expression(paste("Daniell, ", italic(m==29.749))),dw=daniell_design_window,"(b)")

### Figure 319, 3rd plot

fig_318(ocean_wave,24.717,bartlett_priestley_lag_window,expression(paste("Bartlett-Priestley, ", italic(m==24.717))),dw=bartlett_priestley_design_window,"(c)")

### Figure 319, bottom plot

fig_318(ocean_wave,23.666,gaussian_lag_window,expression(paste("Gaussian, ", italic(m==23.666))),"(d)")

### Figure 320 ###

fig_320 <- function(ts,ms,lws,taper=slepian_taper(length(ts),2),delta_t=1/4)
{
    lwe_1 <- ts_to_lag_window_sdf_est(ts,ms[1],taper=taper,lag_window=lws[[1]],center=TRUE,delta_t=delta_t)
    lwe_2 <- ts_to_lag_window_sdf_est(ts,ms[2],taper=taper,lag_window=lws[[2]],center=TRUE,delta_t=delta_t)
    lwe_3 <- ts_to_lag_window_sdf_est(ts,ms[3],taper=taper,lag_window=lws[[3]],center=TRUE,delta_t=delta_t)
    plot(lwe_2$freqs[-1],diff(dB(lwe_2$sdfe)),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-0.55,0.55),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray40",axes=FALSE,
         main="Figure 320")
    lines(lwe_3$freqs[-1],diff(dB(lwe_3$sdfe)))
    lines(lwe_1$freqs[-1],diff(dB(lwe_1$sdfe)),lty="dashed")
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-0.5,0.5,0.5),las=2)
    axis(2,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 320

fig_320(ocean_wave,c(55,24.717,23.666),list(parzen_lag_window,bartlett_priestley_lag_window,gaussian_lag_window))

### Figure 322 ###

fig_322 <- function(ts,y_ats,tag,delta_t=1.7712)
{
    N <- length(ts)
    plot((0:(N-1))*delta_t,ts,
         xlim=c(0,2000),xlab="distance  (m)",
         ylim=c(y_ats[1],y_ats[length(y_ats)]),ylab="height  (feet)",
         typ="l",axes=FALSE,
         main="Figure 322")
    axis(1,at=seq(0,2000,500))
    axis(2,at=y_ats,las=2)
    text(2000,y_ats[length(y_ats)],tag,pos=2)
    box(bty="l")
}

### Figure 322, top plot

fig_322(rough,seq(-50,0,10),"Rough")

### Figure 322, middle plot

fig_322(smooth,seq(-20,0,10),"Smooth")

### Figure 322, bottom plot

fig_322(smooth,seq(-12,-7,1),"Smooth")

### Figure 323 ###

fig_323 <- function(ts,m,tag_1,tag_2,parzen_smooth_p=TRUE,log_smooth_p=FALSE,m_max=60,delta_t=1.7712)
{
    N <- length(ts)
    the_pgram <- pgram(ts,center=TRUE,delta_t=delta_t)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,1/(2*delta_t)),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/meter)")),
         ylim=c(-40,40),yaxs="i",ylab="dB",
         typ="l",lwd=0.25, col="gray",axes=FALSE,
         main=paste("Figure 323",tag_2,sep=""))
    if(log_smooth_p)
    {
        temp <- dB(the_pgram$sdfe[-1])
        cD_Im_r <- compute_cD_and_Im(c(temp[1],temp),N,m_max,delta_t=delta_t)
        log_pgram_sm <- acvs_to_lag_window_sdf_est(cD_Im_r$c_D[1:ceiling((N+1)/2)], m, lag_window=parzen_lag_window, two_sided_p=FALSE, delta_t=delta_t)
        lines(log_pgram_sm$freqs, log_pgram_sm$sdfe, lwd=0.5)
    }
    pgram_sm <- if(parzen_smooth_p) ts_to_lag_window_sdf_est(ts,m,lag_window=parzen_lag_window,center=TRUE,delta_t=delta_t) else gcv_dsp(m,pgram(ts,two_sided=TRUE,start=TRUE,delta_t=delta_t,center=TRUE)$sdfe,delta_t=delta_t,keep_0_p=FALSE)
    lines(pgram_sm$freqs, dB(pgram_sm$sdfe), lwd=1.5)
    x_cc <-  0.06
    y_cc <- -30.0
    lines(c(x_cc,x_cc),y_cc+c(pgram_sm$cc$up,-pgram_sm$cc$down),lwd=0.5)
    lines(x_cc+c(-pgram_sm$cc$width/2,pgram_sm$cc$width/2),c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,0.3,0.1))
    axis(1,at=seq(0,0.3,0.05),label=FALSE)
    axis(1,at=seq(0,0.3,0.01),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,40,20),las=2)
    axis(2,at=seq(-40,40,10),label=FALSE,tcl=-0.25)
    text(0.25/delta_t,40,tag_1,pos=1)
    text(0.49/delta_t,30,tag_2,pos=2)
    box(bty="l")
}

### Figure 323, left-hand column

fig_323(rough,88,"Rough","(a)",delta_t=1.7712)
fig_323(rough,35,"Rough","(c)",log_smooth_p=TRUE,delta_t=1.7712)
fig_323(rough,0.0459,"Rough","(e)",parzen_smooth_p=FALSE,delta_t=1.7712)
    
### Figure 323, right-hand column

fig_323(smooth,20,"Smooth","(b)",delta_t=1.7712)
fig_323(smooth,11,"Smooth","(d)",log_smooth_p=TRUE,delta_t=1.7712)
fig_323(smooth,0.0712,"Smooth","(f)",parzen_smooth_p=FALSE,delta_t=1.7712)

### Table 325, 1st column of numbers

temp <- ts_to_lag_window_sdf_est(rough,88,lag_window=parzen_lag_window,center=TRUE,delta_t=1.7712)
round(temp$cc$width,4)            # 0.0122
round(temp$cc$up+temp$cc$down,2)  # 3.54

### Table 325, 2nd column

temp <- ts_to_lag_window_sdf_est(rough,35,lag_window=parzen_lag_window,center=TRUE,delta_t=1.7712)
round(temp$cc$width,4)            # 0.0302
round(temp$cc$up+temp$cc$down,2)  # 2.22

### Table 325, 3rd column

temp <- gcv_dsp(0.0459,pgram(rough,two_sided=TRUE,start=TRUE,delta_t=1.7712,center=TRUE)$sdfe,delta_t=1.7712,keep_0_p=FALSE)
round(temp$cc$width,4)            # 0.0278
round(temp$cc$up+temp$cc$down,2)  # 2.31

### Table 325, 4th column

temp <- ts_to_lag_window_sdf_est(smooth,20,lag_window=parzen_lag_window,center=TRUE,delta_t=1.7712)
round(temp$cc$width,4)            # 0.0536
round(temp$cc$up+temp$cc$down,2)  # 3.33

### Table 325, 5th column

temp <- ts_to_lag_window_sdf_est(smooth,11,lag_window=parzen_lag_window,center=TRUE,delta_t=1.7712)
round(temp$cc$width,4)            # 0.0964
round(temp$cc$up+temp$cc$down,2)  # 2.46

### Table 325, 6th column

temp <- gcv_dsp(0.0712,pgram(smooth,two_sided=TRUE,start=TRUE,delta_t=1.7712,center=TRUE)$sdfe,delta_t=1.7712,keep_0_p=FALSE)
round(temp$cc$width,4)            # 0.0703
round(temp$cc$up+temp$cc$down,2)  # 2.89

### Figure 325 ###

fig_325 <- function(ts_r,m_r,ts_s,m_s,dB_shift=13,delta_t=1.7712)
{
    parzen_r <- ts_to_lag_window_sdf_est(ts_r,m_r,lag_window=parzen_lag_window,center=TRUE,delta_t=delta_t)
    parzen_s <- ts_to_lag_window_sdf_est(ts_s,m_s,lag_window=parzen_lag_window,center=TRUE,delta_t=delta_t)
    plot(parzen_r$freqs,dB(parzen_r$sdfe),
         xlim=c(0,1/(2*delta_t)),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/meter)")),
         ylim=c(-3,35),yaxs="i",ylab="dB",
         typ="l",lwd=1.5,axes=FALSE,
         main="Figure 325")
    lines(parzen_s$freqs, dB(parzen_s$sdfe)+dB_shift, lwd=0.5)
    j <- which.min(abs(parzen_r$freqs-0.01))
    lines(rep(parzen_r$freqs[j],2), dB(parzen_r$sdfe[j])+c(parzen_r$cc$up,-parzen_r$cc$down),lwd=0.5)
    j <- which.min(abs(parzen_s$freqs-0.01))
    lines(rep(parzen_s$freqs[j],2), dB(parzen_s$sdfe[j])+dB_shift+c(parzen_s$cc$up,-parzen_s$cc$down),lwd=0.5)
    axis(1,at=seq(0,0.3,0.1))
    axis(1,at=seq(0,0.3,0.05),label=FALSE)
    axis(1,at=seq(0,0.3,0.01),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,30,10),las=2)
    axis(2,at=seq(-10,40,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 325

fig_325(rough,35,smooth,11)

### Figure 326 ###

fig_326 <- function(ts,y_lim,tag,y_lab,y_maj_ats,y_min_ats,y_text_lab)
{
    N <- length(ts)
    plot((0:(N-1)),ts,
         xlim=c(0,4000),xlab="time  (min)",
         ylim=y_lim,ylab=y_lab,
         typ="l",axes=FALSE,
         main=paste("Figure 326",tag,sep=""))
    axis(1,at=seq(0,4000,1000))
    axis(2,at=y_maj_ats,las=2)
    axis(2,at=y_min_ats,label=FALSE,tcl=-0.25)
    text(4010,y_text_lab,tag,pos=2)
    box(bty="l")
}

### Figure 326, top plot

fig_326(ac_time_differences,c(-22.5,-20.5),"(a)","time difference",c(-22,-21),seq(-23,-20,0.5),-20.72)

### Figure 326, bottom plot

fig_326(ac_fractional_frequencies,c(-0.75,0.75),"(b)","fractional frequency",seq(-0.5,0.5,0.5),seq(-0.5,0.5,0.5),0.67)

### Figure 327 ###

fig_327 <- function(ts,taper,tag_1,tag_2,log_f_p=FALSE)
{
    dse <- direct_sdf_est(ts,taper(length(ts)),center=TRUE)
    plot(if(log_f_p) dB(dse$freqs[-1]) else dse$freqs[-1],dB(dse$sdfe[-1]),
         xlim=if(log_f_p) range(dB(dse$freqs[-1])) else c(0,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/mimute)")),
         ylim=c(-80,30),yaxs="i",ylab="dB",typ="l",axes=FALSE,
         main=paste("Figure 327",tag_2,sep=""))
    if(log_f_p)
    {
        xs <- c(-28,-8)
        lines(xs,-40-2*xs,lty="solid",col="gray",lwd=1.25)
        ttn <- 2:9
        axis(1, at=dB(c(1/100000,1/10000,1/1000,1/100,1/10,1)), labels=expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
        axis(1, at=dB(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)), labels=FALSE, tcl=-0.25)
    }
    else
    {
        axis(1,at=seq(0,0.5,0.5))
        axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    }
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,30,10),label=FALSE,tcl=-0.25)
    text(0.25,24.25,tag_1,pos=1)
    text(if(log_f_p) dB(0.49) else 0.49,16.25,tag_2,pos=2)
    box(bty="l")
}

### Figure 327, top row

fig_327(ac_time_differences,default_taper,"periodogram","(a)")
fig_327(ac_time_differences,default_taper," ","(b)",log_f_p=TRUE)

### Figure 327, bottom row

fig_327(ac_time_differences,hanning_taper,"Hanning","(c)")
fig_327(ac_time_differences,hanning_taper," ","(d)",log_f_p=TRUE)

### Figure 328 ###

fig_328 <- function(ts,m,tag,log_f_p=FALSE)
{
    the_pgram <- pgram(ts,center=TRUE)
    plot(if(log_f_p) dB(the_pgram$freqs[-1]) else the_pgram$freqs[-1],dB(the_pgram$sdfe[-1]),
         xlim=if(log_f_p) range(dB(the_pgram$freqs[-1])) else c(0,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/mimute)")),
         ylim=c(-50,0),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main=paste("Figure 328",tag,sep=""))
    lwe <- ts_to_lag_window_sdf_est(ts,m,lag_window=gaussian_lag_window,center=TRUE)
    lines(if(log_f_p) dB(lwe$freqs[-1]) else lwe$freqs[-1],dB(lwe$sdfe[-1]))
    if(log_f_p)
    {
        ttn <- 2:9
        axis(1, at=dB(c(1/100000,1/10000,1/1000,1/100,1/10,1)), labels=expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
        axis(1, at=dB(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)), labels=FALSE, tcl=-0.25)
    }
    else
    {
        x_cc <- 0.3
        y_cc <- -45.0
        lines(c(x_cc,x_cc),y_cc+c(lwe$cc$up,-lwe$cc$down),lwd=0.5)
        lines(x_cc+c(-lwe$cc$width/2,lwe$cc$width/2),c(y_cc,y_cc),lwd=0.5)
        axis(1,at=seq(0,0.5,0.5))
        axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    }
    axis(2,at=seq(-50,0,10),las=2)
    text(if(log_f_p) -35.4 else 0.01,-3.125,tag,pos=4)
    box(bty="l")
}

### Figure 328, top row

fig_328(ac_fractional_frequencies,80,"(a)")
fig_328(ac_fractional_frequencies,80,"(b)",log_f_p=TRUE)

### Figure 328, bottom row

fig_328(ac_fractional_frequencies,20,"(c)")
fig_328(ac_fractional_frequencies,20,"(d)",log_f_p=TRUE)

### Figure 329 ###

fig_329 <- function(ts_td,ts_ff,m,tag)
{
    dse <- direct_sdf_est(ts_td,hanning_taper(length(ts_td)),center=TRUE)
    plot(dB(dse$freqs[-1]),dB(dse$sdfe[-1]),
         xlim=range(dB(dse$freqs[-1])),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/mimute)")),
         ylim=c(-80,30),yaxs="i",ylab="dB",
         typ="l",lwd=0.25,col="gray",axes=FALSE,
         main=paste("Figure 329",tag,sep=""))
    lwe <- ts_to_lag_window_sdf_est(ts_ff,m,lag_window=gaussian_lag_window,center=TRUE)
    lines(dB(lwe$freqs[-1]),dB(36*lwe$sdfe[-1]/(10000*4*(sin(pi*lwe$freq[-1]))^2)))
    axis(1, at=dB(c(1/100000,1/10000,1/1000,1/100,1/10,1)), labels=expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0))
    ttn <- 2:9
    axis(1, at=dB(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)), labels=FALSE, tcl=-0.25)
    axis(2,at=seq(-80,20,20),las=2)
    axis(2,at=seq(-80,30,10),label=FALSE,tcl=-0.25)
    text(dB(0.49), 16.25,tag,pos=2)
    box(bty="l")
}

### Figure 329, left-hand plot

fig_329(ac_time_differences,ac_fractional_frequencies,80,"(a)")

### Figure 329, right-hand plot

fig_329(ac_time_differences,ac_fractional_frequencies,20,"(b)")

### Figure 330 ###

fig_330 <- function(ts)
{
    N <- length(ts)
    plot((0:(N-1)),ts,
         xlab="time  (sec)",
         ylim=c(2,12),ylab="altitude  (meters)",
         typ="l",axes=FALSE,
         main="Figure 330")
    axis(1,at=seq(0,2000,500))
    axis(2,at=seq(2,12,5),las=2)
    axis(2,at=seq(2,12,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 330

fig_330(ship_altitude)

### Figure 331 ###

fig_331_top <- function(ts,taper,tag_1,tag_2,lwd=0.5,col="gray")
{
    dse <- direct_sdf_est(ts,taper(length(ts)),center=TRUE)
    plot(dse$freqs,dB(dse$sdfe),
         xlim=c(-0.006,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-70,50),ylab="dB",
         typ="l",lwd=lwd,col=col,axes=FALSE,
         main=paste("Figure 331",tag_2,sep=""))
    cc_xy <- c(0.425,20)
    cc <- dse$cc
    lines(c(cc_xy[1],cc_xy[1]),cc_xy[2]+c(cc$up,-cc$down),lwd=0.5,col=col)
    lines(cc_xy[1]+c(-cc$width/2,cc$width/2),c(cc_xy[2],cc_xy[2]),lwd=0.5,col=col)
    axis(1,at=seq(0,2.0,0.1))
    axis(2,at=seq(-60,40,20),las=2)
    axis(2,at=seq(-70,50,10),label=FALSE,tcl=-0.25)
    text(0.25,50.7,tag_1,pos=1)
    text(0.5,41.2,tag_2,pos=2)
    box(bty="l")
}

fig_331_bot <- function(ts,m,lag_window,tag_1,tag_2)
{
    lwe_pgram <- ts_to_lag_window_sdf_est(ts,m,taper=rectangular_taper,lag_window=lag_window,center=TRUE)
    lwe_hanning <- ts_to_lag_window_sdf_est(ts,m,taper=hanning_taper,lag_window=lag_window,center=TRUE)
    plot(lwe_pgram$freqs,dB(lwe_pgram$sdfe),
         xlim=c(-0.006,0.5),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,22),ylab="dB",
         typ="l",col="gray",axes=FALSE,
         main=paste("Figure 331",tag_2,sep=""))
    cc_xy <- c(0.425,10)
    cc <- lwe_pgram$cc
    lines(c(cc_xy[1],cc_xy[1]),cc_xy[2]+c(cc$up,-cc$down),lwd=0.5,col="gray")
    lines(cc_xy[1]+c(-cc$width/2,cc$width/2),c(cc_xy[2],cc_xy[2]),lwd=0.5,col="gray")
    ##
    lines(lwe_hanning$freqs,dB(lwe_hanning$sdfe))
    cc_xy <- c(0.425,0)
    cc <- lwe_hanning$cc
    lines(c(cc_xy[1],cc_xy[1]),cc_xy[2]+c(cc$up,-cc$down),lwd=0.5)
    lines(cc_xy[1]+c(-cc$width/2,cc$width/2),c(cc_xy[2],cc_xy[2]),lwd=0.5)
    ##
    axis(1,at=seq(0,2.0,0.1))
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(0.25,22.7,tag_1,pos=1)
    text(0.5,17.8,tag_2,pos=2)
    box(bty="l")
}

### Figure 331, top row

fig_331_top(ship_altitude,rectangular_taper,"periodogram","(a)")
fig_331_top(ship_altitude,hanning_taper,"Hanning","(b)",lwd=1,col="black")

### Figure 331, bottom row

fig_331_bot(ship_altitude,80,gaussian_lag_window,"Gaussian","(c)")
fig_331_bot(ship_altitude,185,parzen_lag_window,"Parzen","(d)")

### Figure 333 ###

fig_333 <- function(ts,tag)
{
    N <- length(ts)
    years <- (0:(N-1))/12 + 1958 + 5/24  # first data point is for March 1958
    if(tag == "(a)")
    {
        xs <- years
        ys_light <- ts
        x_reg <- 1:length(ts)
        x2_reg <- x_reg^2
        ys_dark <- fitted(lm(ts ~ x_reg + x2_reg))
        y_lab <- expression(paste(CO[2],"  (ppmv)",sep=""))
        y_ats_big <- seq(320,400,20)
        y_ats_little <- NULL
    }
    else
        if(tag == "(b)")
        {
            xs <- years[-(1:12)]
            ys_light <- diff(ts,lag=12)
            ys_dark <- fitted(lm(ys_light ~ xs))
            y_lab <- "seasonal diff"
            y_ats_big <- 0:4
            y_ats_little <- NULL
        }
    else
    {
        xs <- years[-(1:13)]
        ys_light <- diff(diff(ts,lag=12))
        ys_dark <- NULL
        y_lab <- "seasonal & 1st diffs"
        y_ats_big <- seq(-1,1,1)
        y_ats_little <- seq(-1.5,1.5,0.5)
    }
    plot(xs,ys_light,
         xlim=range(years),xlab="years",
         ylab=y_lab,
         typ="l",col="darkgray",axes=FALSE,
         main=paste("Figure 333",tag,sep=""))
    if(!is.null(ys_dark)) lines(xs,ys_dark)
    axis(1, at=seq(1960,2020,10))
    axis(1, at=seq(1950,2020,2), labels=FALSE, tcl=-0.25)
    axis(2,at=y_ats_big,las=2)
    axis(2,at=y_ats_little,label=FALSE,tcl=-0.25)
    temp <- range(ys_light)
    text(1961.2,temp[1]+diff(temp)*7/8,tag,pos=2)
    box(bty="l")
}


### Figure 333, top, middle and bottom plots

fig_333(co2,"(a)")
fig_333(co2,"(b)")
fig_333(co2,"(c)")

### Figure 335 ###

fig_335a<- function()
{
    freqs     <- (0:(6*1024))/1024
    G1     <- sapply(freqs,function(f) 4*(sin(pi*f/12))^2)
    G12    <- sapply(freqs,function(f) 4*(sin(pi*f))^2)
    G1_G12 <- G1 * G12
    plot(freqs,dB(G1),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-50,10),ylab="dB",
         typ="l",lty="dashed",axes=FALSE,
         main="Figure 335(a)")
    lines(freqs, careful_dB(G12))
    lines(freqs, careful_dB(G1_G12), col="darkgray")
    ##
    axis(1,at=0:6)
    axis(2,at=seq(-50,10,10),las=2)
    text(6,-20,"(a)",pos=2)
    box(bty="l")
}

fig_335b <- function(ts)
{
    ys <- diff(diff(ts,lag=12))
    the_pgram <- pgram(ys,center=TRUE,delta_t=1/12)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-50,10),ylab="dB",
         typ="l",col="darkgray",axes=FALSE,
         main="Figure 335(b)")
    lwe <- ts_to_lag_window_sdf_est(ys,185,lag_window=parzen_lag_window,delta_t=1/12,center=TRUE)
    lines(lwe$freqs,dB(lwe$sdfe))
    cc <- lwe$cc
    cc_x <- 1
    cc_y <- -5
    lines(c(cc_x,cc_x),cc_y+c(cc$up,-cc$down),lwd=0.5)
    lines(cc_x+c(-cc$width/2,cc$width/2),c(cc_y,cc_y),lwd=0.5)
    ##
    axis(1,at=0:6)
    axis(2,at=seq(-50,10,10),las=2)
    text(6,5,"(b)",pos=2)
    box(bty="l")
}

fig_335c_or_d <- function(ts,tag,right_p=FALSE)
{
    ys <- diff(diff(ts,lag=12))
    ys_bar <- mean(ys)
    N <- length(ys)
    var_Y <- sum((ys-ys_bar)^2)/(4*N)
    freqs     <- (0:(6*1024))/1024
    G1     <- sapply(freqs,function(f) 4*(sin(pi*f/12))^2)
    G12    <- sapply(freqs,function(f) 4*(sin(pi*f))^2)
    sq <- if(right_p)
          {
              phi <- 2*sum((ys[-1]-ys_bar)*(ys[-N]-ys_bar))/sum((ys-ys_bar)^2) + 1
              i_var <- var_Y*(1+phi)
              G1*G12*sapply(freqs,function(f) i_var/(12*(1 + phi^2 -2*phi*cos(pi*f/6))))
          }
          else G1*G12*var_Y/12
    plot(freqs,dB(sq),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-40,-10),ylab="dB",
         typ="l",col="darkgray",axes=FALSE,
         main=paste("Figure 335",tag,sep=""))
    lwe <- ts_to_lag_window_sdf_est(ys,185,lag_window=parzen_lag_window,delta_t=1/12,center=TRUE)
    lines(lwe$freqs,dB(lwe$sdfe))
    cc <- lwe$cc
    cc_x <- 1
    cc_y <- -15
    lines(c(cc_x,cc_x),cc_y+c(cc$up,-cc$down),lwd=0.5)
    lines(cc_x+c(-cc$width/2,cc$width/2),c(cc_y,cc_y),lwd=0.5)
    axis(1,at=0:6)
    axis(2,at=seq(-50,10,10),las=2)
    text(6,5,"(b)",pos=2)
    text(6,-10.25,tag,pos=2)
    box(bty="l")
}

### Figure 335, top row

fig_335a()
fig_335b(co2)

### Figure 335, bottom row

fig_335c_or_d(co2,"(c)")
fig_335c_or_d(co2,"(d)",right_p=TRUE)

### Figures 338a and 338b ###

fig_338a <- function(freqs,sdf,tag,log_f_p=FALSE,dB_sdf_p=FALSE)
{
    if(log_f_p)
    {
        xs <- log10(freqs[-1])
        x_lim <- range(xs)
        x_label <- expression(10^-5,10^-4,10^-3,10^-2,10^-1,10^0)
        x_ats_big <- log10(c(1/100000,1/10000,1/1000,1/100,1/10,1))
        ttn <- 2:9
        x_ats_little <- log10(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10))
        temp <- range(xs)
        x_text <- temp[1] + diff(temp) * 0.98
    }
    else
    {
        xs <- freqs[-1]
        x_lim <- c(0,0.5)
        x_label <- TRUE
        x_ats_big <- c(0,0.5)
        x_ats_little <- seq(0,0.5,0.1)
        x_text <- 0.49
    }
    if(dB_sdf_p)
    {
        ys <- dB(sdf[-1])
        y_lim <- c(-15,25)
        y_ats <- seq(-10,20,10)
        y_lab <- "SDF  (dB)"
        y_text <- 20
    }
    else
    {
        ys <- sdf[-1]
        y_lim <- c(0,105)
        y_ats <- seq(0,100,20)
        y_lab <- "SDF"
        y_text <- 91.875
    }
    plot(xs,ys,
         xlim=x_lim,xaxs="i",xlab=expression(italic(f)),
         ylim=y_lim,yaxs="i",ylab=y_lab,
         typ="l",axes=FALSE,
         main=paste("Figure 338a",tag,sep=""))
    axis(1,at=x_ats_big,label=x_label)
    axis(1,at=x_ats_little,label=FALSE,tcl=-0.25)
    axis(2,at=y_ats,las=2)
    text(x_text,y_text,tag,pos=2)
    box(bty="l")
}

fig_338b <- function(freqs,sdf)
{
    xs_vps <- log10(freqs[-1])
    ys_vps <- freqs[-1]*(sdf[-1])/log10(exp(1))
    plot(xs_vps,ys_vps,
         xaxs="i",xlab=expression(log[10](italic(f))),
         ylim=c(0,max(ys_vps)),yaxs="i",ylab="variance preserving SDF",
         typ="l",axes=FALSE,
         main="Figure 338b")
    axis(1, at=log10(c(1/10000,1/1000,1/100,1/10)),labels=(-4):(-1))
    ttn <- 2:9
    axis(1, at=log10(c(ttn/100000,ttn/10000,ttn/1000,ttn/100,ttn/10)),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(0,12,2),las=2)
    axis(2,at=seq(0,12,1),labels=FALSE,tcl=-0.25)
    box(bty="l")
}

N_pad <- 4096
ar4_innov_var <- 0.002
ar4_coeffs <- c(2.7607, -3.8106, 2.6535, -0.9238)
ar4_sdf <- ar_coeffs_to_sdf(ar4_coeffs,ar4_innov_var,N_pad=N_pad)

ar2_innov_var <- 0.0004
ar2_coeffs <- c(-1.15, -0.96)
ar2_sdf <- ar_coeffs_to_sdf(ar2_coeffs,ar2_innov_var,N_pad=N_pad)

ar1_innov_var <- 0.25
ar1_coeffs <- 0.95
ar1_sdf <- ar_coeffs_to_sdf(ar1_coeffs,ar1_innov_var,N_pad=N_pad)
the_sdf <- ar1_sdf$sdf + ar2_sdf$sdf + ar4_sdf$sdf

### Figure 338a, top row

fig_338a(ar4_sdf$freqs,the_sdf,"(a)",log_f=FALSE,dB_sdf=FALSE)
fig_338a(ar4_sdf$freqs,the_sdf,"(b)",log_f=TRUE,dB_sdf=FALSE)

### Figure 338a, bottom row

fig_338a(ar4_sdf$freqs,the_sdf,"(c)",log_f=FALSE,dB_sdf=TRUE)
fig_338a(ar4_sdf$freqs,the_sdf,"(d)",log_f=TRUE,dB_sdf=TRUE)

### Figure 338b

fig_338b(ar4_sdf$freqs,the_sdf)

### Figure 346 ###

fig_346 <- function()
{
    freqs <- seq(0,1.0,0.001)
    hypo_sdf <- sapply(freqs,function(f) 0.25 + 225/(2850 + 3960*cos(pi*f) + 1600 * cos(2*pi*f)))
    plot(freqs,hypo_sdf,
         xlim=c(0,1),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,10),yaxs="i",ylab="SDF",
         typ="l",axes=FALSE,
         main="Figure 346")
    axis(1, at=seq(0,1,0.2))
    axis(1, at=seq(0,1,0.1), labels=FALSE, tcl=-0.25)
    axis(1, at=seq(0,1,0.02), labels=FALSE, tcl=-0.125)
    ##
    axis(2, at=seq(0,10,2), las=2)
    axis(2, at=seq(0,10,1), labels=FALSE, tcl=-0.25)
    axis(2, at=seq(0,10,0.2), labels=FALSE, tcl=-0.125)
    ##
    axis(3, at=seq(0,1,0.2), labels=FALSE)
    axis(3, at=seq(0,1,0.1), labels=FALSE, tcl=-0.25)
    axis(3, at=seq(0,1,0.02), labels=FALSE, tcl=-0.125)
    ##
    axis(4, at=seq(0,10,2), labels=FALSE)
    axis(4, at=seq(0,10,1), labels=FALSE, tcl=-0.25)
    axis(4, at=seq(0,10,0.2), labels=FALSE, tcl=-0.125)
    box(bty="l")
}


### Figure 346

fig_346()
