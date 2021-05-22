### R CODE FOR REPRODUCING CONTENT OF FIGURES AND TABLES IN CHAPTER 11 ...

ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean_wave.txt")
ship_altitude <- scan("http://faculty.washington.edu/dbp/sauts/Data/ship_altitude_2048.txt")
ac_time_differences <- scan("http://faculty.washington.edu/dbp/sauts/Data/maser_deglitched.txt")
ac_fractional_frequencies <- diff(ac_time_differences)*100/6

### functions used to compute content of figures in Chapter 11 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/burg_algorithm.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H_bar.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/cosine_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_lwe.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/gen_cis_for_ar_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_odd.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_windows.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/mt_innov_var.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/mt_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/sim_ar_process.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ts_to_lag_window_sdf_est.R")

### NOTE: to install the sapa library, uncomment the following three statements
###       and execute them:
###
### install.packages("devtools")
### devtools::install_github("wconstan/ifultools")
### devtools::install_github("wconstan/sapa")

library(sapa)

### AR(2) innovation variance and coefficients

ar2_innov_var <- 1
ar2_coeffs    <- c(0.75,-0.5)

### compute ACVS at lags 0 to 1024 using ar_coeffs_to_acvs

ar2_acvs_0_to_1024 <- ar_coeffs_to_acvs(ar2_coeffs,1024,var=ar2_innov_var,process=FALSE)

### compute ACVS at lags 0 to 1024 using Equation 508a

eqn_508a <- c(16/9,8/9)
all.equal(ar2_acvs_0_to_1024[1:2],eqn_508a)  # TRUE
alt_ar2_acvs_0_to_1024 <- rep(0,1025)
alt_ar2_acvs_0_to_1024[1:2] <- eqn_508a
for(i in 3:1025) alt_ar2_acvs_0_to_1024[i] <- sum(ar2_coeffs*alt_ar2_acvs_0_to_1024[(i-1):(i-2)])

### check that two ways of getting ACVS agree

all.equal(ar2_acvs_0_to_1024,alt_ar2_acvs_0_to_1024)  # TRUE

### AR(4) innovation variance and coefficients

ar4_innov_var <- 0.002
ar4_coeffs <- c(2.7607, -3.8106, 2.6535, -0.9238)

### compute ACVS at lags 0 to 256 using ar_coeffs_to_acvs

ar4_acvs_0_to_256 <- ar_coeffs_to_acvs(ar4_coeffs,256,var=ar4_innov_var,process=FALSE)

### compute ACVS at lags 0 to 256 using Equation 508b

eqn_508b <- c(1.523434580,1.091506153,0.054284646,-0.975329449)
all.equal(ar4_acvs_0_to_256[1:4],eqn_508b)  # TRUE
alt_ar4_acvs_0_to_256 <- rep(0,257)
alt_ar4_acvs_0_to_256[1:4] <- eqn_508b
for(i in 5:257) alt_ar4_acvs_0_to_256[i] <- sum(ar4_coeffs*alt_ar4_acvs_0_to_256[(i-1):(i-4)])

### check that two ways of getting ACVS agree

all.equal(ar4_acvs_0_to_256,alt_ar4_acvs_0_to_256)  # TRUE

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 603 ###

fig_603 <- function(ys,tag,y_lab,y_lim=c(-5,5))
{
    M <- length(ys)
    xs_1 <- 0:((M/2)-1)
    xs_2 <- ((M/2)-1):(M-1)
    plot(xs_1,ys[1:(M/2)],
         xlim=c(32,M-32),xlab=if(tag=="(a)") expression(italic(k)) else expression(italic(t)),
         ylim=y_lim,ylab=y_lab,
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 603",tag,sep=""))
    lines(xs_2,ys[(M/2):M],lwd=0.25,col=if(tag=="(a)") "black" else "gray75")
    axis(1,at=seq(0,M,M/4))
    axis(1,at=seq(0,M,M/32),label=FALSE,tcl=-0.25)
    axis(2,at=seq(round(y_lim[1]),round(y_lim[2]),5),las=2)
    axis(2,at=seq(round(y_lim[1]),round(y_lim[2]),1),label=FALSE,tcl=-0.25)
    text(2048,0.9*diff(y_lim)+y_lim[1],tag,pos=2)
    box(bty="l")
}

N <- 1024
circularized_ar2_acvs <- c(ar2_acvs_0_to_1024,rev(ar2_acvs_0_to_1024[-c(1,N+1)]))

S_k_AR2 <- Re(dft(circularized_ar2_acvs))

set.seed(42)
cal_Y_k <- sqrt(S_k_AR2/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
Y_t <- dft(cal_Y_k)
ts_Re <- Re(Y_t)
ts_Im <- Im(Y_t)

### Figure 603, plots from top to bottom

fig_603(S_k_AR2,"(a)",expression(italic(S[k])),y_lim=round(range(S_k_AR2)))
fig_603(ts_Re,"(b)",expression(paste("Re(",italic(Y[t]),")")))
fig_603(ts_Im,"(c)",expression(paste("Im(",italic(Y[t]),")")))

### Figure 607(a) ###

fig_607 <- function(N_prime,tag)
{
    taus <- 0:256
    plot(taus,ar4_acvs_0_to_256,
         xlab=expression(tau),
         ylim=c(-2.5,2.5),ylab="ACVSs",
         pch=16,cex=0.5,axes=FALSE,
         main=paste("Figure 607",tag,sep=""))
    lines(rep(gen_s_prime(N_prime),5)[1:257])
    abline(v=63,lty="dashed")
    axis(1,at=seq(0,256,64))
    axis(1,at=seq(0,256,16),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-2,2,2),las=2)
    axis(2,at=seq(-2,2,1),label=FALSE,tcl=-0.25)
    text(249,2.1,tag,pos=2)
    box(bty="l")
}


gen_s_prime <- function(N_prime,phi=ar4_coeffs,sig2=ar4_innov_var)
{
    M <- 2*N_prime
    varphi <- rep(0,M)
    varphi[1] <- -1
    varphi[2:5] <- phi
    return(Re(inverse_dft(sig2/abs(dft(varphi))^2)))
}

### Figure 607(a)

fig_607(32,"(a)")

### Figure 607(b)

fig_607(64,"(b)")

### Figure 607(c)

fig_607(128,"(c)")

### Table 607, first row

table_607_first_row <- function(N_prime,ar_acvs=ar4_acvs_0_to_256,N=64)
{
    return(sqrt(sum(((gen_s_prime(N_prime)[1:N] - ar_acvs[1:N])/ar_acvs[1])^2)/N))
}

temp <- sapply(2^(5:10),table_607_first_row)
round(temp[1:3],3)  # 0.500 0.106 0.007
round(temp[4:6]*c(10^5,10^9,10^14),2)  # 3.55 1.81 9.21
    
all.equal(temp[6],0)  # TRUE

### Table 607, second row

table_607_second_row <- function(N_prime,N=64)
{
    temp <- gen_s_prime(2*N_prime)
    return(sqrt(sum(((gen_s_prime(N_prime)[1:N] - temp[1:N])/temp[1])^2)/N))
}

temp <- sapply(2^(5:10),table_607_second_row)
round(temp[1:3],3)  # 0.349 0.099 0.007
round(temp[4:6]*c(10^5,10^9,10^14),2)  # 3.55 1.81 3.15
    
all.equal(temp[6],0)  # TRUE

### Figure 620 ###

fig_620 <- function(ys,tag,y_lab,delta_t=1/4)
{
    N <- length(ys)
    xs <- (0:(N-1))*delta_t
    plot(xs,ys,
         xlim=c(0,256),xlab="time  (sec)",
         ylim=c(-1200,1700),ylab=y_lab,
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 620",tag,sep=""))
    axis(1,at=seq(0,256,64))
    axis(2,at=seq(-2000,2000,1000),las=2)
    axis(2,at=seq(-2000,2000,500),label=FALSE,tcl=-0.25)
    text(270,1450,tag,pos=2)
    box(bty="l")
}

N <- length(ocean_wave)
m_parzen <- 150
m_gaussian <- 23.666

acvs_dse <- acvs(ocean_wave,center=TRUE,taper=slepian_taper(N,2))$acvs
acvs_dse_tilde <- c(acvs_dse,0,rev(acvs_dse[-1]))

w_m_parzen <- sapply(0:(N-1),parzen_lag_window,m_parzen)
w_m_parzen_tilde <- c(w_m_parzen,0,rev(w_m_parzen[-1]))
w_m_gaussian <- sapply(0:(N-1),gaussian_lag_window,m_gaussian)
w_m_gaussian_tilde <- c(w_m_gaussian,0,rev(w_m_gaussian[-1]))

S_k_parzen <- Re(dft(acvs_dse_tilde*w_m_parzen_tilde))
S_k_gaussian <- Re(dft(acvs_dse_tilde*w_m_gaussian_tilde))

set.seed(42)
cal_Y_k_parzen <- sqrt(S_k_parzen/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
Y_t_parzen <- dft(cal_Y_k_parzen)
ocean_wave_mean <- mean(ocean_wave)
ts_Re_parzen <- Re(Y_t_parzen)[1:N] + ocean_wave_mean
ts_Im_parzen <- Im(Y_t_parzen)[1:N] + ocean_wave_mean

set.seed(42)
cal_Y_k_gaussian <- sqrt(S_k_gaussian/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
Y_t_gaussian <- dft(cal_Y_k_gaussian)
ts_Re_gaussian <- Re(Y_t_gaussian)[1:N] + ocean_wave_mean
ts_Im_gaussian <- Im(Y_t_gaussian)[1:N] + ocean_wave_mean

### Figure 620, plots (a), (b), (c), (d) and (e)

fig_620(ocean_wave,"(a)",expression(italic(X[t])))
fig_620(ts_Re_parzen,"(b)",expression(paste("Re(",italic(Y[t]),")")))
fig_620(ts_Re_gaussian,"(c)",expression(paste("Re(",italic(Y[t]),")")))
fig_620(ts_Im_parzen,"(d)",expression(paste("Im(",italic(Y[t]),")")))
fig_620(ts_Im_gaussian,"(e)",expression(paste("Im(",italic(Y[t]),")")))

### Figure 621 ###

fig_621 <- function(q_x,q_y,tag="(a)",x_lab="Gaussian quantiles (IID)",little_at=seq(-5,5,0.5),big_at=seq(-6,6,1))
{
    plot(q_x,q_y,
         xlab=x_lab,
         ylim=c(-1200,1700),ylab="empirical quantiles",
         typ="p",pch=".",axes=FALSE,
         main=paste("Figure 621",tag,sep=""))
    abline(lm(q_y ~ q_x),lwd=0.5,col="gray40")
    axis(1,at=big_at)
    axis(1,at=little_at, labels=FALSE, tcl=-0.25)
    axis(2,at=seq(-2000,2000,1000),las=2)
    axis(2,at=seq(-2000,2000,100),label=FALSE,tcl=-0.25)
    temp <- range(q_x)
    text(temp[1]+0.125*diff(temp),1600,tag,pos=2)
    box(bty="l")
}

ocean_wave_centered <- ocean_wave - mean(ocean_wave)
temp <- qqnorm(ocean_wave_centered,plot=FALSE)
qq_ow_x <- temp$x[order(temp$x)]

temp <- qqnorm(ocean_wave,plot=FALSE)
qq_ow_y <- temp$y[order(temp$x)]

N <- length(ocean_wave)
m_parzen <- 150
acvs_dse <- acvs(ocean_wave,center=TRUE,taper=slepian_taper(N,2))$acvs
acvs_dse_tilde <- c(acvs_dse,0,rev(acvs_dse[-1]))
w_m_parzen <- sapply(0:(N-1),parzen_lag_window,m_parzen)
w_m_parzen_tilde <- c(w_m_parzen,0,rev(w_m_parzen[-1]))
S_k_parzen <- Re(dft(acvs_dse_tilde*w_m_parzen_tilde))

### WARNING: takes a minute or so to compute m_prime_parzen

set.seed(42)
N_rep <- 50000
m_prime_parzen <- rep(0,N)
for(n in 1:N_rep)
{
    cal_Y_k <- sqrt(S_k_parzen/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N))
    Y_t <- dft(cal_Y_k)
    temp_Re <- Re(Y_t)[1:N]
    temp_Im <- Im(Y_t)[1:N]
    m_prime_parzen <- m_prime_parzen + sort(temp_Re) + sort(temp_Im)
}
m_prime_parzen <- m_prime_parzen/(2*N_rep)

### Figure 621, left-hand plot

fig_621(qq_ow_x,qq_ow_y)

### Figure 621, right-hand plot

fig_621(m_prime_parzen,qq_ow_y,"(b)","Gaussian quantiles (correlated)",little_at=seq(-1500,1500,100),big_at=seq(-1200,1200,600))

### Figures 623a and 623b ###

fig_623a <- function(lwe,upper,lower)
{
    plot(lwe$freqs,dB(lwe$sdfe),
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 623a")
    lines(lwe$freqs,dB(upper),lwd=0.5,col="gray40")
    lines(lwe$freqs,dB(lower),lwd=0.5,col="gray40")
    cc <- lwe$cc
    x_cc <- 0.16
    y_cc <- 20
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.0,80,expression(paste("Gaussian, ", italic(m==23.666))),pos=1)
    box(bty="l")
}

fig_623b <- function(lwe,upper,lower)
{
    widths_dB <- dB(upper)-dB(lower)
    plot(lwe$freqs,widths_dB,
         xlim=c(0,2.0),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(3.2,6.0),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 623b")
    abline(h=mean(widths_dB), lwd=0.5, col="gray40")
    cc <- lwe$cc
    abline(h=cc$up+cc$down, lwd=0.5, col="gray40", lty="dashed")
    half_bandwidth <- cc$width/2
    lines(c(0,half_bandwidth),rep(3.6,2),lwd=2)
    lines(c(2.0-half_bandwidth,2.0),rep(3.6,2),lwd=2)
    axis(1,at=seq(0,2.0,0.5))
    axis(1,at=seq(0,2.0,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(3,6,1),las=2)
    axis(2,at=seq(3,6,0.2),label=FALSE,tcl=-0.25)
    box(bty="l")
}

N <- length(ocean_wave)
m_gaussian <- 23.666
the_taper <- slepian_taper(N,2)
ow_delta_t <- 1/4
gau_lwe <- ts_to_lag_window_sdf_est(ocean_wave,m_gaussian,taper=the_taper,lag_window=gaussian_lag_window,center=TRUE,delta_t=ow_delta_t)


### NOTE: evaluation of the following commented-out R code creates
###       upper_CI and lower_CI, but is a bit time-consuming (about
###       7 minutes on a 2017-vintage MacBook Pro); upper_CI and lower_CI
###       can be restored by evaluating the two load forms below
###       the commented out code
###
###         S_k <- c(gau_lwe$sdfe/ow_delta_t,rev(gau_lwe$sdfe[-c(1,N+1)]/ow_delta_t))
###         N_reps <- 50000
###         N_freqs <- length(gau_lwe$freqs)
###         lots_of_sdfes <- matrix(nrow=2*N_reps,ncol=N_freqs)
###         ow_mean <- mean(ocean_wave)
###         set.seed(42)
###         for(n in 1:N_reps)
###         {
###             cal_Y_k <- sqrt(S_k/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
###             Y_t <- dft(cal_Y_k)
###             ts_Re <- Re(Y_t[1:N]) + ow_mean
###             ts_Im <- Im(Y_t[1:N]) + ow_mean
###             lots_of_sdfes[2*n-1,] <- ts_to_lag_window_sdf_est(ts_Re,m=m_gaussian,taper=the_taper,lag_window=gaussian_lag_window,delta_t=ow_delta_t,center=TRUE)$sdfe
###             lots_of_sdfes[2*n,] <- ts_to_lag_window_sdf_est(ts_Im,m=m_gaussian,taper=the_taper,lag_window=gaussian_lag_window,delta_t=ow_delta_t,center=TRUE)$sdfe
###         }
###         
###         upper_CI <- rep(NA,N_freqs)
###         lower_CI <- rep(NA,N_freqs)
###         upper_index <- round(2*N_reps*0.975)
###         lower_index <- round(2*N_reps*0.025)
###         for(k in 1:N_freqs)
###         {
###             temp <- sort(lots_of_sdfes[,k])
###             upper_CI[k] <- temp[upper_index]
###             lower_CI[k] <- temp[lower_index]
###         }

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/upper_CI.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/lower_CI.Rdata"))

### Figure 623a

fig_623a(gau_lwe,upper_CI,lower_CI)

### Figure 623b

fig_623b(gau_lwe,upper_CI,lower_CI)

### Figures 625a, 625b and 626 ###

fig_625a <- function(ar_est,upper_boot,lower_boot)
{
    plot(ar_est$freqs,dB(ar_est$sdfe),
         xlim=c(0,2),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-40,80),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 625a")
    lines(ar_est$freqs,dB(upper_boot),lwd=0.5,col="gray40")
    lines(ar_est$freqs,dB(lower_boot),lwd=0.5,col="gray40")
    axis(1,at=seq(0,2,0.5))
    axis(1,at=seq(0,2,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-40,80,20),las=2)
    axis(2,at=seq(-40,80,10),label=FALSE,tcl=-0.25)
    text(1.0,80,"Burg, AR(25)",pos=1)
    box(bty="l")
}

fig_625b <- function(ar_est,upper_boot,lower_boot)
{
    widths_dB_boot <- dB(upper_boot)-dB(lower_boot)
    plot(ar_est$freqs,widths_dB_boot,
         xlim=c(0,2),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(0,70),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 625b")
    widths_dB_theo <- dB(ar_est$upper)-dB(ar_est$lower)
    lines(ar_est$freqs,widths_dB_theo,lwd=0.5, col="gray40")
    axis(1,at=seq(0,2,0.5))
    axis(1,at=seq(0,2,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,70,10),las=2)
    axis(2,at=seq(0,70,5),label=FALSE,tcl=-0.25)
    box(bty="l")
}


fig_626 <- function(ar_est,upper_boot,lower_boot,upper_Z,lower_Z)
{
    widths_dB_boot <- dB(upper_boot)-dB(lower_boot)
    plot(ar_est$freqs,widths_dB_boot,
         xlim=c(0,2),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(2,9.5),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 626")
    widths_dB_Z <- dB(upper_Z)-dB(lower_Z)
    lines(ar_est$freqs,widths_dB_Z,lwd=0.5, col="gray40")
    axis(1,at=seq(0,2,0.5))
    axis(1,at=seq(0,2,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(2,10,1),las=2)
    axis(2,at=seq(2,10,0.5),label=FALSE,tcl=-0.25)
    box(bty="l")
}

N <- length(ocean_wave)
ow_delta_t <- 1/4
burg_goodies_25 <- burg_algorithm(ocean_wave,25,center=TRUE)
burg_ar_est_with_CIs <- gen_cis_for_ar_sdf(N,burg_goodies_25$coeffs,burg_goodies_25$innov_var,delta_t=ow_delta_t)

### NOTE: evaluation of the following commented-out R code creates
###       upper_CI_ar_boot, lower_CI_ar_boot, upper_CI_ar_Z and lower_CI_ar_Z,
###       but is a bit time-consuming (about 20 minutes on a 2017-vintage MacBook
###       Pro); upper_CI_ar_boot etc_ can be restored by evaluating
###       the four load forms below the commented out code
###
###        pe_normalized <- c(c(ocean_wave[1]-mean(ocean_wave),sapply(1:24,function(k) burg_goodies_25$forward_pe[[k]][1]))*sqrt(burg_goodies_25$innov_var/burg_goodies_25$pev[1:25]),burg_goodies_25$forward_pe[[25]])
###        
###        N_reps_ar <- 100000
###        N_freqs_ar <- length(burg_ar_est_with_CIs$freqs)
###        lots_of_ar_sdfes <- matrix(nrow=N_reps_ar,ncol=N_freqs_ar)
###        LD_stuff_burg_25_pev_unity <- step_down_LD_recursions(burg_goodies_25$coeffs,var=1, proc=FALSE)
###        set.seed(42)
###        for(n in 1:N_reps_ar)
###        {
###            a_ts <- sim_ar_process(N,LD_stuff=LD_stuff_burg_25_pev_unity,innovations=sample(pe_normalized,N,replace=TRUE))
###            burg_goodies <- burg_algorithm(a_ts,25,center=TRUE)
###            lots_of_ar_sdfes[n,] <- ar_coeffs_to_sdf(burg_goodies$coeffs,burg_goodies$innov_var,N_pad=N,delta_t=ow_delta_t)$sdfe[-c(1,N/2+1)]
###        }
###        
###        upper_CI_ar_boot <- rep(NA,N_freqs_ar)
###        lower_CI_ar_boot <- rep(NA,N_freqs_ar)
###        upper_index <- round(N_reps_ar*0.975)
###        lower_index <- round(N_reps_ar*0.025)
###        for(k in 1:N_freqs_ar)
###        {
###            temp <- sort(lots_of_ar_sdfes[,k])
###            upper_CI_ar_boot[k] <- temp[upper_index]
###            lower_CI_ar_boot[k] <- temp[lower_index]
###        }
###        
###        LD_stuff_burg_25 <- step_down_LD_recursions(burg_goodies_25$coeffs,var=burg_goodies_25$innov_var,proc=FALSE)
###        set.seed(42)
###        for(n in 1:N_reps_ar)
###        {
###            a_ts <- sim_ar_process(N,LD_stuff=LD_stuff_burg_25)
###            burg_goodies <- burg_algorithm(a_ts,25,center=TRUE)
###            lots_of_ar_sdfes[n,] <- ar_coeffs_to_sdf(burg_goodies$coeffs,burg_goodies$innov_var,N_pad=N,delta_t=ow_delta_t)$sdfe[-c(1,N/2+1)]
###        }
###        
###        upper_CI_ar_Z <- rep(NA,N_freqs_ar)
###        lower_CI_ar_Z <- rep(NA,N_freqs_ar)
###        for(k in 1:N_freqs_ar)
###        {
###            temp <- sort(lots_of_ar_sdfes[,k])
###            upper_CI_ar_Z[k] <- temp[upper_index]
###            lower_CI_ar_Z[k] <- temp[lower_index]
###        }

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/upper_CI_ar_boot.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/lower_CI_ar_boot.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/upper_CI_ar_Z.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/lower_CI_ar_Z.Rdata"))

### Figure 625a

fig_625a(burg_ar_est_with_CIs,upper_CI_ar_boot,lower_CI_ar_boot)

### Figure 625b

fig_625b(burg_ar_est_with_CIs,upper_CI_ar_boot,lower_CI_ar_boot)

### Figure 626

fig_626(burg_ar_est_with_CIs,upper_CI_ar_boot,lower_CI_ar_boot,upper_CI_ar_Z,lower_CI_ar_Z)

### Figure 627 ###

fig_627 <- function(ys,tag,y_lab)
{
    N <- length(ys)
    xs <- 0:(N-1)
    plot(xs,ys,
         xlim=c(0,N),xlab="time  (sec)",
         ylim=c(2,12),ylab=y_lab,
         typ="l",lwd=0.25,axes=FALSE,
         main=paste("Figure 627",tag,sep=""))
    axis(1,at=seq(0,2048,512))
    axis(2,at=seq(2,12,5),las=2)
    axis(2,at=seq(2,12,1),label=FALSE,tcl=-0.25)
    text(N,11,tag,pos=2)
    box(bty="l")
}

N <- length(ship_altitude)
m_gaussian <- 80

acvs_dse <- acvs(ship_altitude,center=TRUE,taper=cosine_taper(N))$acvs
acvs_dse_tilde <- c(acvs_dse,0,rev(acvs_dse[-1]))
S_k_dse <- Re(dft(acvs_dse_tilde))
w_m_tau <- sapply(0:(N-1),gaussian_lag_window,m_gaussian)
w_m_tau_tilde <- c(w_m_tau,0,rev(w_m_tau[-1]))
S_k <- Re(dft(acvs_dse_tilde*w_m_tau_tilde))

set.seed(7)
cal_Y_k <- sqrt(S_k/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
Y_t <- dft(cal_Y_k)
ship_altitude_mean <- mean(ship_altitude)
ts_Re <- Re(Y_t) + ship_altitude_mean
ts_Im <- Im(Y_t) + ship_altitude_mean

set.seed(7)
cal_Y_k_dse <- sqrt(S_k_dse/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
Y_t_dse <- dft(cal_Y_k_dse)
ts_dse_Re <- Re(Y_t_dse) + ship_altitude_mean
ts_dse_Im <- Im(Y_t_dse) + ship_altitude_mean

### Figure 627, plots (a), (b), (c), (d) and (e)

fig_627(ship_altitude,"(a)",expression(italic(X[t])))
fig_627(ts_Re[1:N],"(b)",expression(paste("Re(",italic(Y[t]),")")))
fig_627(ts_dse_Re[1:N],"(c)",expression(paste("Re(",italic(Y[t]),")")))
fig_627(ts_Im[1:N],"(d)",expression(paste("Im(",italic(Y[t]),")")))
fig_627(ts_dse_Im[1:N],"(e)",expression(paste("Im(",italic(Y[t]),")")))

### Figure 629 ###

fig_629 <- function(list_with_CIs)
{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    ## par(mar=c(5.1,4.1,4.1,2.1))  # default
    par(mar=c(5.1,5.1,4.1,2.1))  # default
    plot_ci <- function(stuff,x,sym,cex)
    {
        eps <- 0.00014
        points(x,stuff[[2]],pch=sym,cex=cex)
        if(stuff[[2]] > stuff[[3]]) lines(rep(x,2),c(stuff[[1]],stuff[[3]]))
        else
        {
            lines(rep(x,2),c(stuff[[2]]-eps,stuff[[1]]))
            lines(rep(x,2),c(stuff[[2]]+eps,stuff[[3]]))
        }
        oops <- 0.5
        lines(c(x-oops,x+oops),rep(stuff[[1]],2))
        lines(c(x-oops,x+oops),rep(stuff[[3]],2))
    }
    plot(1:5,1:5,
         xlim=c(2.5,44.5),xaxs="i",xlab=" ",
         ylim=c(0.0195,0.0240),yaxs="i",ylab=" ",
         typ="n",axes=FALSE,
         main="Figure 629")
    mtext("innovation variance",2,3.75)
    xs <- c(seq(5,25,5),seq(32,42,5))
    syms <- c(rep(1,5),rep(15,3))
    cexs <- c(rep(1,5),rep(0.83,3))
    ##cexs <- rep(1,8)
    for(n in 1:length(list_with_CIs)) plot_ci(list_with_CIs[[n]],xs[n],syms[n],cexs[n]) 
    axis(1,at=seq(5,25,5))
    axis(1,at=seq(32,42,5), label=c("L","G","N"))
    axis(1,at=15, label="multitaper",line=2,tcl=0.0001)
    axis(1,at=37, label="Burg",line=2,tcl=0.0001)
    axis(2,at=seq(0.020,0.024,0.001),las=2)
    axis(2,at=seq(0.020,0.024,0.0005),labels=FALSE,tcl=-0.25)
    box(bty="l")
}

N_ff <- length(ac_fractional_frequencies)

### NOTE: evaluation of the following commented-out R code creates
###       ff_5, ff_10, ff_15, ff_20, ff_25 and ff_burg_42
###       but is QUITE time-consuming (calculation times on
###       a 2017-vintage MacBook Pro for creating these six variables
###       individually were, respectively, 30 minutes, 1 hour,
###       1.5 hours, 2 hours, 2.5 hours and 1 hour); ff_5 etc.
###       can be restored by evaluating the six load forms below
###       the commented out code.
###
###        do_it_mt_fig_629 <- function(ts,the_tapers,center=TRUE,N_reps=50000)
###        {
###            if(center) ts <- ts - mean(ts)
###            N <- length(ts)
###            ts_mt_innov_var <- mt_innov_var(ts,the_tapers,center=FALSE)
###            sdfe_mt <- mt_sdf_est(ts,the_tapers=the_tapers,center=FALSE,pad=2)$sdfe
###            S_k <- c(sdfe_mt,rev(sdfe_mt[-c(1,N+1)]))
###            lots_of_mt_innov_vars <- rep(NA,2*N_reps)
###            for(n in 1:N_reps)
###            {
###                cal_Y_k <- sqrt(S_k/(2*N)) * complex(real=rnorm(2*N),imag=rnorm(2*N)) 
###                Y_t <- dft(cal_Y_k)
###                ts_Re <- Re(Y_t[1:N])
###                ts_Im <- Im(Y_t[1:N])
###                lots_of_mt_innov_vars[2*n-1] <- mt_innov_var(ts_Re,the_tapers,center=TRUE)
###                lots_of_mt_innov_vars[2*n]   <- mt_innov_var(ts_Im,the_tapers,center=TRUE)
###            }
###            upper_index <- round(2*N_reps*0.975)
###            lower_index <- round(2*N_reps*0.025)
###            temp <- sort(lots_of_mt_innov_vars)
###            upper_CI <- temp[upper_index]
###            lower_CI <- temp[lower_index]
###            return(list(lower_CI=lower_CI,
###                        mt_innov_var=ts_mt_innov_var,
###                        upper_CI=upper_CI,
###                        S_k=S_k,
###                        lots_of_mt_innov_vars=lots_of_mt_innov_vars))
###        }
###        
###        do_it_burg_fig_629 <- function(ts,p=42,center=TRUE,N_reps=100000)
###        {
###            if(center) ts <- ts - mean(ts)
###            N <- length(ts)
###            burg_p <- burg_algorithm(ts,p=p)
###            ## unnormalized observed prediction errors
###            uope <- c(ts[1],sapply(1:(p-1),function(k) burg_p$forward_pe[[k]][1]),burg_p$forward_pe[[p]])
###            ## normalized observed prediction errors
###            nope <- uope/sqrt(c(burg_p$pev,rep(burg_p$innov_var,N-(p+1))))
###            LD_stuff_burg_p <- step_down_LD_recursions(burg_p$coeffs,burg_p$innov_var, proc=FALSE)
###            lots_of_ar_innov_vars_boot <- rep(NA,N_reps)
###            lots_of_ar_innov_vars_Z <- rep(NA,N_reps)
###            for(n in 1:N_reps)
###            {
###                ts_boot <- sim_ar_process(N,LD_stuff=LD_stuff_burg_p,innovations=sample(nope,N,replace=TRUE))
###                lots_of_ar_innov_vars_boot[n] <- burg_algorithm(ts_boot,center=center,p=p)$innov_var
###                ts_Z <- sim_ar_process(N,LD_stuff=LD_stuff_burg_p)
###                lots_of_ar_innov_vars_Z[n] <- burg_algorithm(ts_Z,center=center,p=p)$innov_var
###            }
###            upper_index <- round(N_reps*0.975)
###            lower_index <- round(N_reps*0.025)
###            temp <- sort(lots_of_ar_innov_vars_boot)
###            upper_CI_boot <- temp[upper_index]
###            lower_CI_boot <- temp[lower_index]
###            temp <- sort(lots_of_ar_innov_vars_Z)
###            upper_CI_Z <- temp[upper_index]
###            lower_CI_Z <- temp[lower_index]
###            return(list(ar_innov_var=burg_p$innov_var,
###                        lower_CI_boot=lower_CI_boot,
###                        upper_CI_boot=upper_CI_boot,
###                        lower_CI_Z=lower_CI_Z,
###                        upper_CI_Z=upper_CI_Z,
###                        lots_of_ar_innov_vars_boot=lots_of_ar_innov_vars_boot,
###                        lots_of_ar_innov_vars_Z=lots_of_ar_innov_vars_Z
###                        ))
###        }
###        
###        tapers_sine_K25 <- t(as.matrix(taper("sine",N_ff,n.taper=25)))
###        
###        ff_5  <- do_it_mt_fig_629(ac_fractional_frequencies,tapers_sine_K25[,1:5])
###        ff_10 <- do_it_mt_fig_629(ac_fractional_frequencies,tapers_sine_K25[,1:10])
###        ff_15 <- do_it_mt_fig_629(ac_fractional_frequencies,tapers_sine_K25[,1:15])
###        ff_20 <- do_it_mt_fig_629(ac_fractional_frequencies,tapers_sine_K25[,1:20])
###        ff_25 <- do_it_mt_fig_629(ac_fractional_frequencies,tapers_sine_K25)
###        ff_burg_42 <- do_it_burg_fig_629(ac_fractional_frequencies,42,center=TRUE)

load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ff_5.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ff_10.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ff_15.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ff_20.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ff_25.Rdata"))
load(url("http://faculty.washington.edu/dbp/sauts/Rdata/ff_burg_42.Rdata"))

### Figure 629

fig_629(list(ff_5,ff_10,ff_15,ff_20,ff_25,
             list(ff_burg_42$ar_innov_var-2*ff_burg_42$ar_innov_var*sqrt(2/N_ff),
                  ff_burg_42$ar_innov_var,
                  ff_burg_42$ar_innov_var+2*ff_burg_42$ar_innov_var*sqrt(2/N_ff)),
             list(ff_burg_42$lower_CI_Z,
                  ff_burg_42$ar_innov_var,
                  ff_burg_42$upper_CI_Z),
             list(ff_burg_42$lower_CI_boot,
                  ff_burg_42$ar_innov_var,
                  ff_burg_42$upper_CI_boot)))


### Figure 635 ###

fig_635 <- function(tag)
{
    rhotau <- seq(-1,1,0.002)
    ys <- if(tag=="(a)") (6/pi)*asin(rhotau/2)
          else
          {
              dsq <- c(0.8157660, 0.1773910,0.0066847, 0.0001343, 0.0000169, 0.0000073)
              sapply(1:length(rhotau),function(n) sum(dsq*rhotau[n]^(1:6)))
          }
    plot(rhotau,ys,
         xlim=c(-1,1),xaxs="i",xlab=expression(paste(rho[X],"  (correlation at input)")),
         ylim=c(-1,1),yaxs="i",ylab=expression(paste(rho[Y],"  (correlation at output)")),
         typ="l",axes=FALSE,
         main=paste("Figure 635",tag,sep=""))
    abline(a=0,b=1,lty="dotted")   
    axis(1,at=seq(-1,1,1))
    axis(1,at=seq(-0.5,0.5,1))
    axis(1,at=seq(-1,1,0.1),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(-1,1,1),las=2)
    axis(2,at=seq(-0.5,0.5,1),las=2)
    axis(2,at=seq(-1,1,0.1),labels=FALSE,tcl=-0.25)
    text(-0.675,0.85,tag,pos=2)
    box(bty="l")
}

### Figure 635, left-hand plot

fig_635("(a)")

### Figure 635, right-hand plot

fig_635("(b)")
