### R CODE FOR REPRODUCING CONTENT OF FIGURES AND TABLES IN CHAPTER 10 ...

ocean_wave <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean-wave.txt")
ocean_noise <- scan("http://faculty.washington.edu/dbp/sauts/Data/ocean-noise-128.txt")
spring_tides <- read.table("http://faculty.washington.edu/dbp/sauts/Data/spring-tides.txt")
neap_tides <- read.table("http://faculty.washington.edu/dbp/sauts/Data/neap-tides.txt")
M_tidal_frequencies <- scan("http://faculty.washington.edu/dbp/sauts/Data/M-tidal-frequencies.txt")
non_M_tidal_frequencies <- scan("http://faculty.washington.edu/dbp/sauts/Data/non-M-tidal-frequencies.txt")
L_60_sdf <- read.table("http://faculty.washington.edu/dbp/sauts/Data/Portsmouth-sdf-L-60.txt")
L_102_sdf <- read.table("http://faculty.washington.edu/dbp/sauts/Data/Portsmouth-sdf-L-102.txt")
Willamette_River_all <- read.table("http://faculty.washington.edu/dbp/sauts/Data/Willamette-River-395.txt")

### functions used to compute content of figures in Chapter 10 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/acvs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/AICC_given_ar_coeffs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_coeffs_to_sdf_single_freqs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/ar_least_squares.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/burg_algorithm.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_H_bar.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/B_U.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/compute_slepian_eigenvalue.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/create_tapered_series.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/cumulative_pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/daniell_design_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/direct_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dirichlet_kernel.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/discrepancy_measures.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_lwe.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/do_crisscross_mt.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/F_2_nu_alpha.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/F_2_nu_percentage_point.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/fejer_kernel.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/find_ar_3_dB_width.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/find_max_peaks_ar.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/find_max_peaks_dse.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/find_pgram_3_dB_width.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/fisher_g_statistic.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/get_dse_freqs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/hanning_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_even.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/is_odd.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/J_multiple_freq.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/lag_windows.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/mt_sdf_est.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/pgram_single_freqs.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/rectangular_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/reshaped_mt_sdfe.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/siegel_T_statistic.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/spec_window.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/step_down_LD_recursions.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/svd_harmonic_analysis.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/thomson_F_test.R")

### NOTE: to install the sapa library, uncomment the following three statements
###       and execute them:
###
### install.packages("devtools")
### devtools::install_github("wconstan/ifultools")
### devtools::install_github("wconstan/sapa")

library(sapa)

### following functions used for Figures 530, 531, 533 and 534 and Table 533

km_cos_1 <- function(t) {0.9 * cos(2*pi*(t+1)/7.5)}
km_cos_2 <- function(t) {0.9 * cos(2*pi*(t+1)/5.3 + pi/2)}
km_cos_3 <- function(t) {cos(2*pi*(t+1)/3)}

### following function used for Figures 533 and 534

gen_Fejer_snippet <- function(N,freq_center,y_max)
{
    freqs <- seq(0,0.1,0.001)
    kernel <- sapply(freqs,fejer_kernel_single_freq,N)
    gg <- kernel >= 0.05*kernel[1]
    freqs_gg <- freqs[gg]
    freqs_2_sided <- c(-rev(freqs_gg[-1]),freqs_gg) + freq_center
    kernel_gg <- kernel[gg]
    kernel_2_sided <- c(rev(kernel_gg[-1]),kernel_gg)*y_max/N
    return(list(freqs=freqs_2_sided,
                kernel=kernel_2_sided))
}

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 513 ###

fig_513 <- function(tides,main=main)
{
    x_lim <- round(range(tides$V1))
    plot(tides$V1,tides$V2,
         xlim=x_lim,xaxs="i",xlab="hours  (from beginning of 1977)",
         ylim=c(-2.5,2.5),ylab="height  (meters)",
         typ="l",axes=FALSE,
         main=main)
    axis(1,at=seq(x_lim[1],x_lim[2],10))
    axis(2,at=c(-2.5,2.5),las=2)
    axis(2,at=0,las=2)
    box(bty="l")
}

### Figure 513, top plot

fig_513(spring_tides,"Figure 513, Top Plot")

### Figure 513, bottom plot

fig_513(neap_tides,"Figure 513, Bottom Plot")

### Figure 521 ###

fig_521 <- function(M_freqs,non_M_freqs,sdfe_60,sdfe_102,x_lim=c(0,0.5))
{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    layout(rbind(1,2), heights=c(1,3))
    ## top plot ...
    par(mar=c(0.1,4.5,2.0,1.0))
    plot(abs(M_freqs),sign(M_freqs)*0.8,
         xlim=x_lim,xaxs="i",xlab=" ",
         ylim=c(-0.8,0.8),xaxs="i",ylab=" ",
         typ="h",axes=FALSE,
         main="Figure 521")
    lines(abs(non_M_freqs),sign(non_M_freqs)*0.4,typ="h")
    abline(h=0)
    ## bottom plot ...
    par(mar=c(5.0,4.5,1.0,1.0))
    plot(sdfe_60$V1,sdfe_60$V2,
         xlim=x_lim,xaxs="i",xlab=expression(paste(italic(f),"  (cycles/minute)")),
         ylim=c(-40,0),yaxs="i",ylab="dB",
         typ="n",axes=FALSE)
    poly_x <- c(-0.1,sdfe_102$V1,0.6)
    poly_y <- c(-50,sdfe_102$V2,-50)
    polygon(poly_x,poly_y,col="grey90",lwd=0.01)
    lines(sdfe_60$V1,sdfe_60$V2,lwd=1)
    lines(sdfe_102$V1, sdfe_102$V2, lwd=0.25,col="gray40")
    cc <- do_crisscross_lwe(sapply(0:(365*24-1),parzen_lag_window,400))
    x_left <- 0.47
    x_mid  <- x_left + cc$width/2
    y_cc   <- -10.0
    lines(c(x_left,x_left+cc$width), c(y_cc,y_cc), lwd=0.5)
    lines(c(x_mid,x_mid), c(y_cc-cc$down,y_cc+cc$up), lwd=0.5)
    axis(1,at=seq(x_lim[1],x_lim[2],length.out=6))
    axis(2,at=seq(-40,0,10),las=2)
}

### Figure 521

fig_521(M_tidal_frequencies,non_M_tidal_frequencies,L_60_sdf,L_102_sdf)

### Figure 530 ###

fig_530 <- function()
{
    xs_lines <- seq(-1,17,0.01)
    plot(xs_lines,km_cos_1(xs_lines),
         xlim=c(-0.5,15.5),xaxs="i",xlab=expression(italic(t)),
         ylim=c(-3,3),ylab=expression(italic(y[t])),
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 530")
    lines(xs_lines,km_cos_2(xs_lines),lty="dashed")
    lines(xs_lines,km_cos_3(xs_lines),lwd=1.5)
    xs_points <- 0:15
    km_ts <- km_cos_1(xs_points) + km_cos_2(xs_points) + km_cos_3(xs_points)
    points(xs_points,km_ts,pch=16,cex=0.75)
    axis(1,at=seq(0,15,5))
    axis(1,at=seq(0,15,1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-3,3,3),las=2)
    axis(2,at=seq(-3,3,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 530

fig_530()

### Figure 531 ###

fig_531 <- function(tag,main=main)
{
    xs_points <- 0:15
    km_ts <- km_cos_1(xs_points) + km_cos_2(xs_points) + km_cos_3(xs_points)
    pad_factor <- if(tag=="(a)") 1 else {if(tag=="(b)") 2 else 64}
    the_pgram <- pgram(km_ts,center=FALSE,pad_factor=pad_factor)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xlab=expression(italic(f)),
         ylim=c(-40,20),ylab="periodogram  (dB)",
         typ=if(pad_factor > 2) "l" else "b",pch=16,axes=FALSE,
         main=main)
    abline(v=c(1/7.5,1/5.3,1/3),lwd=1,lty="dashed")
    axis(1,at=seq(0,0.5,0.1))
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=12,tag,pos=2)
    box(bty="l")
}

### Figure 531, top, middle and bottom plots

fig_531("(a)","Figure 531(a)")
fig_531("(b)","Figure 531(b)")
fig_531("(c)","Figure 531(c)")

### Figures 533 and 534 ###

fig_533 <- function(ts,tag,main,snippet_p=FALSE)
{
    the_pgram <- pgram(ts,center=FALSE,pad_factor=64)
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xlab=expression(italic(f)),
         ylim=c(-40,20),ylab="periodogram  (dB)",
         typ="l",axes=FALSE,
         main=main)
    abline(v=c(1/7.5,1/5.3,1/3),lwd=1,lty="dashed")
    if(snippet_p)
    {
        snippet  <- gen_Fejer_snippet(length(ts),0.4,convert_from_dB(-25))
        lines(snippet$freqs, dB(snippet$kernel),lwd=2)
    }
    axis(1,at=seq(0,0.5,0.1))
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    text(x=0.5,y=12,tag,pos=2)
    box(bty="l")
}

km_ts_128 <- km_cos_1(0:127) + km_cos_2(0:127) + km_cos_3(0:127)
km_ts_64 <- km_ts_128[1:64]
km_ts_32 <- km_ts_128[1:32]
km_ts <- km_ts_128[1:16]
set.seed(2)
noisy_km_ts_1 <- km_ts + rnorm(16,sd=0.1)
noisy_km_ts_2 <- km_ts + rnorm(16,sd=sqrt(0.1))
rn_128 <- rnorm(128)
noisy_km_ts_3 <- km_ts     + rn_128[1:16]
noisy_km_ts_4 <- km_ts_32  + rn_128[1:32]
noisy_km_ts_5 <- km_ts_64  + rn_128[1:64]
noisy_km_ts_6 <- km_ts_128 + rn_128

### Figure 533, top, middle and bottom plots

fig_533(noisy_km_ts_1,"(a)","Figure 533(a)")
fig_533(noisy_km_ts_2,"(b)","Figure 533(b)")
fig_533(noisy_km_ts_3,"(c)","Figure 533(c)",snippet_p=TRUE)

### Figure 534, top, middle and bottom plots

fig_533(noisy_km_ts_4,"(a)","Figure 534(a)",snippet_p=TRUE)
fig_533(noisy_km_ts_5,"(b)","Figure 534(b)",snippet_p=TRUE)
fig_533(noisy_km_ts_6,"(c)","Figure 534(c)",snippet_p=TRUE)

### Table 533 ###

find_all_peaks <- function(ts)
{
    pgram_obj <- pgram(ts,center=FALSE,pad_factor=64)
    the_pgram <- pgram_obj$sdfe
    pg_extend <- c(the_pgram[2],the_pgram,the_pgram[length(the_pgram)-1])
    the_freqs <- pgram_obj$freqs
    freq_extend <- c(NA,the_freqs,NA)
    N <- length(pg_extend)
    temp_1 <- pg_extend[-c(1,N)] > pg_extend[-c(N-1,N)]
    temp_2 <- pg_extend[-c(1,N)] > pg_extend[-c(1,2)]
    peaks <- pg_extend[-c(1,N)] > pg_extend[-c(N-1,N)] & pg_extend[-c(1,N)] > pg_extend[-c(1,2)]
    return(list(freqs=the_freqs[peaks],sdfe=the_pgram[peaks]))
}

find_3_biggest_peaks <- function(pgram_obj)
{
    temp <- find_all_peaks(pgram_obj)
    gg <- sort(rev(order(temp$sdfe))[1:3])
    return(temp$freqs[gg])
}

### Table 533, column 2

freqs_3 <- c(1/7.5,1/5.3,1/3)
round(freqs_3,3)  # 0.133 0.189 0.333

### Table 533, columns 3, 4, 5 and 6

round(find_3_biggest_peaks(km_ts),3)                  # 0.117 0.181 0.338
round(find_3_biggest_peaks(noisy_km_ts_1),3)          # 0.117 0.181 0.337
round(find_3_biggest_peaks(noisy_km_ts_2),3)          # 0.117 0.178 0.338
temp <- round(find_3_biggest_peaks(noisy_km_ts_3),3)
c(temp[2],"-",temp[3])                                # "0.147" "-" "0.355"

### Figure 536 ###

fig_536 <- function(ts,tag,main,taper=default_taper(length(ts)),v_lines=fs,sw_peak_dB=10)
{
    the_dse <- direct_sdf_est(ts,taper,center=FALSE,pad_factor=4)
    plot(the_dse$freqs,dB(the_dse$sdfe),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-126,20),yaxs="i",ylab="dB",
         typ="l",col="gray40",axes=FALSE,
         main=main)
    abline(v=v_lines,lty="dashed")
    abline(h=-100,lty="dashed")
    the_sw <- spec_window_central_lobe(spec_window(taper,pad_factor=16,fix=TRUE,first=TRUE,last=TRUE),freq=0.05,sw_mult=convert_from_dB(sw_peak_dB),dB_cutoff=-24)
    dB_adjust <- sw_peak_dB - max(dB(the_sw$sw))
    lines(the_sw$freqs,dB(the_sw$sw)+ dB_adjust,lwd=2)
    axis(1,at=seq(0,0.5,0.1))
    axis(2,at=seq(-120,20,20),las=2)
    axis(2,at=seq(-120,20,10),label=FALSE,tcl=-0.25)
    text(x=0.481,y=-3,tag,pos=2)
    box(bty="l")
}

N <- 256
set.seed(42)
phis <- runif(3,-pi,pi)
eps <- rnorm(N,sd=sqrt(10^{-10}))
Ds <- c(0.0316, 1, 0.0001)
fs <- c(0.2943, 0.3333, 0.3971)/2
times <- 0:(N-1)
the_ts <- Ds[1]*cos(2*pi*fs[1]*times+phis[1]) + Ds[2]*cos(2*pi*fs[2]*times+phis[2]) + Ds[3]*cos(2*pi*fs[3]*times+phis[3]) + eps

### Figure 536, top, middle and bottom plots

fig_536(the_ts,"(a)","Figure 536(a)")
fig_536(the_ts,"(b)","Figure 536(b)",slepian_taper(N,2))
fig_536(the_ts,"(c)","Figure 536(c)",slepian_taper(N,4))

### Figure 538a ###

fig_538a <- function(N=64,f_1=0.0725,sig2_eps=10^{-4},freqs=seq(0.0565002,0.0884998,0.0000002))
{
    freqs_plus_f_1  <- freqs + f_1
    freqs_minus_f_1 <- freqs - f_1
    fej_kernel_1 <- sapply(freqs_plus_f_1,fejer_kernel_single_freq,N)
    fej_kernel_2 <- sapply(freqs_minus_f_1,fejer_kernel_single_freq,N)
    dir_kernel_1 <- sapply(freqs_plus_f_1,dirichlet_kernel,N)
    dir_kernel_2 <- sapply(freqs_minus_f_1,dirichlet_kernel,N)
    curve_1 <- sig2_eps + 0.25*(fej_kernel_1 + fej_kernel_2)
    phi <- 5*pi/12
    curve_2 <- sig2_eps + 0.25*(fej_kernel_1 + fej_kernel_2 + 2*N*dir_kernel_1*dir_kernel_2*cos(2*pi*(N-1)*f_1+2*phi))
    plot(freqs,curve_1,
         xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,18),yaxs="i",ylab="expected value of periodogram",
         typ="l",axes=FALSE,
         main="Figure 538a")
    lines(freqs,curve_2,col="gray40",lwd=0.5)
    abline(v=freqs[which.max(curve_1)],lwd=0.5,col="gray40",lty="dashed")
    abline(v=freqs[which.max(curve_2)],lwd=0.5,col="gray40")
    axis(1,at=seq(0.0575,0.0875,0.005))
    axis(1,at=seq(0.0575,0.0875,0.001),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,20,5),las=2)
    axis(2,at=seq(0,20,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 538a (WARNING: takes a minute or so)

fig_538a()

### Figures 538b and 557 ###

fig_538b <- function(main,N=64,N_reps=50,the_seed=42,f_1=0.0725,f_2=0.0721472,the_taper=hanning_taper(N),x_lab="periodogram peak frequency",y_lab="Hanning peak frequency",ar_p=FALSE,ar_method=ar_fbls,p=16,J=10)
{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    ## par(mar=c(5.1,4.1,4.1,2.1))  # default
    par(mar=c(5.1,6.1,4.1,2.1))  # default
    gen_ts <- function(N,f_1=0.0725,sig2_eps=10^{-4},phi=5*pi/12)
    {
        return(cos(2*pi*f_1*(0:(N-1))+phi) + rnorm(N,sd=sqrt(sig2_eps)))
    }
    set.seed(the_seed)
    ts_lots <- sapply(rep(N,N_reps),gen_ts)
    peaks_dse <- find_max_peaks_dse(ts_lots,the_taper=the_taper,center=FALSE,J=J)
    xs <- if(ar_p) find_max_peaks_ar(ts_lots,p,J=J,ar_method=ar_method,pad=512,center=FALSE) else find_max_peaks_dse(ts_lots,center=FALSE,J=J)
    plot(as.vector(xs),as.vector(peaks_dse),
         xlim=c(0.0721,0.0729),xlab=x_lab,
         ylim=c(0.0723,0.0727),ylab=" ",
         typ="p",pch=16,cex=0.25,axes=FALSE,
         main=main)
    mtext(y_lab,2,4.5)
    abline(h=f_1,col="gray40",lwd=0.5,lty="dashed")
    abline(v=f_1,col="gray40",lwd=0.5,lty="dashed")
    abline(v=f_2,col="gray40",lwd=0.5)
    axis(1,at=seq(0.0721,0.0729,0.0002))
    axis(1,at=seq(0.0721,0.0729,0.0001),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0.0721,0.0729,0.0002),las=2)
    axis(2,at=seq(0.0721,0.0729,0.0001),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 538b

fig_538b("Figure 538b")

### Figure 557

fig_538b("Figure 557",ar_p=TRUE,f_2=NULL,x_lab="FBLS AR(16) peak frequency")

### Table 542, column 1

(Ns <- 2^(4:12))  # 16 32 64 128 256 512 1024 2048 4096

### Table 542, column 2

(Ms <- (Ns-2)/2)  # 7 15 31 63 127 255 511 1023 2047

### Table 542, column 4

col_4 <- sapply(Ms,siegel_t_approximate_lam_0p6_via_M)
round(col_4[1:2],3)  # 0.253 0.146
round(col_4[3:7],4)  # 0.0861 0.0515 0.0310 0.0187 0.0113
round(col_4[8:9],5)  # 0.00686 0.00415

### Table 542, column 3 (WARNING: takes a minute or so)

col_3 <- rep(NA,9)
left_bracket <- c(0.2,0.14,0.08,col_4[4:9])
right_bracket  <- c(col_4[1:3],0.06,0.04,0.02,0.02,0.01,0.005)

for(m in 1:9) col_3[m] <- siegel_t_exact_via_M(Ms[m],c(left_bracket[m],right_bracket[m]),low=1.1)
round(col_3[1:2],3)  # 0.225 0.140
round(col_3[3:7],4)  # 0.0861 0.0523 0.0315 0.0190 0.0114
round(col_3[8:9],5)  # 0.00688 0.00416

### Table 542, column 5

col_5 <- sapply(Ms,function(M) {temp <- siegel_c_and_beta(M); qchisq(0.95,0,temp$beta,lower.tail=TRUE)*temp$c})
round(col_5[1:2],3)  # 0.220 0.138
round(col_5[3:7],4)  # 0.0852 0.0520 0.0315 0.0190 0.0114
round(col_5[8:9],5)  # 0.00688 0.00416

### Table 542, column 7

col_7 <- sapply(Ms,siegel_t_approximate_lam_0p6_via_M,alpha=0.01)
round(col_7[1:2],3)  # 0.318 0.173
round(col_7[3:7],4)  # 0.0971 0.0552 0.0316 0.0181 0.0104
round(col_7[8:9],5)  # 0.00598 0.00344

### Table 542, column 6 (WARNING: takes a minute or so)

col_6 <- rep(NA,9)
left_bracket <- c(0.2,0.16,0.08,col_7[4:7],0.005,0.003)
right_bracket  <- c(col_7[1:3],0.06,0.04,0.02,0.02,col_7[8:9])

for(m in 1:9) col_6[m] <- siegel_t_exact_via_M(Ms[m],c(left_bracket[m],right_bracket[m]),alpha=0.01,low=1.1)
round(col_6[1:2],3)  # 0.268 0.164
round(col_6[3:7],4)  # 0.0969 0.0562 0.0322 0.0184 0.0104
round(col_6[8:9],5)  # 0.00595 0.00341

### Table 542, column 8

col_8 <- sapply(Ms,function(M) {temp <- siegel_c_and_beta(M,alpha=0.01); qchisq(0.99,0,temp$beta,lower.tail=TRUE)*temp$c})
round(col_8[1:2],3)  # 0.273 0.165
round(col_8[3:7],4)  # 0.0974 0.0564 0.0323 0.0184 0.0105
round(col_8[8:9],5)  # 0.00597 0.00341

### Figure 551 ###

fig_551 <- function(the_pgram,f_hats)
{
    plot(the_pgram$freqs,dB(the_pgram$sdfe),
         xlim=c(0,0.5),xlab=expression(italic(f)),
         ylim=c(-40,20),ylab="periodogram  (dB)",
         typ="l",axes=FALSE,
         main="Figure 551")
    abline(v=f_hats,lwd=1,lty="dashed")
    axis(1,at=seq(0,0.5,0.1))
    axis(2,at=seq(-40,20,20),las=2)
    axis(2,at=seq(-40,20,10),label=FALSE,tcl=-0.25)
    box(bty="l")
}

R_3 <- function(parms,ts=noisy_km_ts_6)
{
    N <- length(ts)
    return(ts - as.vector(cbind(sapply(parms[7:9],function(f) cos(2*pi*f*(0:(N-1)))),
                                sapply(parms[7:9],function(f) sin(2*pi*f*(0:(N-1))))) %*% parms[1:6]))
}

km_ts_128 <- km_cos_1(0:127) + km_cos_2(0:127) + km_cos_3(0:127)
set.seed(2)
rn_128 <- rnorm(160)[33:160]
noisy_km_ts_6 <- km_ts_128 + rn_128
f_hats <- as.vector(find_max_peaks_dse(noisy_km_ts_6,N_peaks=3,center=FALSE))
N <- length(noisy_km_ts_6)
X_reg <-cbind(sapply(f_hats,function(f) cos(2*pi*f*(0:(N-1)))),
              sapply(f_hats,function(f) sin(2*pi*f*(0:(N-1)))))
initial_parms <- c((2*as.vector(t(X_reg) %*% noisy_km_ts_6)/N),f_hats)
temp <- nlm(function(...) sum(R_3(...)^2),initial_parms,gradtol = 1e-8,steptol = 1e-8)
R_euls <- R_3(temp$est)
pg_R <- pgram(R_euls,center=FALSE,pad_factor=8)

### Figure 551

fig_551(pg_R,temp$est[7:9])

### Figure 568a ###

fig_568a <- function(ts)
{
    plot(ts$V1,ts$V2,
         xlim=c(1950,1985),xaxs="i",xlab="time  (years)",
         ylim=c(7.5,12.5),yaxs="i",ylab="log(flow)",
         typ="l",axes=FALSE,
         main="Figure 568a")
    abline(h=mean(ts$V2),lty="dashed",col="gray40")
    abline(v=c(1952.792,1981.542),lty="dashed",col="gray40")
    axis(1,at=seq(1950,1984,17))
    axis(1,at=seq(1950,1984,1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(8,12,2),las=2)
    axis(2,at=seq(8,12,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 568a

fig_568a(Willamette_River_all)

### Figure 568b ###

fig_568b <- function(ts,delta_t=1/12)
{
    the_pgram <- pgram(ts$V2,center=TRUE,pad=1024/length(ts$V2),delta_t=delta_t)
    plot(the_pgram$freqs[-1],dB(the_pgram$sdfe[-1]),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-50,10),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main="Figure 568b")
    x_cc <- 5.75
    y_cc <- -8
    lines(c(x_cc,x_cc),y_cc+c(the_pgram$cc$up,-the_pgram$cc$down),lwd=0.5)
    lines(x_cc+c(-the_pgram$cc$width/2,the_pgram$cc$width/2),c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,6,1))
    axis(2,at=seq(-50,10,10),las=2)
    box(bty="l")
}

### Figure 568b

fig_568b(Willamette_River_all)

### Table 570 ###

wr_centered <- Willamette_River_all$V2 - mean(Willamette_River_all$V2)
N <- length(wr_centered)
delta_t <- 1/12

### Table 570, first row

f_hat <- as.vector(find_max_peaks_dse(wr_centered,center=FALSE,delta_t=1/12))
A_hat <- 2*sum(wr_centered*cos(2*pi*f_hat*(0:(N-1))*delta_t))/N
B_hat <- 2*sum(wr_centered*sin(2*pi*f_hat*(0:(N-1))*delta_t))/N
SS <- sum((wr_centered - A_hat*cos(2*pi*f_hat*(0:(N-1))*delta_t) - B_hat*sin(2*pi*f_hat*(0:(N-1))*delta_t))^2)

round(f_hat,5)                  #  1.00319
round(A_hat,4)                  # -0.3033
round(B_hat,4)                  #  0.8443
round((A_hat^2 + B_hat^2)/4,4)  #  0.2012
round(SS,5)                     # 87.40164

### Table 570, second row

X_row_2 <- cbind(cos(2*pi*f_hat*(0:(N-1))*delta_t),sin(2*pi*f_hat*(0:(N-1))*delta_t))
lm_row_2 <- lm(wr_centered ~ X_row_2 - 1)

round(f_hat,5)                    #  1.00319
round(coef(lm_row_2)[1],4)        # -0.3029
round(coef(lm_row_2)[2],4)        #  0.8447
round(sum(coef(lm_row_2)^2)/4,4)  #  0.2013
round(sum(resid(lm_row_2)^2),5)   # 87.40157

### Table 570, third row

SS_3 <- function(parms,ts=wr_centered)
{
    N <- length(ts)
    return(sum((ts - as.vector(cbind(cos(2*pi*parms[3]*(0:(N-1))*delta_t),sin(2*pi*parms[3]*(0:(N-1))*delta_t)) %*% parms[1:2]))^2))
}
    
initial_parms <- c(coef(lm_row_2),f_hat)
nlm_row_3 <- nlm(SS_3,initial_parms,gradtol = 1e-8,steptol = 1e-8)

round(nlm_row_3$est[3],5)             #  1.00326
round(nlm_row_3$est[1],4)             # -0.3086
round(nlm_row_3$est[2],4)             #  0.8427
round(sum(nlm_row_3$est[1:2]^2)/4,4)  #  0.2013
round(SS_3(nlm_row_3$est),5)          # 87.39881

### Table 570, fourth row

f_hats <- as.vector(find_max_peaks_dse(wr_centered,center=FALSE,delta_t=1/12,N_peaks=2))
A_1_hat <- 2*sum(wr_centered*cos(2*pi*f_hats[1]*(0:(N-1))*delta_t))/N
B_1_hat <- 2*sum(wr_centered*sin(2*pi*f_hats[1]*(0:(N-1))*delta_t))/N
A_2_hat <- 2*sum(wr_centered*cos(2*pi*f_hats[2]*(0:(N-1))*delta_t))/N
B_2_hat <- 2*sum(wr_centered*sin(2*pi*f_hats[2]*(0:(N-1))*delta_t))/N
SS <- sum((wr_centered - A_1_hat*cos(2*pi*f_hats[1]*(0:(N-1))*delta_t) - B_1_hat*sin(2*pi*f_hats[1]*(0:(N-1))*delta_t) - A_2_hat*cos(2*pi*f_hats[2]*(0:(N-1))*delta_t) - B_2_hat*sin(2*pi*f_hats[2]*(0:(N-1))*delta_t))^2)

round(f_hats[1],5)                  #  1.00319
round(A_1_hat,4)                    # -0.3033
round(B_1_hat,4)                    #  0.8443
round((A_1_hat^2 + B_1_hat^2)/4,4)  #  0.2012
round(f_hats[2],5)                  #  2.00339
round(A_2_hat,4)                    #  0.0365
round(B_2_hat,4)                    #  0.1937
round((A_2_hat^2 + B_2_hat^2)/4,5)  #  0.00971
round(SS,5)                         # 79.58391

### Table 570, fifth row

X_row_5 <- cbind(sapply(f_hats,function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(f_hats,function(f) sin(2*pi*f*(0:(N-1))*delta_t)))
lm_row_5 <- lm(wr_centered ~ X_row_5 - 1)

round(f_hats[1],5)                        #  1.00319
round(coef(lm_row_5)[1],4)                # -0.303
round(coef(lm_row_5)[3],4)                #  0.8451
round(sum(coef(lm_row_5)[c(1,3)]^2)/4,4)  #  0.2015
round(f_hats[2],5)                        #  2.00339
round(coef(lm_row_5)[2],4)                #  0.0364 
round(coef(lm_row_5)[4],4)                #  0.1956 
round(sum(coef(lm_row_5)[c(2,4)]^2)/4,5)  #  0.00989
round(sum(resid(lm_row_5)^2),5)           # 79.58305

### Table 570, sixth row

SS_6 <- function(parms,ts=wr_centered)
{
    N <- length(ts)
    return(sum((ts - as.vector(cbind(sapply(parms[5:6],function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(parms[5:6],function(f) sin(2*pi*f*(0:(N-1))*delta_t))) %*% parms[1:4]))^2))
}
    
initial_parms <- c(coef(lm_row_5),f_hats)
nlm_row_6 <- nlm(SS_6,initial_parms,gradtol = 1e-8,steptol = 1e-8)

round(nlm_row_6$est[5],5)                #  1.00333
round(nlm_row_6$est[1],4)                # -0.3147
round(nlm_row_6$est[3],4)                #  0.8411 
round(sum(nlm_row_6$est[c(1,3)]^2)/4,4)  #  0.2016
round(nlm_row_6$est[6],5)                #  2.00244
round(nlm_row_6$est[2],4)                #  0.0553
round(nlm_row_6$est[4],4)                #  0.1915
round(sum(nlm_row_6$est[c(2,4)]^2)/4,5)  #  0.00994 
round(SS_6(nlm_row_6$est),5)             # 79.54593

### Figures 571a, 571b, 572, 573 and 574 ###

fig_571a <- function(xs,ys,tag,main,y_lim=c(7.5,12.5),y_lab="fitted log(flow)",h_line=NULL)
{
    plot(xs,ys,
         xlim=c(1950,1985),xaxs="i",xlab="time  (years)",
         ylim=y_lim,yaxs="i",ylab=y_lab,
         typ="l",axes=FALSE,
         main=main)
    abline(h=h_line,lty="dashed",col="gray40")
    axis(1,at=seq(1950,1984,17))
    axis(1,at=seq(1950,1984,1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(ceiling(y_lim[1]),floor(y_lim[2]),2),las=2)
    axis(2,at=seq(ceiling(y_lim[1]),floor(y_lim[2]),1),label=FALSE,tcl=-0.25)
    text(1983.5,1.5+mean(y_lim),tag,pos=2)
    box(bty="l")
}

fig_571b <- function(ts,main,delta_t=1/12)
{
    the_pgram <- pgram(ts,center=FALSE,pad=1024/length(ts),delta_t=delta_t)
    plot(the_pgram$freqs[-1],dB(the_pgram$sdfe[-1]),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-50,10),yaxs="i",ylab="dB",
         typ="l",axes=FALSE,
         main=main)
    x_cc <- 5.75
    y_cc <- -8
    lines(c(x_cc,x_cc),y_cc+c(the_pgram$cc$up,-the_pgram$cc$down),lwd=0.5)
    lines(x_cc+c(-the_pgram$cc$width/2,the_pgram$cc$width/2),c(y_cc,y_cc),lwd=0.5)
    axis(1,at=seq(0,6,1))
    axis(2,at=seq(-50,10,10),las=2)
    box(bty="l")
}

fig_574 <- function(resids,v_lines,delta_t=1/12)
{
    the_c_pgram <- cumulative_pgram(resids,delta_t=delta_t)
    plot(the_c_pgram$freqs,the_c_pgram$c_pgram,
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(0,1),yaxs="i",ylab="cumulative periodogram",
         typ="l",axes=FALSE,
         main="Figure 574")
    lines(the_c_pgram$L_u[,1],the_c_pgram$L_u[,2],lwd=0.5,lty="dashed")
    lines(the_c_pgram$L_l[,1],the_c_pgram$L_l[,2],lwd=0.5,lty="dashed")
    abline(v=v_lines,lwd=0.5,lty="dashed")
    abline(v=the_c_pgram$f_D_stat,lwd=0.5)
    axis(1,at=seq(0,6,1))
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(0,1,0.1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

X_bar <- mean(Willamette_River_all$V2)
wr_centered <- Willamette_River_all$V2 - X_bar
N <- length(wr_centered)
delta_t <- 1/12
f_hats <- as.vector(find_max_peaks_dse(wr_centered,center=FALSE,delta_t=1/12,N_peaks=2))

X_ECLS_1 <- cbind(cos(2*pi*f_hats[1]*(0:(N-1))*delta_t),sin(2*pi*f_hats[1]*(0:(N-1))*delta_t))
lm_ECLS_1 <- lm(wr_centered ~ X_ECLS_1 - 1)

resids_1 <- function(parms,ts=wr_centered)
{
    N <- length(ts)
    return(ts - as.vector(cbind(cos(2*pi*parms[3]*(0:(N-1))*delta_t),sin(2*pi*parms[3]*(0:(N-1))*delta_t)) %*% parms[1:2]))
}
    
initial_parms <- c(coef(lm_ECLS_1),f_hats[1])
nlm_EULS_1 <- nlm(function(...) sum(resids_1(...)^2),initial_parms,gradtol = 1e-8,steptol = 1e-8)

X_ECLS_2 <- cbind(sapply(f_hats,function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(f_hats,function(f) sin(2*pi*f*(0:(N-1))*delta_t)))
lm_ECLS_2 <- lm(wr_centered ~ X_ECLS_2 - 1)

resids_2 <- function(parms,ts=wr_centered)
{
    N <- length(ts)
    return(ts - as.vector(cbind(sapply(parms[5:6],function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(parms[5:6],function(f) sin(2*pi*f*(0:(N-1))*delta_t))) %*% parms[1:4]))
}

initial_parms <- c(coef(lm_ECLS_2),f_hats)
nlm_EULS_2 <- nlm(function(...) sum(resids_2(...)^2),initial_parms,gradtol = 1e-8,steptol = 1e-8)

### Figure 571a, top plot

fig_571a(Willamette_River_all$V1,Willamette_River_all$V2-resids_1(nlm_EULS_1$est),"(a)","Figure 571a(a)",h_line=X_bar)

### Figure 571a, bottom plot

fig_571a(Willamette_River_all$V1,resids_1(nlm_EULS_1$est),"(b)","Figure 571a(b)",y_lim=c(-2.5,2.5))

### Figure 571b

fig_571b(resids_1(nlm_EULS_1$est),"Figure 571b")

### Figure 572, top plot

fig_571a(Willamette_River_all$V1,Willamette_River_all$V2-resids_2(nlm_EULS_2$est),"(a)","Figure 572(a)",h_line=X_bar)

### Figure 572, bottom plot

fig_571a(Willamette_River_all$V1,resids_2(nlm_EULS_2$est),"(b)","Figure 572(b)",y_lim=c(-2.5,2.5))

### Figure 573

fig_571b(resids_2(nlm_EULS_2$est),"Figure 573")

### Figure 574

fig_574(resids_2(nlm_EULS_2$est),f_hats)

### Figure 575 ###

fig_575_top <- function(ts,which_tapers="dpss",K_tapers=5,bandwidth=4,pad_factor=1024/N_ts,delta_t=1/12)
{
    N_ts <- length(ts)
    the_tapers <- t(as.matrix(taper(which_tapers,N_ts,n.taper=K_tapers,bandwidth=bandwidth)))
    F_test_results <- thomson_F_test(ts,the_tapers,pad_factor=pad_factor,delta_t=delta_t,center=TRUE)
    plot(F_test_results$freqs[-1],F_test_results$F_test[-1],
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(0,53),yaxs="i",ylab=expression(paste(italic(F),"-test")),
         typ="l",axes=FALSE,
         main="Figure 575(a)")
    abline(h=sapply(c(0.01,0.001), function(alpha) F_2_nu_percentage_point(2*K_tapers-2,alpha)),lwd=0.5,lty="longdash")
    abline(h=F_2_nu_percentage_point(2*K_tapers-2,1/N_ts),lwd=0.5,lty="dotted")
    axis(1,at=seq(0,6,1))
    axis(2,at=seq(0,50,10),las=2)
    text(5.9,46.5,"(a)",pos=2)
    box(bty="l")
}

fig_575_bottom <- function(ts,K_tapers=5,which_tapers="dpss",bandwidth=4,pad_factor=1024/N_ts,delta_t=1/12,reshape_indices=c(86,171))
{
    N_ts <- length(ts)
    the_tapers <- t(as.matrix(taper(which_tapers,N_ts,n.taper=K_tapers,bandwidth=bandwidth)))
    the_mt_sdfe <- mt_sdf_est(ts,the_tapers,center=TRUE,pad_factor=pad_factor,delta_t=delta_t)
    plot(the_mt_sdfe$freqs[-1],dB(the_mt_sdfe$sdfe[-1]),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-50,10),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,col="gray40",axes=FALSE,
         main="Figure 575(b)")
    x_cc <- 5.75
    y_cc <- -8
    cc <- the_mt_sdfe$cc
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
    ## reshaping ...
    F_test_results <- thomson_F_test(ts,the_tapers,pad_factor=pad_factor,delta_t=delta_t,center=TRUE)
    for (i in 1:length(reshape_indices))
    {
        temp <- reshaped_mt_sdfe(reshape_indices[i],F_test_results,the_mt_sdfe)
        lines(temp$freqs,dB(temp$sdfe),lwd=2)
    }
    axis(1,at=seq(0,6,1))
    axis(2,at=seq(-50,10,10),las=2)
    text(5.9,2.8,"(b)",pos=2)
    box(bty="l")
}

### Figure 575, top plot

fig_575_top(Willamette_River_all$V2)

### Figure 575, bottom plot

fig_575_bottom(Willamette_River_all$V2)

### Figure 576 ###

fig_576 <- function(ts,which_ar_method=burg_algorithm,p_max=150,center=TRUE,delta_t=1/12)
{
    blpc <-which_ar_method(ts,p=p_max,center=center)$blpc
    AICCs <- sapply(1:p_max,function(p) AICC_given_ar_coeffs(ts,blpc[[p]],center=TRUE)$AICC)
    plot(1:p_max,AICCs,
         xlim=c(0,p_max),xaxs="i",xlab=expression(italic(p)),
         ylab=expression(paste("AICC(",italic(p),")")),
         typ="l",axes=FALSE,
         main="Figure 576")
    abline(v=which.min(AICCs),col="gray40",lwd=0.5)
    axis(1,at=seq(0,p_max,25))
    axis(1,at=seq(0,p_max,5),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,1000,50),las=2)
    box(bty="l")
}

### Figure 576

fig_576(Willamette_River_all$V2)

### Figure 577 ###

fig_577_spec <- function(ts,p,main,which_ar_method=burg_algorithm,center=TRUE,delta_t=1/12,N_pad=65536)
{
    ar_results <-which_ar_method(ts,p,center=center)
    the_ar_sdf <- ar_coeffs_to_sdf(ar_results$coeffs,ar_results$innov_var,delta_t=delta_t,N_pad=N_pad)
    plot(the_ar_sdf$freqs,dB(the_ar_sdf$sdf),
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=c(-40,30),yaxs="i",ylab="dB",
         typ="l",col="gray40",axes=FALSE,
         main=main)
    axis(1,at=seq(0,6,1))
    axis(2,at=seq(-40,30,10),las=2)
    text(5.8,19,paste("AR(",p,")",sep=""),pos=2)
    box(bty="l")
}


fig_577_int_spec <- function(ts,p,which_ar_method=burg_algorithm,center=TRUE,delta_t=1/12,N_pad=65536)
{
    ar_results <-which_ar_method(ts,p,center=center)
    the_ar_sdf <- ar_coeffs_to_sdf(ar_results$coeffs,ar_results$innov_var,delta_t=delta_t,N_pad=N_pad)
    N_f <- length(the_ar_sdf$freqs)
    delta_f <- the_ar_sdf$freqs[2]
    two_sided_sdf <- c(rev(the_ar_sdf$sdfe[-1]),the_ar_sdf$sdfe)
    ## burg_27_two_sided_freqs <- c(-rev(burg_27_sdf$freqs[-1]),burg_27_sdf$freqs)
    int_spec_half <- delta_f*cumsum(two_sided_sdf)[N_f:(2*N_f-1)]
    y_lim <- range(int_spec_half)
    plot(the_ar_sdf$freqs,int_spec_half,
         xlim=c(0,6),xaxs="i",xlab=expression(paste(italic(f),"  (cycles/year)")),
         ylim=y_lim,yaxs="i",ylab="integrated spectrum",
         typ="l",col="gray40",axes=FALSE,
         main="Figure 577, Bottom Plot")
    N_ts <- length(ts)
    f_hat <- as.vector(find_max_peaks_dse(ts,center=TRUE,delta_t=1/12))
    A_hat <- 2*sum(wr_centered*cos(2*pi*f_hat*(0:(N_ts-1))*delta_t))/N_ts
    B_hat <- 2*sum(wr_centered*sin(2*pi*f_hat*(0:(N_ts-1))*delta_t))/N_ts
    jump <- (A_hat^2 + B_hat^2)/4
    f_peak <- as.vector(find_max_peaks_ar(ts,p,ar_method=which_ar_method,N_pad=N_pad,center=center,delta_t=delta_t))
    j_f_peak <- which.min(abs(the_ar_sdf$freqs-f_peak))
    k_spread <- which.min(abs(int_spec_half[j_f_peak+(1:500)]-int_spec_half[j_f_peak-(1:500)] - jump))
    hlines_for_ispec <- c(int_spec_half[j_f_peak-k_spread],int_spec_half[j_f_peak+k_spread])
    lines(c(4.5,4.5),c(hlines_for_ispec[2],hlines_for_ispec[2]-jump*0.4),lwd=0.5)
    lines(c(4.5,4.5),c(hlines_for_ispec[1],hlines_for_ispec[1]+jump*0.4),lwd=0.5)
    diddle_1 <- 0.00675
    text(4.5,hlines_for_ispec[2]-diddle_1,"<",srt=270)
    text(4.5,hlines_for_ispec[1]+diddle_1,"<",srt=90)
    diddle_2 <-  0.02
    text(4.5,hlines_for_ispec[1]+jump/2+diddle_2,round(jump,4),pos=1)
    abline(h=hlines_for_ispec,col="gray40",lwd=0.5)
    axis(1,at=seq(0,6,1))
    axis(2,at=seq(0,1,0.05),las=2)
    text(5.8,y_lim[1]+123*diff(y_lim)/146,paste("AR(",p,")",sep=""),pos=2)
    box(bty="l")
}

### Figure 577, top two plots

fig_577_spec(Willamette_River_all$V2,150,"Figure 577, Top Plot")
fig_577_spec(Willamette_River_all$V2,27,"Figure 577, Middle Plot")

### Figure 577, bottom plot

fig_577_int_spec(Willamette_River_all$V2,27)

### Table 578 ###

peaks_150 <- find_max_peaks_ar(Willamette_River_all$V2,150,N_peaks=2,ar_method=burg_algorithm,N_pad=65536,center=TRUE,delta_t=1/12)
ar_results_150 <-burg_algorithm(Willamette_River_all$V2,150,center=TRUE)
width_peak_1_150 <- find_ar_3_dB_width(peaks_150[1],ar_results_150$innov_var,ar_results_150$coeffs,delta_t=1/12)$width_3dB
width_peak_2_150 <- find_ar_3_dB_width(peaks_150[2],ar_results_150$innov_var,ar_results_150$coeffs,delta_t=1/12)$width_3dB
S_hat_peaks_150 <- dB(ar_coeffs_to_sdf_single_freqs(peaks_150,innov_var=ar_results_150$innov_var,coeffs=ar_results_150$coeffs,delta_t=1/12))

peaks_27 <- find_max_peaks_ar(Willamette_River_all$V2,27,N_peaks=2,ar_method=burg_algorithm,N_pad=65536,center=TRUE,delta_t=1/12)
ar_results_27 <-burg_algorithm(Willamette_River_all$V2,27,center=TRUE)
width_peak_1_27 <- find_ar_3_dB_width(peaks_27[1],ar_results_27$innov_var,ar_results_27$coeffs,delta_t=1/12)$width_3dB
width_peak_2_27 <- find_ar_3_dB_width(peaks_27[2],ar_results_27$innov_var,ar_results_27$coeffs,delta_t=1/12)$width_3dB
S_hat_peaks_27 <- dB(ar_coeffs_to_sdf_single_freqs(peaks_27,innov_var=ar_results_27$innov_var,coeffs=ar_results_27$coeffs,delta_t=1/12))

peaks_pgram <- find_max_peaks_dse(Willamette_River_all$V2,N_peaks=2,center=TRUE,delta_t=1/12)
width_peak_1_pgram <- find_pgram_3_dB_width(Willamette_River_all$V2,peaks_pgram[1],delta_t=1/12)$width_3dB
width_peak_2_pgram <- find_pgram_3_dB_width(Willamette_River_all$V2,peaks_pgram[2],delta_t=1/12,N_pad=65536)$width_3dB
S_hat_peaks_pgram <- dB(pgram_single_freqs(peaks_pgram,ts=Willamette_River_all$V2,delta_t=1/12)$sdfe)

### Table 578, 1st row

round(peaks_150[1],4)        #  1.0055
round(width_peak_1_150,4)    #  9e-04
round(S_hat_peaks_150[1],1)  #  21.6
round(peaks_150[2],4)        #  2.0055
round(width_peak_2_150,4)    #  0.0091
round(S_hat_peaks_150[2],1)  # -1.1

### Table 578, 2nd row

round(peaks_27[1],4)        #  1.0042
round(width_peak_1_27,4)    #  0.0126
round(S_hat_peaks_27[1],1)  # 10.3
round(peaks_27[2],4)        #  1.9757
round(width_peak_2_27,4)    #  0.1225
round(S_hat_peaks_27[2],1)  # -9.8

### Table 578, 3rd row

round(peaks_pgram[1],4)        #  1.0032
round(width_peak_1_pgram,4)    #  0.0265
round(S_hat_peaks_pgram[1],1)  #  8.2
round(peaks_pgram[2],4)        #  2.0034
round(width_peak_2_pgram,4)    #  0.0269
round(S_hat_peaks_pgram[2],1)  # -5

### Figures 579, 580a and 580b ###

fig_579 <- function(years,ts)
{
    plot(years,ts,
         xlim=c(1950,1985),xaxs="i",xlab="time  (years)",
         ylim=c(-2.5,2.5),yaxs="i",ylab="imaginary part",
         typ="l",axes=FALSE,
         main="Figure 579")
    abline(h=0,lty="dashed",col="gray40")
    axis(1,at=seq(1950,1984,17))
    axis(1,at=seq(1950,1984,1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-2,2,2),las=2)
    axis(2,at=seq(-2,2,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

fig_580a <- function(eig)
{
    plot(1:length(eig),eig,
         xlim=c(1,26.5),xlab=expression(italic(k)),
         ylab=expression(lambda[k]),
         typ="p",axes=FALSE,
         main="Figure 580a")
    axis(1,at=seq(0,30,5))
    axis(1,at=seq(0,30,1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,120,20),las=2)
    box(bty="l")
}

fig_580b <- function(the_roots,mark_largest=NULL,span=1.01)
{
    unit_circle_angles <- seq(0,2*pi,pi/1000)
    x_unit_circle <- cos(unit_circle_angles)
    y_unit_circle <- sin(unit_circle_angles)
    plot(x_unit_circle,y_unit_circle,
         xlim=span*c(-1,1),xlab=" ",
         ylim=span*c(-1,1),ylab=" ",
         typ="l",axes=FALSE,
         main="Figure 580b")
    lines(c(-1,1), c(0,0))
    lines( c(0,0),c(-1,1))
   for(an_angle in c(30,60,120,150,210,240,300,330))
    {
        lines(c(0.5*cos(2*pi*an_angle/360),cos(2*pi*an_angle/360)),
              c(0.5*sin(2*pi*an_angle/360),sin(2*pi*an_angle/360)))
    }
    points(Re(the_roots), Im(the_roots))
    if(!is.null(mark_largest))
    {
        i_mark <- head(rev(order(abs(the_roots))),n=mark_largest)
        points(Re(the_roots[i_mark]),Im(the_roots[i_mark]), pch=16)
    }
    mtext(expression(270*degree),1,-0.5)
    mtext(expression(180*degree),2,-0.5,las=2)
    mtext(expression(90*degree),3,-0.5)
    mtext(expression(0*degree),4,-0.5,las=2)
}

hilbert_irs <- function(k) if(is_odd(k)) 2/(pi*k) else 0

hanning_lag_window <- function(tau,m) if(abs(tau) > m) 0 else 1 - 0.5 + 0.5*cos(pi*tau/m)

wr_centered <- Willamette_River_all$V2-mean(Willamette_River_all$V2)
the_filter <- sapply((-25):25,hanning_lag_window,26)*sapply((-25):25,hilbert_irs)
imag_345  <- convolve(wr_centered,rev(the_filter),type="filter")
years_345 <- Willamette_River_all$V1[26:370]
real_345  <- wr_centered[26:370]
z_345 <- complex(real=real_345,imag=imag_345)
svd_results <- svd_harmonic_analysis(z_345,p=3,p_prime=27)

### Figure 579

fig_579(years_345,imag_345)

### Figure 580a

fig_580a(svd_results$eigen)
    
### Figure 580b

fig_580b(svd_results$roots,mark=2)

### Figure 581 ###

fig_581 <- function(ts,delta_t=1)
{
    N <- length(ts)
    M <- floor((N-1)/2)
    temp <- pgram(ts,delta_t=delta_t)
    the_pgram <- temp$sdfe[2:(M+1)]
    the_freqs <- temp$freqs[2:(M+1)]
    pgram_rescaled <- the_pgram/sum(the_pgram)
    plot(the_freqs,pgram_rescaled,
         xlim=c(0,0.5),xaxs="i",xlab=expression(paste(italic(f[k]),"  (Hz)")),
         ylim=c(-0.003,0.15),yaxs="i",ylab="rescaled periodogram",
         typ="p",pch=16,axes=FALSE,
         main="Figure 581")
    g_F <- g_F_exact_via_M(M)
    abline(h=c(g_F,0.6*g_F), col="gray40", lwd=0.5)
    gg <- pgram_rescaled > 0.6*g_F
    segments(the_freqs[gg],rep(0.6*g_F,sum(gg)),the_freqs[gg],pgram_rescaled[gg], col="gray40", lwd=0.5)
    axis(1,at=seq(0,0.5,0.1))
    axis(1,at=seq(0,0.5,0.01),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,0.15,0.05),las=2)
    axis(2,at=seq(0,0.15,0.01),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 581

fig_581(ocean_noise)

### Figure 582 ###

fig_582_top <- function(ts,which_tapers="dpss",K_tapers=5,bandwidth=4,pad_factor=4,delta_t=1/4,reshape_indices=c(86,171))
{
    N_ts <- length(ts)
    the_tapers <- t(as.matrix(taper(which_tapers,N_ts,n.taper=K_tapers,bandwidth=bandwidth)))
    F_test_results <- thomson_F_test(ts,the_tapers,pad_factor=pad_factor,delta_t=delta_t,center=TRUE)
    plot(F_test_results$freqs[-1],F_test_results$F_test[-1],
         xlim=c(0,2),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(0,90),yaxs="i",ylab=expression(paste(italic(F),"-test")),
         typ="l",axes=FALSE,
         main="Figure 582(a)")
    abline(h=sapply(c(0.01,0.001), function(alpha) F_2_nu_percentage_point(2*K_tapers-2,alpha)),lwd=0.5,lty="longdash")
    abline(h=F_2_nu_percentage_point(2*K_tapers-2,1/N_ts),lwd=0.5,lty="dotted")
    axis(1,at=seq(0,2,0.5))
    axis(1,at=seq(0,2,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,100,20),las=2)
    axis(2,at=seq(0,100,10),label=FALSE,tcl=-0.25)
    text(1.97,79.2,"(a)",pos=2)
    box(bty="l")
}

fig_582_bottom <- function(ts,K_tapers=5,which_tapers="dpss",bandwidth=4,pad_factor=4,delta_t=1/4,reshape_indices=1505)
{
    N_ts <- length(ts)
    the_tapers <- t(as.matrix(taper(which_tapers,N_ts,n.taper=K_tapers,bandwidth=bandwidth)))
    the_mt_sdfe <- mt_sdf_est(ts,the_tapers,center=TRUE,pad_factor=pad_factor,delta_t=delta_t)
    plot(the_mt_sdfe$freqs[-1],dB(the_mt_sdfe$sdfe[-1]),
         xlim=c(0,2),xaxs="i",xlab=expression(paste(italic(f),"  (Hz)")),
         ylim=c(-20,70),yaxs="i",ylab="dB",
         typ="l",lwd=0.5,col="gray40",axes=FALSE,
         main="Figure 582(b)")
    x_cc <- 1.92
    y_cc <- 43
    cc <- the_mt_sdfe$cc
    lines(c(x_cc,x_cc),y_cc+c(cc$up,-cc$down),lwd=0.5)
    lines(x_cc+c(-cc$width/2,cc$width/2),c(y_cc,y_cc),lwd=0.5)
    ## reshaping ...
    F_test_results <- thomson_F_test(ts,the_tapers,pad_factor=pad_factor,delta_t=delta_t,center=TRUE)
    for (i in 1:length(reshape_indices))
    {
        temp <- reshaped_mt_sdfe(reshape_indices[i],F_test_results,the_mt_sdfe)
        lines(temp$freqs,dB(temp$sdfe),lwd=2)
    }
    axis(1,at=seq(0,2,0.5))
    axis(1,at=seq(0,2,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(0,100,20),las=2)
    axis(2,at=seq(0,100,10),label=FALSE,tcl=-0.25)
    text(1.97,59.2,"(b)",pos=2)
    box(bty="l")
}

### Figure 582, top plot

fig_582_top(ocean_wave)

### Figure 582, bottom plot

fig_582_bottom(ocean_wave)

### Table 583 ###

wr_uncentered <- Willamette_River_all$V2
X_bar <- mean(wr_uncentered)
wr_centered <- wr_uncentered - X_bar
N <- length(wr_centered)
delta_t <- 1/12
f_hats <- as.vector(find_max_peaks_dse(wr_centered,center=FALSE,delta_t=1/12,N_peaks=2))

### Table 583, 1st row

X_rows_1_2 <- cbind(cos(2*pi*f_hats[1]*(0:(N-1))*delta_t),sin(2*pi*f_hats[1]*(0:(N-1))*delta_t))
lm_row_1 <- lm(wr_uncentered ~ X_rows_1_2)

round(coef(lm_row_1)[1],5)       #  9.82573 
round(sum(resid(lm_row_1)^2),5)  # 87.40153

### Table 583, 2nd row

lm_row_2 <- lm(wr_centered ~ X_row_2 - 1)

round(X_bar,5)                   #  9.82542
round(sum(resid(lm_row_2)^2),5)  # 87.40157

### Table 583, 3rd row

SS_3 <- function(parms,ts=wr_uncentered)
{
    N <- length(ts)
    return(sum((ts - as.vector(cbind(rep(1,N),cos(2*pi*parms[4]*(0:(N-1))*delta_t),sin(2*pi*parms[4]*(0:(N-1))*delta_t)) %*% parms[1:3]))^2))
}
    
X_table_570_row_2 <- cbind(cos(2*pi*f_hats[1]*(0:(N-1))*delta_t),sin(2*pi*f_hats[1]*(0:(N-1))*delta_t))
lm_table_570_row_2 <- lm(wr_centered ~ X_table_570_row_2 - 1)
initial_parms <- c(X_bar,coef(lm_table_570_row_2),f_hats[1])
nlm_row_3 <- nlm(SS_3,initial_parms,gradtol = 1e-8,steptol = 1e-8)

round(nlm_row_3$est[1],5)     #  9.82576
round(SS_3(nlm_row_3$est),5)  # 87.39877

### Table 583, 4th row

SS_4 <- function(parms,ts=wr_centered)
{
    N <- length(ts)
    return(sum((ts - as.vector(cbind(cos(2*pi*parms[3]*(0:(N-1))*delta_t),sin(2*pi*parms[3]*(0:(N-1))*delta_t)) %*% parms[1:2]))^2))
}
    
initial_parms <- c(coef(lm_table_570_row_2),f_hats[1])
nlm_row_4 <- nlm(SS_4,initial_parms,gradtol = 1e-8,steptol = 1e-8)

round(X_bar,5)                #  9.82542
round(SS_4(nlm_row_4$est),5)  # 87.39881

### Table 583, 5th row

X_rows_5_6 <- cbind(sapply(f_hats,function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(f_hats,function(f) sin(2*pi*f*(0:(N-1))*delta_t)))
lm_row_5 <- lm(wr_uncentered ~ X_rows_5_6)

round(coef(lm_row_5)[1],5)       #  9.82564 
round(sum(resid(lm_row_5)^2),5)  # 79.58303

### Table 583, 6th row

lm_row_6 <- lm(wr_centered ~ X_rows_5_6 - 1)

round(X_bar,5)                   #  9.82542
round(sum(resid(lm_row_6)^2),5)  # 79.58305

### Table 583, 7th row

SS_7 <- function(parms,ts=wr_uncentered)
{
    N <- length(ts)
    return(sum((ts - as.vector(cbind(rep(1,N),sapply(parms[6:7],function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(parms[6:7],function(f) sin(2*pi*f*(0:(N-1))*delta_t))) %*% parms[1:5]))^2))
}

X_table_570_row_5 <- cbind(sapply(f_hats,function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(f_hats,function(f) sin(2*pi*f*(0:(N-1))*delta_t)))
lm_table_570_row_5 <- lm(wr_centered ~ X_table_570_row_5 - 1)
initial_parms <- c(X_bar,coef(lm_table_570_row_5),f_hats)
nlm_row_7 <- nlm(SS_7,initial_parms,gradtol = 1e-8,steptol = 1e-8)

round(nlm_row_7$est[1],5)     #  9.82566
round(SS_7(nlm_row_7$est),5)  # 79.54591

### Table 583, 8th row

SS_8 <- function(parms,ts=wr_centered)
{
    N <- length(ts)
    return(sum((ts - as.vector(cbind(sapply(parms[5:6],function(f) cos(2*pi*f*(0:(N-1))*delta_t)),sapply(parms[5:6],function(f) sin(2*pi*f*(0:(N-1))*delta_t))) %*% parms[1:4]))^2))
}

initial_parms <- c(coef(lm_table_570_row_5),f_hats)
nlm_row_8 <- nlm(SS_8,initial_parms,gradtol = 1e-8,steptol = 1e-8)

round(X_bar,5)                #  9.82542
round(SS_8(nlm_row_8$est),5)  # 79.54593
