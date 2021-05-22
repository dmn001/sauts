### R code for reproducing content of figures in Chapter 5 ...

earth_100 <- scan("http://faculty.washington.edu/dbp/sauts/Data/earth_100.txt")

### functions used to compute content of figures in Chapter 5 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/slepian_taper.R")

### Figure 139 ###

fig_139 <- function(W,y_lab)
{
    plot(W[,1],W[,2],
         xlim=c(-1/2,1/2),xaxs="i",xlab=expression(italic(t)),
         ylim=c(-1,1),ylab=y_lab,
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 139")
    axis(1,at=seq(-0.5,0.5,0.5),labels=c("-1/2","0","1/2"))
    axis(2,at=c(-1,1),las=2)
    box(bty="l")
}

W_0 <- cbind(c(-1/2,1/2),c(1,1))
W_1 <- cbind(c(-1/2,0,0,1/2),c(-1,-1,1,1))
W_2 <- cbind(c(-1/2,-1/4,-1/4,1/4,1/4,1/2),c(-1,-1,1,1,-1,-1))
W_3 <- cbind(c(-1,-1/4,-1/4,0,0,1/4,1/4,1),c(1,1,-1,-1,1,1,-1,-1))
W_4 <- cbind(c(-1/2,-3/8,-3/8,-1/8,-1/8,1/8,1/8,3/8,3/8,1/2),c(1,1,-1,-1,1,1,-1,-1,1,1))
W_5 <- cbind(c(-1/2,-3/8,-3/8,-1/8,-1/8,0,0,1/8,1/8,3/8,3/8,1/2),c(-1,-1,1,1,-1,-1,1,1,-1,-1,1,1))
W_6 <- cbind(c(-1/2,-3/8,-3/8,-1/4,-1/4,-1/8,-1/8,1/8,1/8,1/4,1/4,3/8,3/8,1/2),c(-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1))
W_7 <- cbind(c(-1/2,-3/8,-3/8,-1/4,-1/4,-1/8,-1/8,0,0,1/8,1/8,1/4,1/4,3/8,3/8,1/2),c(1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1))
W_8 <- cbind(c(-1/2,-7/16,-7/16,-5/16,-5/16,-3/16,-3/16,-1/16,-1/16,1/16,1/16,3/16,3/16,5/16,5/16,7/16,7/16,1/2),c(1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1))

### Figure 139, plots from top tp bottom

fig_139(W_0,expression(W[0]))
fig_139(W_1,expression(W[1]))
fig_139(W_2,expression(W[2]))
fig_139(W_3,expression(W[3]))
fig_139(W_4,expression(W[4]))
fig_139(W_5,expression(W[5]))
fig_139(W_6,expression(W[6]))
fig_139(W_7,expression(W[7]))
fig_139(W_8,expression(W[8]))

### Figures 140 and 141 ###

fig_140 <- function(xs,ys,y_lab,main="Figure 140")
{
    plot(xs,ys,
         xlim=c(-1/2,1/2),xaxs="i",xlab=expression(italic(t)),
         ylim=round(range(ys)),ylab=y_lab,
         typ="l",lwd=0.5,axes=FALSE,
         main=main)
    axis(1,at=seq(-0.5,0.5,0.25),labels=c("-1/2","-1/4","0","1/4","1/2"))
    axis(2,at=round(c(range(ys),mean(range(ys)))),las=2)
    box(bty="l")
}

ts <- seq(-1/2,1/2,length=1025)
cos_in <- cos(2*pi*ts)
cos_out <- 2*sqrt(2)*cos(2*pi*ts)/pi

### Figure 140, plots from top tp bottom

fig_140(c(-1/2,-1/4,-1/4,1/4,1/4,1/2),c(-1,-1,1,1,-1,-1),expression(g[p](t)))
fig_140(c(-1/2,-1/8,-1/8,1/8,1/8,1/2),c(0,0,4,4,0,0),expression(h[p](t)))
fig_140(c(-1/2,-1/8,1/8,1/2),c(-1,1,1,-1),expression(paste(g[p],"*",h[p](t))))

### Figure 141, plots from top tp bottom

fig_140(ts,cos_in,expression(g[p](t)),main="Figure 141")
fig_140(c(-1/2,-1/8,-1/8,1/8,1/8,1/2),c(0,0,4,4,0,0),expression(h[p](t)),main="Figure 141")
fig_140(ts,cos_out,expression(paste(g[p],"*",h[p](t))),main="Figure 141")

### Figure 146 ###

fig_146 <- function(xs,ys,y_lab)
{
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    ## par(mar=c(5.1,4.1,4.1,2.1)) is default (bottom, left, top, right)
    par(mar=c(5.1,5.1,4.1,3.1))
    plot(xs,ys,
         xlim=c(0,1/2),xaxs="i",xlab=expression(italic(f)),
         ylab=y_lab,
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 146")
    axis(1,at=seq(0,0.5,0.125),labels=c("0","1/8","1/4","3/8","1/2"))
    axis(2,at=c(0,1),las=2)
    box(bty="l")
}

### Figure 146, plots from top to bottom

fig_146(c(0,3/8,3/8,1/2),c(1,1,0,0),expression(paste("|",G^(LP),(f),"|"^2)))
fig_146(c(0,1/8,1/8,1/2),c(0,0,1,1),expression(paste("|",G^(HP),(f),"|"^2)))
fig_146(c(0,1/8,1/8,3/8,3/8,1/2),c(0,0,1,1,0,0),expression(paste("|",G^(BP),(f),"|"^2)))

### Figures 149 and 150 ###

fig_149 <- function(lp_filt,bottom_p=FALSE,ts=earth_100-mean(earth_100),main="Figure 149")
{
    if(!bottom_p)
    {
        N_dft  <- 1024
        freqs  <- (0:(N_dft/2))/N_dft
        sgf_lpf <- (abs(fft(c(lp_filt,rep(0,N_dft-3))))^2)[1:(1+N_dft/2)]
        sgf_res <- (abs(fft(c(-lp_filt[1],1-lp_filt[2],-lp_filt[3],rep(0,N_dft-3))))^2)[1:(1+N_dft/2)]
        plot(freqs,sgf_lpf,
             xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
             ylim=c(0,4.05),yaxs="i",ylab=" ",
             typ="l",axes=FALSE,
             main=main)
        lines(freqs,sgf_res,lwd=0.5,col="gray40")
        axis(1,at=seq(0,0.5,0.5))
        axis(1,at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
        axis(2,at=seq(0,4,1),las=2)
    }
    else
    {
        N <- length(ts)
        ts_lp <- lp_filt[1]*ts[1:(N-2)] + lp_filt[2]*ts[2:(N-1)] + lp_filt[3]*ts[3:N]
        res <- ts[2:(N-1)]-ts_lp
        plot(1:(N-2),ts_lp,
             xlim=c(0,100),xlab=expression(italic(t)),
             ylim=c(-30,30),ylab=" ",
             typ="l",lwd=0.5,col="gray40",axes=FALSE,
             main=main)
        points(0:(N-1),ts,pch=3,cex=0.2)
        points(1:(N-2),res,pch=20,cex=0.25)
        axis(1,at=seq(0,100,25))
        axis(2,at=seq(-30,30,15),las=2)
        axis(2,at=seq(-30,30,5),label=FALSE,tcl=-0.25)
    }
    box(bty="l")
}

### Figure 149, top plot

fig_149(c(1/4,1/2,1/4))

### Figure 149, bottom plot

fig_149(c(1/4,1/2,1/4),bottom_p=TRUE)

### Figure 150, top plot

fig_149(2*c(1/4,1/2,1/4),main="Figure 150")

### Figure 150, bottom plot

fig_149(2*c(1/4,1/2,1/4),bottom_p=TRUE,main="Figure 150")

### Figures 153a and 153b ###

fig_153a <- function(sgf,tag,main="Figure 153a")
{
    freqs <- seq(0,0.5,length=length(sgf))
    plot(freqs,dB(sgf),
         xlim=c(0,1/2),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="squared gain  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main=main)
    abline(v=0.1, lty="dotted")
    axis(1,at=c(0,0.5),labels=c("0","1/2"))
    axis(1,at=seq(0,0.5,0.1),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),labels=FALSE,tcl=-0.25)
    text(0.5,10,tag,pos=2)
    box(bty="l")
}

g_I_u <- function(W,K)
{
    results <- rep(0,2*K+1)
    mid <- K+1
    results[mid] <- 2*W
    for(k in 1:K)
    {
        results[mid+k] <- sin(2*pi*W*k)/(pi*k)
        results[mid-k] <- results[mid+k]
    }
    return(results)
}

triangular_cf <- function(K)
{
    results <- rep(0,2*K+1)
    mid <- K+1
    results[mid] <- 1
    for(k in 1:K)
    {
        results[mid+k] <- 1 - k/(K+1)
        results[mid-k] <- results[mid+k]
    }
    results
}

N_fft <- 1024
GK32 <- (abs(fft(c(rep(0,N_fft-65),g_I_u(0.1,32))))^2)[1:((N_fft/2)+1)]
freqs <- seq(0,0.5,length.out=length(GK32))
N_freqs <- length(freqs)
local_mins <- c(FALSE,GK32[-c(1,N_freqs)] < GK32[-c(1,2)] & GK32[-c(1,N_freqs)] < GK32[-c(N_freqs-1,N_freqs)],FALSE) & freqs > 0.1
GK32[local_mins] <- convert_from_dB(-100)
gc <- triangular_cf(32) * g_I_u(0.1,32)
gc <- gc/(sum(gc))
GK32c <- (abs(fft(c(rep(0,N_fft-65),gc)))^2)[1:((N_fft/2)+1)]

### Figure 153a

fig_153a(GK32,expression(paste("|",G[32],"(",italic(f),")|"^2)))

### Figure 153b

fig_153a(GK32c,expression(paste("|",G[32]^(C),"(",italic(f),")|"^2)),main="Figure 153b")

### Figure 154 ###

fig_154_upper <- function(tf,tag)
{
    freqs <- seq(0,0.5,length=length(tf))
    plot(freqs,tf,
         xlim=c(0,1/2),xaxs="i",xlab=expression(italic(f)),
         ylim=c(0,1.5),yaxs="i",ylab="transfer function",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 154")
    abline(v=0.1, lty="dotted")
    axis(1,at=c(0,0.5),labels=c("0","1/2"))
    axis(1,at=seq(0,0.5,0.1),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(0,1.5,0.5),las=2)
    text(0.5,1.25,tag,pos=2)
    box(bty="l")
}


fig_154_lower <- function(sgf,tag)  # identical to fig_153a
{
    freqs <- seq(0,0.5,length=length(sgf))
    plot(freqs,dB(sgf),
         xlim=c(0,1/2),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="squared gain  (dB)",
         typ="l",lwd=0.5,axes=FALSE,
         main="Figure 154")
    abline(v=0.1, lty="dotted")
    axis(1,at=c(0,0.5),labels=c("0","1/2"))
    axis(1,at=seq(0,0.5,0.1),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),labels=FALSE,tcl=-0.25)
    text(0.5,10,tag,pos=2)
    box(bty="l")
}

G_delta <- function(f,delta=0.02,W=0.1)
{
    f <- f %% 1
    if(f > 0.5) f <- f - 1
    f_abs <- abs(f)
    if(f_abs <= W-delta)
        1
    else
        if(f_abs <= W+delta)
            0.5*(1 + cos(pi*(f_abs - W + delta)/(2*delta)))
    else
        0
}

N_fft <- 1024
Gupper <- sapply((0:(N_fft/2))/N_fft,G_delta)

###

g_delta <- function(u,delta=0.02,W=0.1)
{
    if(u==0)
        2*W
    else
    {
        if(abs(u)==1/(4*delta))
            delta*sin(pi*W/(2*delta))
        else
            sin(2*pi*W*u)*cos(2*pi*delta*u)/(pi*(u-16*u^3*delta^2))
    }
}


us <- (-32):32
gus <- sapply(us,g_delta)
Glower <- (abs(fft(c(rep(0,N_fft-65),gus)))^2)[1:((N_fft/2)+1)]
freqs <- seq(0,0.5,length.out=length(GK32))
N_freqs <- length(freqs)
local_mins <- c(FALSE,Glower[-c(1,N_freqs)] < Glower[-c(1,2)] & Glower[-c(1,N_freqs)] < Glower[-c(N_freqs-1,N_freqs)],FALSE) & freqs > 0.1
Glower[local_mins] <- convert_from_dB(-100)


### Figure 154, upper plot

fig_154_upper(Gupper,expression(paste(G[delta],"(",italic(f),")")))

### Figure 154, lower plot

fig_154_lower(Glower,expression(paste("|",G[32],","[delta],"(",italic(f),")|"^2)))

### Figures 156a and 156b ###

fig_156a <- function(sgf_1,sgf_2,tag,main="Figure 156a")
{
    freqs <- seq(0,0.5,length=length(sgf_1))
    plot(freqs,dB(sgf_1),
         xlim=c(0,1/2),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-60,20),yaxs="i",ylab="squared gain  (dB)",
         typ="l",lwd=0.5,col="gray40",axes=FALSE,
         main=main)
    lines(freqs,dB(sgf_2))
    abline(v=0.1, lty="dotted")
    axis(1,at=c(0,0.5),labels=c("0","1/2"))
    axis(1,at=seq(0,0.5,0.1),labels=FALSE,tcl=-0.25)
    axis(2,at=seq(-60,20,20),las=2)
    axis(2,at=seq(-60,20,10),labels=FALSE,tcl=-0.25)
    text(0.5,10,tag,pos=2)
    box(bty="l")
}

gen_sq_gain_func <- function(K)
{
    twoKp1 <- 2*K+1
    temp <- slepian_taper(twoKp1,twoKp1*0.1)
    gK <- temp/sum(temp)
    N_fft <- 1024
    temp <- (abs(fft(c(rep(0,N_fft-twoKp1),gK)))^2)[1:((N_fft/2)+1)]
    freqs <- seq(0,0.5,length=length(temp))
    N_freqs <- length(freqs)
    local_mins <- c(FALSE,temp[-c(1,N_freqs)] < temp[-c(1,2)] & temp[-c(1,N_freqs)] < temp[-c(N_freqs-1,N_freqs)],FALSE) & freqs > 0.1
    temp[local_mins] <- convert_from_dB(-100)
    return(temp)
}

g_I_u <- function(W,K)
{
    results <- rep(0,2*K+1)
    mid <- K+1
    results[mid] <- 2*W
    for(k in 1:K)
    {
        results[mid+k] <- sin(2*pi*W*k)/(pi*k)
        results[mid-k] <- results[mid+k]
    }
    results
}

gen_sq_gain_func_cf <- function(K,delta)
{
    twoKp1 <- 2*K+1
    temp <- slepian_taper(twoKp1,twoKp1*delta)*g_I_u(0.1,K)
    gK <- temp/sum(temp)
    N_fft <- 1024
    temp <- (abs(fft(c(rep(0,N_fft-twoKp1),gK)))^2)[1:((N_fft/2)+1)]
    freqs <- seq(0,0.5,length=length(temp))
    N_freqs <- length(freqs)
    local_mins <- c(FALSE,temp[-c(1,N_freqs)] < temp[-c(1,2)] & temp[-c(1,N_freqs)] < temp[-c(N_freqs-1,N_freqs)],FALSE) & freqs > 0.1
    temp[local_mins] <- convert_from_dB(-100)
    return(temp)
}

GK8 <- gen_sq_gain_func(8)
GK32 <- gen_sq_gain_func(32)

GK32d04 <- gen_sq_gain_func_cf(32,0.04)
GK32d01 <- gen_sq_gain_func_cf(32,0.01)


### Figure 156a

fig_156a(GK8,GK32," ")

### Figure 156b

fig_156a(GK32d04,GK32d01,expression(paste("|",G[32]^(D),"(",italic(f),")|"^2)),main="Figure 156b")
