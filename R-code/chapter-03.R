### R CODE FOR REPRODUCING CONTENT OF FIGURES IN CHAPTER 3 ...

### functions used to compute content of figures in Chapter 3 ...

source("http://faculty.washington.edu/dbp/sauts/R-code/dB.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dirichlet_kernel.R")
source("http://faculty.washington.edu/dbp/sauts/R-code/dft.R")

### BEGINNING OF CODE TO REPRODUCE CONTENT OF FIGURES/TABLES

### Figure 51a ###

fig_51a <- function(m=4,tag=expression(italic(m==4)))
{
    xs <- seq(-pi,pi,length=1301)
    plot(xs,sapply(xs,function(t){ns <- 1:m;1 + 2 * sum(0.9^ns*cos(ns*t))}),
         xlim=c(-pi,pi),xaxs="i",xlab=expression(italic(t)),
         ylim=c(-1,20),ylab=" ",
         typ="l",col="gray40",axes=FALSE,
         main="Figure 51a")
    lines(xs,sapply(xs,function(x,phi=0.9) (1-phi^2)/(1+phi^2-2*phi*cos(x))))
    axis(1,at=seq(-pi,pi,pi),labels=expression(paste("-",pi,sep=""),0,pi))
    axis(2,at=seq(0,20,5),las=2)
    text(pi*0.9,20,tag,pos=2)
    box(bty="l")
}

### Figure 51a, plots from top to bottom

fig_51a()
fig_51a(8,expression(italic(m==8)))
fig_51a(16,expression(italic(m==16)))
fig_51a(32,expression(italic(m==32)))

### Figure 51b ###

fig_51b <- function()
{
    xs <- (-32):32
    plot(xs,dB(0.9^(2*abs(xs))),
         xlim=c(-32.75,32.75),xaxs="i",xlab=expression(italic(f[n])),
         ylim=c(-30,0),ylab="power spectrum  (dB)",
         typ="p",pch=20,cex=0.5,axes=FALSE,
         main="Figure 51b")
    axis(1,at=seq(-32,32,32),labels=expression(paste("-16/",pi,sep=""),0,paste("16/",pi,sep="")))
    axis(2,at=seq(-30,0,10),las=2)
    box(bty="l")
}

### Figure 51b

fig_51b()

### Figure 54 ###

fig_54 <- function(top_p)
{
    amps <- c(1,0.5,0.05,0.1,0.1,0.025,0.01,0.5,0.03,0.02,0.01,0.005)
    mus  <- c(-0.4,0.25,0.75,-3.0,-1.8,-2.2,-3.7,1.6,2.5,2.8,3.4,3.8)
    sds  <- c(0.25,0.2,0.02,0.5,0.5,0.1,0.1,0.8,0.1,0.1,0.1,0.1)
    ts <- seq(-4,4,length=4001)
    g <- sapply(ts,function(t) sum(amps*mapply(dnorm,t,mus,sds)))
    plot(ts,g, #xs,dB(0.9^(2*abs(xs))),
         xlim=c(-4,4),xaxs="i",xlab=" ",
         ylim=c(0,1.7),ylab=" ",
         typ="n",axes=FALSE,
         main="Figure 54")
    if(top_p) lines(ts,g)
    else
    {
        gg <- ts > -1 & ts <= 1
        ts_1_cycle <- ts[gg] 
        gp <- g[gg]
        lines(ts_1_cycle,gp)
        lines(ts_1_cycle-4,gp)
        lines(ts_1_cycle-2,gp)
        lines(ts_1_cycle+2,gp)
        lines(ts_1_cycle+4,gp)
        abline(v=c(-3,3),lty="dotted")
    }
    abline(v=c(-1,1),lty="dotted")
    axis(1, at=seq(-4,4,2), labels=expression(italic(-2*T),italic(-T),0,italic(T),italic(2*T)))
    box(bty="l")
}

### Figure 54, top plot

fig_54(TRUE)

### Figure 54, bottom plot

fig_54(FALSE)

### Figure 56 ###

fig_56 <- function(sigma=1,tag=expression(sigma==1))
{
    t_or_f <- seq(-4,4,0.01)
    g <- sapply(t_or_f,dnorm,0,sigma)
    plot(t_or_f,g,
         xlim=c(-4,4),xaxs="i",xlab=expression(paste(italic(t),"  (light curves) or ",italic(f),"  (dark)")),
         ylim=c(0,1),yaxs="i",ylab=" ",
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 56")
    G <- sapply(t_or_f,function(x,sd=sigma) exp(-2*(pi*x*sd)^2))
    lines(t_or_f,G)
    axis(1,at=seq(-4,4,4))
    axis(1,at=seq(-4,4,1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    text(3,0.9,tag,pos=4)
    box(bty="l")
}

### Figure 56, top plot

fig_56()

### Figure 56, middle plot

fig_56(2,expression(sigma==2))

### Figure 56, bottom plot

fig_56(4,expression(sigma==4))

### Figure 59a ###

fig_59a <- function(a=1,tag=expression(a==1))
{
    t_or_f <- seq(-4,4,0.01)
    g <- sapply(t_or_f,function(x) sqrt(abs(a))*dnorm(a*x))
    plot(t_or_f,g,
         xlim=c(-4,4),xaxs="i",xlab=expression(paste(italic(t),"  (light curves) or ",italic(f),"  (dark)")),
         ylim=c(0,1),yaxs="i",ylab=" ",
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 59a")
    G <- sapply(t_or_f,function(x) exp(-2*(pi*x/a)^2)/sqrt(abs(a)))
    lines(t_or_f,G)
    axis(1,at=seq(-4,4,4))
    axis(1,at=seq(-4,4,1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    text(3,0.9,tag,pos=4)
    box(bty="l")
}

### Figure 59a, top plot

fig_59a()

### Figure 59a, middle plot

fig_59a(2,expression(a==2))

### Figure 59a, bottom plot

fig_59a(4,expression(a==4))

### Figure 59b ###

fig_59b <- function()
{
    xs <- seq(-4,4,0.01)
    plot(xs,exp(-abs(xs)),
         xlim=c(-4,4),xaxs="i",xlab=expression(italic(t)),
         ylim=c(0,1),yaxs="i",ylab=" ",
         typ="l",lty="longdash",axes=FALSE,
         main="Figure 59b")
    polygon(c(-1,-1,1,1),c(0,1,1,0),col="grey",lwd=1.0)
    lines(xs,exp(-abs(xs)), lty="longdash", lwd=1.0)
    ## lines(c(-1,-1,1,1),c(0,1,1,0), lwd=1.0)
    axis(1,at=seq(-4,4,1))
    axis(2,at=c(0,1),las=2)
    box(bty="l")
}

### Figure 59b

fig_59b()

### Figure 64 ###

fig_64 <- function(x_or_ts,psi)
{
    plot(x_or_ts,psi,
         xlim=c(-10,10),xaxs="i",xlab=expression(paste(italic(x)," or ",italic(t))),
         ylab=" ",
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 64")
    abline(v=c(-1,1),lty="dashed")
    axis(1,at=seq(-10,10,5))
    axis(1,at=seq(-10,10,1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    box(bty="l")
}

### following approximation to psi is based on Equations (64) and (65)

N <- 1001
delta_N <- 2/N
xs_minus_1_to_plus_1 <- seq(-1+delta_N/2,1-delta_N/2,delta_N)

big_guy <- diag(rep(4/pi,N))
for(i in 1:(N-1))
    for(j in (i+1):N)
    {
        big_guy[i,j] <- sin(8*(i-j)/N)/(2*pi*(i-j)/N)
        big_guy[j,i] <- big_guy[i,j]
    }

temp <- eigen(big_guy,TRUE)
the_eigenvector <- temp$vectors[,1]
sum(the_eigenvector^2)  # 1
(the_eigenvalue <- temp$values[1])  # 498.4407
ys_minus_1_to_plus_1 <- the_eigenvector*sqrt(the_eigenvalue)
(lam_0 <- delta_N * sum(ys_minus_1_to_plus_1^2))  # 0.9958856

M <- 4505
xs_1_to_10 <- xs_minus_1_to_plus_1[N] + delta_N*(1:M)
ys_1_to_10 <- rep(0,M)
for(i in 1:M)
{
    ys_1_to_10[i] <- delta_N*sum(ys_minus_1_to_plus_1*sin(4*(xs_1_to_10[i]-xs_minus_1_to_plus_1))/(pi*(xs_1_to_10[i]-xs_minus_1_to_plus_1)))/lam_0   
}

xs <- c(-rev(xs_1_to_10),xs_minus_1_to_plus_1,xs_1_to_10)
ys <- c( rev(ys_1_to_10),ys_minus_1_to_plus_1,ys_1_to_10)
delta_N * sum(ys^2)  # 0.9997629 - close to unity, as should be the case ...

### Figure 64

fig_64(xs,ys)

### Figure 68 ###

fig_68_top <- function(t,left_p=TRUE)
{
    plot(1:2,1:2,
         xlim=c(-2,2),xlab=expression(italic(u)),
         ylim=c(0,1.2),yaxs="i",ylab=paste("t =",t),
         typ="n",axes=FALSE,
         main="Figure 68")
    ys_g <- c(0,0.75,0.75,0)
    if(left_p)
    {
        xs_h <- c(-0.25,-0.25,0.25,0.25)
        ys_h <- c(0,1,1,0)
        lines(xs_h + t, ys_h, col="gray", lwd=2.0)
        xs_g <- c(-0.75,-0.75,0.75,0.75)
        lines(xs_g,ys_g,lwd=0.75)
    }
    else if(t>-1 && t<1)
    {
        xs_h <- c(rep(max(-0.25+t,-0.75),2),rep(min(0.25+t,0.75),2))
        polygon(xs_h, ys_g,col="grey")
    }
    axis(1,at=seq(-2,2,2))
    axis(1,at=seq(-2,2,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

fig_68_bot <- function()
{
    xs <- c(-1,-0.5,0.5,1)
    ys <- c(0,0.375,0.375,0)
    plot(xs,ys,
         xlim=c(-1.25,1.25),xaxs="i",xlab=expression(italic(t)),
         ylim=c(0,1.2),yaxs="i",ylab="g*h(t)",
         typ="l",axes=FALSE,
         main="Figure 68")
    axis(1,at=seq(-1,1,1))
    axis(1,at=seq(-1.25,1.25,0.25),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 68, top plot, left-hand column

fig_68_top(-1.25)
fig_68_top(-1)
fig_68_top(-0.75)
fig_68_top(-0.5)
fig_68_top(-0.25)
fig_68_top(0)
fig_68_top(0.25)
fig_68_top(0.5)
fig_68_top(0.75)
fig_68_top(1)
fig_68_top(1.25)

### Figure 68, top plot, right-hand column

fig_68_top(-1.25,left_p=FALSE)
fig_68_top(-1,left_p=FALSE)
fig_68_top(-0.75,left_p=FALSE)
fig_68_top(-0.5,left_p=FALSE)
fig_68_top(-0.25,left_p=FALSE)
fig_68_top(0,left_p=FALSE)
fig_68_top(0.25,left_p=FALSE)
fig_68_top(0.5,left_p=FALSE)
fig_68_top(0.75,left_p=FALSE)
fig_68_top(1,left_p=FALSE)
fig_68_top(1.25,left_p=FALSE)

### Figure 68, bottom plot

fig_68_bot()

### Figure 69 ###

fig_69 <- function(n_conv)
{
    xs <- seq(-2.5,2.5,0.001)
    ys <- sapply(xs,switch(n_conv,
                           function(x) if(abs(x) < 0.5) 1 else 0,
                           function(x) if(abs(x) < 1) 1 - abs(x) else 0,
                           function(x) if(abs(x) <= 0.5) 0.75 - x^2 else (if(abs(x) < 1.5) (3-2*abs(x))^2/8 else 0),
                           function(x) if(abs(x) <= 1) 2/3 - x^2 + abs(x)^3/2 else (if(abs(x) < 2) (2-abs(x))^3/6 else 0)))
    plot(xs,ys,
         xlim=c(-2,2),xlab=expression(italic(t)),
         ylim=c(0,1.2),yaxs="i",ylab=" ",
         typ="l",axes=FALSE,
         main="Figure 69")
    axis(1,at=seq(-2,2,2))
    axis(1,at=seq(-2,2,1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    box(bty="l")
}

### Figure 69, top row

fig_69(1)
fig_69(1)
fig_69(2)

### Figure 69, middle row

fig_69(1)
fig_69(2)
fig_69(3)

### Figure 69, bottom row

fig_69(1)
fig_69(3)
fig_69(4)


### Figure 70 ###

fig_70 <- function(sigma,h_and_conv=FALSE,tag=" ")
{
    plot(1:2,1:2,
         xlim=c(-4,4),xlab=expression(italic(t)),
         ylim=if(h_and_conv) c(-6.1,6.1) else c(0,4.1),yaxs="i",ylab=if(h_and_conv) "h(t) and g*h(t)" else "g(t)",
         typ="n",axes=FALSE,
         main="Figure 70")
    xs <- seq(-6,6,0.01)
    if(h_and_conv)
    {
        h <- sapply(xs,function(t) 5*cos(pi*t/3 + 0.5) + cos(6*pi*t+1.1))
        lines(xs, h, col="gray40", lwd=0.5)
        gh <- sapply(xs,function(t) exp(-(sigma*pi/3)^2/2)*5*cos(pi*t/3 + 0.5) + exp(-(sigma*6*pi)^2/2)*cos(6*pi*t+1.1))
        lines(xs, gh)
        axis(2,at=seq(-5,5,5))
        axis(2,at=seq(-6,6,1),label=FALSE,tcl=-0.25)
    }
    else
    {
        g  <- sapply(xs,dnorm,sd=sigma)
        polygon(xs,g,col="grey")
        axis(2,at=seq(0,4,4),las=2)
        axis(2,at=seq(0,4,1),label=FALSE,tcl=-0.25)
        text(4.1,1.2,tag,pos=2)
    }
    axis(1,at=seq(-4,4,2))
    axis(1,at=seq(-4,4,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 70, plots from top to bottom

fig_70(0.1,tag=expression(sigma == 0.1))
fig_70(0.1,h_and_conv=TRUE)
fig_70(0.25,tag=expression(sigma == 0.25))
fig_70(0.25,h_and_conv=TRUE)
fig_70(0.625,tag=expression(sigma == 0.625))
fig_70(0.625,h_and_conv=TRUE)

### Figure 71 ###

fig_71 <- function(delta,h_and_conv=FALSE,tag=" ")
{
    plot(1:2,1:2,
         xlim=c(-4,4),xlab=expression(italic(t)),
         ylim=if(h_and_conv) c(-6.1,6.1) else c(0,4.1),yaxs="i",ylab=if(h_and_conv) "h(t) and r*h(t)" else "r(t)",
         typ="n",axes=FALSE,
         main="Figure 71")
    xs <- seq(-6,6,0.01)
    if(h_and_conv)
    {
        h <- sapply(xs,function(t) 5*cos(pi*t/3 + 0.5) + cos(6*pi*t+1.1))
        lines(xs, h, col="gray40", lwd=0.5)
        sinc <- function(u) {if(u==0) 1 else sin(pi*u)/(pi*u)}
        rh <- sapply(xs,function(t) 5*sinc(delta/3)*cos(pi*t/3 + 0.5) + sinc(6*delta)*cos(6*pi*t + 1.1))
        lines(xs, rh)
        axis(2,at=seq(-5,5,5))
        axis(2,at=seq(-6,6,1),label=FALSE,tcl=-0.25)
    }
    else
    {
        r  <- sapply(xs,function(t) if(abs(t) <= delta) 1/(2*delta) else 0)
        polygon(xs,r,col="grey")
        axis(2,at=seq(0,4,4),las=2)
        axis(2,at=seq(0,4,1),label=FALSE,tcl=-0.25)
        text(4.1,1.2,tag,pos=2)
    }
    axis(1,at=seq(-4,4,2))
    axis(1,at=seq(-4,4,1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 71, plots from top to bottom

fig_71(1/6,tag=expression(delta == 1/6))
fig_71(1/6,h_and_conv=TRUE)
fig_71(1/4,tag=expression(delta == 1/4))
fig_71(1/4,h_and_conv=TRUE)

### Figure 73 ###

fig_73 <- function()
{
    xs <- seq(-4,4,0.01)
    plot(xs,dnorm(xs),
         xlim=c(-3,3),xlab=expression(italic(t)),
         ylim=c(0,0.45),yaxs="i",ylab=" ",
         typ="l",lwd=1.5,axes=FALSE,
         main="Figure 73")
    e_half_width <- sqrt(2*pi)/2
    v_half_width <- 2*sqrt(3)/2
    a_half_width <- 2*sqrt(pi)/2
    abline(v=e_half_width*c(-1,1),lwd=0.5)
    abline(v=v_half_width*c(-1,1),lty="dotted")
    abline(v=a_half_width*c(-1,1),lty="longdash")
    axis(1,at=seq(-3,3,1))
    axis(2,at=seq(0,0.4,0.4),las=2)
    axis(2,at=seq(0,0.4,0.1),label=FALSE,tcl=-0.25)
    box(bty="l")
}

### Figure 73

fig_73()

### Figure 77 ###

fig_77 <- function(m,tag)
{
    freqs <- seq(-0.5,0.5,0.001)
    plot(freqs,sapply(freqs,dirichlet_kernel,2*m+1),
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-0.275,1.025),yaxs="i",ylab=" ",
         typ="l",axes=FALSE,
         main="Figure 77")
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(-0.5,1,0.1),label=FALSE,tcl=-0.25)
    text(0.48,0.8,tag,pos=2)
    box(bty="l")
}

### Figure 77, plots from top to bottom

fig_77(4,tag=expression(italic(m == 4)))
fig_77(16,tag=expression(italic(m == 16)))
fig_77(64,tag=expression(italic(m == 64)))

### Figure 78 ###

fig_78 <- function(m,tag)
{
    freqs <- seq(-0.5,0.5,0.001)
    Gp_func <- function(f) exp(-10000*(f-0.23)^2) + exp(-10000*(f+0.23)^2) + exp(-10000*(f-0.27)^2) + exp(-10000*(f+0.27)^2)
    gt_seq <- function(t) 2*sqrt(pi)*exp(-(pi*t)^2/10000)*(cos(pi*0.46*t) + cos(pi*0.54*t))/100
    Gpm_func <- function(f,m)
    {
        ts <- (-m):m
        gts <- sapply(ts,gt_seq)
        sum(gts * cos(2*pi*f*ts))
    }
    plot(freqs,sapply(freqs,Gp_func),
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-0.125,1.025),yaxs="i",ylab=" ",
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 78")
    lines(freqs,sapply(freqs,Gpm_func,m))
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(-0.5,1,0.1),label=FALSE,tcl=-0.25)
    text(0.48,0.9,tag,pos=2)
    box(bty="l")
}

### Figure 78, plots from top to bottom

fig_78(4,tag=expression(italic(m == 4)))
fig_78(16,tag=expression(italic(m == 16)))
fig_78(64,tag=expression(italic(m == 64)))

### Figure 79 ###

fig_79 <- function(m,tag)
{
    freqs <- seq(-0.5,0.5,0.001)
    Gp_func <- function(f) if(abs(f) <= 0.25) 1 else 0
    gt_seq <- function(t) if(t == 0) 0.5 else c(0,1,0,-1)[(abs(t)%%4)+1]/(pi*abs(t))
    Gpm_func <- function(f)
    {
        ts <- (-m):m
        gts <- sapply(ts,gt_seq)
        sum(gts * rep(1,2*m+1) * cos(2*pi*f*ts))
    }
    plot(freqs,sapply(freqs,Gp_func),
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-0.125,1.125),yaxs="i",ylab=" ",
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 79")
    lines(freqs,sapply(freqs,Gpm_func))
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(-0.5,1.1,0.1),label=FALSE,tcl=-0.25)
    text(0.48,0.9,tag,pos=2)
    box(bty="l")
}

### Figure 79, plots from top to bottom

fig_79(4,expression(italic(m == 4)))
fig_79(16,expression(italic(m == 16)))
fig_79(64,expression(italic(m == 64)))

### Figure 80 ###

fig_80 <- function(m,tag)
{
    freqs <- seq(-0.5,0.5,0.001)
    plot(freqs,sapply(freqs,dirichlet_kernel,m)^2,
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-0.275,1.025),yaxs="i",ylab=" ",
         typ="l",axes=FALSE,
         main="Figure 80")
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(-0.5,1,0.1),label=FALSE,tcl=-0.25)
    text(0.48,0.8,tag,pos=2)
    box(bty="l")
}

### Figure 80, plots from top to bottom

fig_80(4,tag=expression(italic(m == 4)))
fig_80(16,tag=expression(italic(m == 16)))
fig_80(64,tag=expression(italic(m == 64)))

### Figure 81 ###

fig_81 <- function(m,tag)
{
    freqs <- seq(-0.5,0.5,0.001)
    Gp_func <- function(f) if(abs(f) <= 0.25) 1 else 0
    gt_seq <- function(t) if(t == 0) 0.5 else c(0,1,0,-1)[(abs(t)%%4)+1]/(pi*abs(t))
    Gpm_func <- function(f)
    {
        ts <- (-m):m
        gts <- sapply(ts,gt_seq)
        sum(gts * (1-abs(ts)/m) * cos(2*pi*f*ts))
    }
    plot(freqs,sapply(freqs,Gp_func),
         xlim=c(-0.5,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-0.125,1.125),yaxs="i",ylab=" ",
         typ="l",col="gray40",lwd=0.5,axes=FALSE,
         main="Figure 81")
    lines(freqs,sapply(freqs,Gpm_func))
    axis(1,at=seq(-0.5,0.5,0.5))
    axis(1,at=seq(-0.5,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=c(0,1),las=2)
    axis(2,at=seq(-0.5,1.1,0.1),label=FALSE,tcl=-0.25)
    text(0.48,0.9,tag,pos=2)
    box(bty="l")
}

### Figure 81, plots from top to bottom

fig_81(4,expression(italic(m == 4)))
fig_81(16,expression(italic(m == 16)))
fig_81(64,expression(italic(m == 64)))

### Figure 83a ###

fig_83a <- function(bottom_plot=FALSE,colors=c("gray75","gray25"))
{
    plot(1:2,1:2,
         xlim=c(-1.5,1.5),xlab=expression(italic(f)),
         ylim=c(0,1.1),yaxs="i",ylab=" ",
         typ="n",axes=FALSE,
         main="Figure 83a")
    hbox_x <- c(-0.0625,-0.0625,0,0)
    box_y <- c(0,0.25,0.25,0)
    polygon(hbox_x-0.5,box_y,col=colors[1],border=NA)
    box_x <- c(-0.0625,-0.0625,0.0625,0.0625)
    polygon(box_x+0.75,2*box_y,col=colors[2])
    polygon(box_x-0.5,box_y)
    polygon(box_x,4*box_y)
    abline(v=c(-0.5,0.5),lty="dashed")
    if(bottom_plot)
    {
        polygon(hbox_x-1.5,box_y,col=colors[1],border=NA)
        polygon(hbox_x+0.5,box_y,col=colors[1],border=NA)
        polygon(hbox_x+1.5,box_y,col=colors[1],border=NA)
        polygon(box_x-1.25,2*box_y,col=colors[2])
        polygon(box_x-0.25,2*box_y,col=colors[2])
        polygon(box_x-1.5,box_y)
        polygon(box_x+0.5,box_y)
        polygon(box_x+1.5,box_y)
        polygon(box_x-1.0,4*box_y)
        polygon(box_x+1.0,4*box_y)
    }
    axis(1,at=seq(-1.5,1.5,0.5),label=expression(-3*f[N],-2*f[N],-f[N],0,f[N],2*f[N],3*f[N]))
    box(bty="l")
}

### Figure 83a, top plot

fig_83a()

### Figure 83a, bottom plot

fig_83a(bottom_plot=TRUE)

### alterative version of top plot of Figure 83a (flashier colors!)

fig_83a(colors=c("red","blue"))

### alterative version of bottom plot of Figure 83a (flashier colors!)

fig_83a(bottom_plot=TRUE,colors=c("red","blue"))

### Figure 83b ###

fig_83b <- function(row)
{
    xs <- seq(-1,9,0.001)
    plot(xs, cos((1+2*pi*row)*xs),
         xlim=c(0,8), xlab=expression(italic(t)),
         ylim=c(-1,1), ylab=" ",
         type='l', col="gray40", lwd=0.5, axes=FALSE,
         main="Figure 83b")
    lines(xs, cos(xs), col="black")
    points(0:8,cos(0:8))
    axis(1, at=0:8, cex.axis=1.2)
}

### Figure 83b, plots from top to bottom

fig_83b(1)
fig_83b(2)
fig_83b(3)

### Figure 89 ###

fig_89 <- function(k,N,tag)
{
    xs <- 0:(N-1)
    temp <- (N-1)*0.05
    x_range <- c(-temp,N-1+temp)
    W <- 4/N
    A_mat<- diag(2*W,N)
    tpW <- 2*pi*W
    for(i in 1:(N-1))
        for(j in (i+1):N)
            A_mat[i,j] <- A_mat[j,i] <- sin(tpW*(j-i))/(pi*(j-i))
    plot(xs,-eigen(A_mat,TRUE)$vectors[,k+1],
         xlim=x_range,xaxs="i",xlab=expression(italic(t)),
         ylim=c(-0.44,0.44),yaxs="i",ylab=" ",
         type='p', pch=20, cex=0.2, axes=FALSE,
         main="Figure 89")
    abline(h=0,lwd=0.5)
    axis(1, at=c(0,N-1))
    axis(2,at=seq(-0.4,0.4,0.4),las=2)
    axis(2,at=seq(-0.4,0.4,0.1),label=FALSE,tcl=-0.25)
    text(x_range[2],-0.44+0.9*0.88,tag,pos=2)
    box(bty="l")
}

### Figure 89, left-hand column

fig_89(0,32,expression(italic(k == 0)))
fig_89(1,32,expression(italic(k == 1)))
fig_89(2,32,expression(italic(k == 2)))
fig_89(3,32,expression(italic(k == 3)))

### Figure 89, right-hand column

fig_89(0,99,expression(italic(k == 0)))
fig_89(1,99,expression(italic(k == 1)))
fig_89(2,99,expression(italic(k == 2)))
fig_89(3,99,expression(italic(k == 3)))

### Figure 90 ###

fig_90 <- function(k,N,tag)
{
    N_dft <- 4096
    xs <- (0:(N_dft/2))/N_dft
    W <- 4/N
    A_mat<- diag(2*W,N)
    tpW <- 2*pi*W
    for(i in 1:(N-1))
        for(j in (i+1):N)
            A_mat[i,j] <- A_mat[j,i] <- sin(tpW*(j-i))/(pi*(j-i))
    ys <- dB(abs(dft(c(-eigen(A_mat,TRUE)$vectors[,k+1],rep(0,N_dft-N)))[1:(N_dft/2+1)])^2)
    null_fixer <- function(x,tol=5.0,level=-120.0)
    {
        N <- length(x)
        r <- 2:(N-1)
        x[(x[r-1] > x[r]) & (x[r+1] > x[r]) & (x[r-1] > (x[r]+tol))] <- level
        return(x)
    }
    plot(xs,null_fixer(ys,0.25),
         xlim=c(0,0.5),xaxs="i",xlab=expression(italic(f)),
         ylim=c(-104,24),yaxs="i",ylab=" ",
         type='l', lwd=0.25, axes=FALSE,
         main="Figure 90")
    abline(v=4/N)
    axis(1, at=seq(0,0.5,0.5))
    axis(1, at=seq(0,0.5,0.1),label=FALSE,tcl=-0.25)
    axis(2,at=seq(-100,20,40),las=2)
    axis(2,at=seq(-100,20,10),label=FALSE,tcl=-0.25)
    text(0.5,10,tag,pos=2)
    box(bty="l")
}

### Figure 90, left-hand column

fig_90(0,32,expression(italic(k == 0)))
fig_90(1,32,expression(italic(k == 1)))
fig_90(2,32,expression(italic(k == 2)))
fig_90(3,32,expression(italic(k == 3)))

### Figure 90, right-hand column

fig_90(0,99,expression(italic(k == 0)))
fig_90(1,99,expression(italic(k == 1)))
fig_90(2,99,expression(italic(k == 2)))
fig_90(3,99,expression(italic(k == 3)))

### Figure 91 ###

fig_91 <- function(N,tag)
{
    xs <- 0:(N-1)
    W <- 4/N
    A_mat<- diag(2*W,N)
    tpW <- 2*pi*W
    for(i in 1:(N-1))
        for(j in (i+1):N)
            A_mat[i,j] <- A_mat[j,i] <- sin(tpW*(j-i))/(pi*(j-i))
    plot(xs,eigen(A_mat,TRUE,TRUE)$values,
         xlim=c(0,N-1),xlab=expression(italic(k)),
         ylim=c(0,1),ylab=" ",
         type='p', pch=20, cex=0.5, axes=FALSE,
         main="Figure 91")
    abline(v=8,lwd=0.5)
    axis(1, at=c(0,8,N-1))
    axis(2,at=c(0,1),las=2)
    text(N-1,0.8,tag,pos=2)
    box(bty="l")
}

### Figure 91, top plot

fig_91(32,expression(italic(N == 32)))

### Figure 91, bottom plot

fig_91(99,expression(italic(N == 99)))
