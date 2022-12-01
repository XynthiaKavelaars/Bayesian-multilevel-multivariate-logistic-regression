#### 2.  Figure 2. Superiority regions of various decision rules ####
#### 2.1 Figure 2a Single rule ####
setEPS()
postscript("Application/Plots/Figure2a_Single.eps", family="Times")
par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk7), 
     cex.lab=2.00, cex.axis=2.00, las=1)
rect(xleft=0,ybottom=par("usr")[3],xright=par("usr")[2],ytop=par("usr")[4],
     col = "gray90", lwd=1, border="black")
abline(h=0,lty=2)
abline(v=0,lty=2)
dev.off()

#### 2.2 Figure 2b Any rule ####
setEPS()
postscript("Application/Plots/Figure2b_Any.eps", family="Times")
par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk7), 
     cex.lab=2.00, cex.axis=2.00, las=1)
rect(xleft=0,ybottom=par("usr")[3],xright=par("usr")[2],ytop=par("usr")[4],
     #density=10, angle = 45,
     col = "gray90", lwd=1, border="black")
rect(xleft=par("usr")[1], ybottom=0,xright=par("usr")[2],ytop=par("usr")[4],
     #density=10, angle = -45,
     col = "gray90", lwd=1, border="black")

abline(h=0,lty=2)
abline(v=0,lty=2)
dev.off()

#### 2.3 Figure 2c All rule ####
setEPS()
postscript("Application/Plots/Figure2c_All.eps", family="Times")
par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk7), 
     cex.lab=2.00, cex.axis=2.00, las=1)
rect(xleft=0,ybottom=0,xright=par("usr")[2],ytop=par("usr")[4],#density=10, 
     col = "gray90", lwd=1, border="black")
abline(h=0,lty=2)
abline(v=0,lty=2)

dev.off()

#### 2.4 Figure 2d Compensatory rule ####
weights <- function(a1, delta1){(-1*a1*delta1)/(1-a1)}
delta1 <- seq(-1,1,0.1)
dist <- 0.6

setEPS()
postscript("Application/Plots/Figure2d_Compensatory.eps", family="Times")
par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk), 
     cex.lab=2.00, cex.axis=2.00, las=1)
polygon(x=c(-0.25,0.25,1,1,0.25,-0.25),
        y=c(par("usr")[4],par("usr")[3],par("usr")[3],par("usr")[4],par("usr")[4],par("usr")[4]),
        col = "gray90", lwd=1, border="black")
rect(xleft=par("usr")[1], ybottom=par("usr")[3], xright=par("usr")[2], ytop=par("usr")[4])
abline(h=0,lty=2)
abline(v=0,lty=2)

dev.off()


