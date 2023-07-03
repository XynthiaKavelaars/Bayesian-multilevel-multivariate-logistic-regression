#### 2.  Figure 2. Superiority regions of various decision rules ####
#### 2.1 Figure 2a Single rule ####
setEPS()
postscript("Application/Plots/Decision_rules.eps", family="Times", width=170/25, height = 225/25, pointsize = )
layout(matrix(1:8, nrow=4,byrow=TRUE), height=c(10,1,10,1), width=c(1,1))

par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk7), 
     cex.lab=2.00, cex.axis=2.00, las=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
rect(xleft=0,ybottom=par("usr")[3],xright=par("usr")[2],ytop=par("usr")[4],
     col = "gray90", lwd=1, border="black")
abline(h=0,lty=2)
abline(v=0,lty=2)

par(mgp=c(3.5,1,0), mar=c(6,6,4,2))

plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk7), 
     cex.lab=2.00, cex.axis=2.00, las=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
rect(xleft=0,ybottom=par("usr")[3],xright=par("usr")[2],ytop=par("usr")[4],
    col = "gray90", lwd=1, border="black")
rect(xleft=par("usr")[1], ybottom=0,xright=par("usr")[2],ytop=par("usr")[4],
     col = "gray90", lwd=1, border="black")
abline(h=0,lty=2)
abline(v=0,lty=2)


par(mar=rep(0,01,4))
plot(NULL,
     xlim=c(0,1), xaxt="n",
     ylim=c(0,1), yaxt="n", bty="n")
text(x=0.5, y=0.1, pos = 3, labels=expression("Single (Indep6)"), cex=2.00)


par(mar=rep(0,01,4))
plot(NULL,
     xlim=c(0,1), xaxt="n",
     ylim=c(0,1), yaxt="n", bty="n")
text(x=0.5, y=0.1, pos = 3, labels="Any", cex=2.00)


par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk7), 
     cex.lab=2.00, cex.axis=2.00, las=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
rect(xleft=0,ybottom=0,xright=par("usr")[2],ytop=par("usr")[4],#density=10, 
     col = "gray90", lwd=1, border="black")
abline(h=0,lty=2)
abline(v=0,lty=2)


par(mgp=c(3.5,1,0), mar=c(6,6,4,2))
plot(NULL,
     xlim=c(-1,1), xaxs="i", xlab=expression(delta^Indep6), 
     ylim=c(-1,1), yaxs="i", ylab=expression(delta^Strk), 
     cex.lab=2.00, cex.axis=2.00, las=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
polygon(x=c(-0.25,0.25,1,1,0.25,-0.25),
        y=c(par("usr")[4],par("usr")[3],par("usr")[3],par("usr")[4],par("usr")[4],par("usr")[4]),
        col = "gray90", lwd=1, border="black")
rect(xleft=par("usr")[1], ybottom=par("usr")[3], xright=par("usr")[2], ytop=par("usr")[4])
abline(h=0,lty=2)
abline(v=0,lty=2)

par(mar=rep(0,01,4))
plot(NULL,
     xlim=c(0,1), xaxt="n",
     ylim=c(0,1), yaxt="n", bty="n")
text(x=0.5, y=0.1, pos = 3, labels="All", cex=2.00)

par(mar=rep(0,01,4))
plot(NULL,
     xlim=c(0,1), xaxt="n",
     ylim=c(0,1), yaxt="n", bty="n")
text(x=0.5, y=0.1, pos = 3, labels="Compensatory", cex=2.00)

dev.off()


