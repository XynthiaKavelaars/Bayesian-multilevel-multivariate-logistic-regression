#### Treatment effects ####
load("Workspaces/Application3.RData")

pReject.3 <- cbind.data.frame(Types[,c("Methods", "Populations")], 
                              do.call(rbind, sapply(Decisions.3, function(x) unlist(x[["Bias"]]), simplify = FALSE, USE.NAMES = TRUE)),
                              do.call(rbind, sapply(Decisions.3, function(x) unlist(x[["Pop"]]), simplify = FALSE, USE.NAMES = TRUE)),
                              do.call(rbind, sapply(Decisions.3, function(x) c(x[["Decision"]][grep(".LS", names(x[["Decision"]]))]), simplify = FALSE, USE.NAMES = TRUE)))

pReject.3[pReject.3 == "MvB"] <- "Reference"
#pReject.3 <- pReject3[which(pReject3[,"Methods"] != "Analytical" & pReject3[,"Populations"] %in% c("Trial", "Lo", "Hi")),]

#### 1 Make Latex table ####

tab.pReject.3 <- data.frame(matrix(NA, nrow = nrow(pReject.3), ncol = 10))

# Layout table, with brackets etc.
for(i in 1:nrow(pReject.3)){
  tab.pReject.3[i,] <-
    unlist(c(pReject.3[i,c("Methods")], NA,
             paste0("(", paste0(sprintf("%6.3f", pReject.3[i,grep("BiasM", colnames(pReject.3))]), collapse = ", "), ")", collapse = ""),
             paste0("(", paste0(sprintf("%5.3f", pReject.3[i,grep("PopM", colnames(pReject.3))]), collapse = ", "), ")", collapse = ""),
             if(max(pReject.3[i,grep("PopM", colnames(pReject.3))]) > (1 - Alpha / 4) 
                & min(pReject.3[i,grep("PopM", names(pReject.3))]) < (Alpha / 4)){"$\\bm{<} \\& \\bm>}$"
             } else if(max(pReject.3[i,grep("PopM", colnames(pReject.3))]) > (1 - Alpha / 4)){"$\\bm{<}$"
             } else if(min(pReject.3[i,grep("PopM", names(pReject.3))]) < (Alpha / 4)){"$\\bm{>}$"
             }  else {"-"},
             if(min(pReject.3[i,grep("PopM", colnames(pReject.3))]) > (1 - Alpha / 2)){"$\\bm{<}$"
             } else if(max(pReject.3[i,grep("PopM", names(pReject.3))]) < (Alpha / 2)){"$\\bm{>}$"
             } else {"-"},
             NA,
             sprintf("%6.3f", pReject.3[i,grep("BiasW", colnames(pReject.3))]),
             sprintf("%5.3f", pReject.3[i,grep("PopW", colnames(pReject.3))]),
             if(pReject.3[i,grep("PopW", colnames(pReject.3))] > (1 - Alpha / 2)){"$\\bm{<}$"
             } else if(pReject.3[i,grep("PopW", colnames(pReject.3))] < (Alpha / 2)){"$\\bm{>}$"} else {"-"}))
}

#### 1.1 Range of blood pressure ####
tab.pReject.Range.3 <- tab.pReject.3[pReject.3[,"Methods"] %in% c("mB", "mLR") & pReject.3[,"Populations"] %in% c("Trial", "Lo", "Hi"),]
AddToRow.pReject.Range.3 <- list()
rows <- sapply(Populations[which(Populations %in% Types.Range[,"Populations"])], function(pop) sum(Types.Range[,"Populations"] == pop))
AddToRow.pReject.Range.3[["pos"]] <- as.list(c(cumsum(c(0, 0, rows[-length(rows)]))), nrow(tab.pReject.Range.3))
AddToRow.pReject.Range.3[["command"]] <- c(paste0("Method  & & \\multicolumn{1}{l}{$\\bm{\\delta}(Bp)$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[1:2], collapse = ""), 
                                                  " & & \\multicolumn{1}{l}{$\\delta(\\bm{w},Bp)$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[3], collapse = ""), " \\\\\n ", collapse = ""),

paste0("\\midrule \n \\multicolumn{4}{l}{\\textbf{", names.Populations[which(Populations %in% Types.Range[,"Populations"])], "} (",
                                                  #"} & \\multicolumn{3}{l}{", 
                                                  range.Populations[[3]][which(Populations %in% Types.Range[,"Populations"])], ") 
                                                  } & \\multicolumn{", ncol(tab.pReject.Range.3) - 4, "}{l}{$n_{H+A} = ", nRange[["RSBP_C"]][which(Populations %in% Types.Range[,"Populations"]),"nE"], "$, $n_{A} = ", 
       nRange[["RSBP_C"]][which(Populations %in% Types.Range[,"Populations"]),"nC"],"$} \\\\\n " , c("", rep("\\midrule \n", length(which(Populations %in% Types.Range[,"Populations"]))-1))
                                                  ),
paste0("\\midrule \n 
        \multicolumn{", ncol(tab.pReject.Range.3),"}{l}{mB = Multivariate Bernoulli analysis} \\\\\n
        \multicolumn{", ncol(tab.pReject.Range.3),"}{l}{mLR= Multivariate logistic regression} \\\\\n")
)

Align.pReject.Range.3 <- paste0("ll", paste0(rep("r", ncol(tab.pReject.Range.3) - 1), collapse = ""))

print(xtable(tab.pReject.Range.3, 
             caption="Average and conditional average treatment effects (ATE and CATE respectively) in the IST-data, by interval of blood pressure.",
             label="tab:Application.Range_rbsp", align=Align.pReject.Range.3),
      add.to.row=AddToRow.pReject.Range.3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)                                   


#### 1.2 Values of blood pressure ####
names.Values.3 <- paste0(rep(c("-", "+"), each = length(Values[["RSBP_C"]]) / 2), c((length(Values[["RSBP_C"]]) / 2):1, 1:(length(Values[["RSBP_C"]]) / 2)), " SD")# (" , 

tab.pReject.Value.3 <- cbind.data.frame(names.Values.3, subset(tab.pReject.3[tab.pReject.3[,"X1"] %in% c("Value"),-1]))

AddToRow.pReject.Value.3 <- list()
AddToRow.pReject.Value.3[["pos"]] <- as.list(c(0))
AddToRow.pReject.Value.3[["command"]] <- c(paste0("Value  & & \\multicolumn{1}{l}{$\\bm{\\delta} (\\bm{Bp})$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[1:2], collapse = ""), 
                                                  " & & \\multicolumn{1}{l}{$\\delta(\\bm{w}, Bp)$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[3], collapse = ""), " \\\\\n \\midrule \n"))

Align.pReject.Value.3 <- paste0("ll", paste0(rep("r", ncol(tab.pReject.Value.3) - 1), collapse = ""))

print(xtable(tab.pReject.Value.3, 
             caption="Conditional average treatment effects in the IST-data, by value of blood pressure.",
             label="tab:Application.Value_rbsp", align=Align.pReject.Value.3),
       add.to.row=AddToRow.pReject.Value.3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE) 




#### 2 Make Figure ####
StandardizedData <- cbind(Data[,c("TRT")], (Data[,"RSBP"] - mean(Data[,"RSBP"])) / sd(Data[,"RSBP"]), Y12)
colnames(StandardizedData) <- c("TRT", "RSBP", colnames(Y12))


xPoints <- seq(-3,3,1)
thetaMvb.E <- lapply(1:(length(xPoints)-1), function(x) {
  ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 1)
  phi <- rdirichlet(1e4, colSums(Y12[ind,]) + Alpha0)
  theta <- cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)]))
  return(theta)
})

thetaMvb.C <- lapply(1:(length(xPoints)-1), function(x) {
  ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 0)
  phi <- rdirichlet(1e4, colSums(Y12[ind,]) + Alpha0)
  theta <- cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)]))
  return(theta)
})

delta.Mvb <- lapply(1:(length(xPoints)-1), function(x) {
  delta <- thetaMvb.E[[x]] - thetaMvb.C[[x]]
  mu <- colMeans(delta)
  se <- apply(delta, 2, function(x) sd(x))
  return(list(mu=mu, se=se))
})


thetaEmp.E <-  lapply(1:(length(xPoints)-1), function(x) {
    theta <- do.call(rbind, lapply(1:nIt, function(i){
      ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 1)
      phi <- exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]) / rowSums(exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]))
      theta <- colMeans(cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)])))
    }))
    theta    })

thetaEmp.C <-  lapply(1:(length(xPoints)-1), function(x) {
 theta <- do.call(rbind, lapply(1:nIt, function(i){
    ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 0)
    phi <- exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]) / rowSums(exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]))
    theta <- colMeans(cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)])))
 }))
theta    })


delta.Emp <- lapply(1:(length(xPoints)-1), function(x) {
  delta <- thetaEmp.E[[x]] - thetaEmp.C[[x]]
  mu <- colMeans(delta)
  se <- apply(delta, 2, function(x) sd(x))
  return(list(mu=mu, se=se))
})

#thetaAna.E <- lapply(1:(length(xPoints)-1), function(x) {
#  theta <- EstimateThetaAnalytical(BetaPG = Pars.3[["bDrawPG"]], X = X3, Trt = 1, RangeX = xPoints[x+0:1] * sdBp)
#theta
#})

#thetaAna.C <- lapply(1:(length(xPoints)-1), function(x) {
#  theta <- EstimateThetaAnalytical(BetaPG = Pars.3[["bDrawPG"]], X = X3, Trt = 0, RangeX = xPoints[x+0:1] * sdBp)
#theta
#})

#delta.Ana <- lapply(1:(length(xPoints)-1), function(x) {
#  delta <- do.call(rbind, Map("-", thetaAna.E[[x]], thetaAna.C[[x]]))
#  mu <- colMeans(delta)
#  se <- apply(delta, 2, function(x) sd(x))
#  return(list(mu=mu, se=se))
#})

save(thetaEmp.E, thetaEmp.C, #thetaAna.E, #thetaAna.C, 
     delta.Emp, # delta.Ana, 
     file = "Workspaces/Application_figure.RData")

load("Workspaces/Application_figure.RData")

xError1 <- xPoints[-length(xPoints)] + 2/4 * diff(xPoints)
xError2 <- xPoints[-length(xPoints)] + 2/4 * diff(xPoints)
xError3 <- xPoints[-length(xPoints)] + 2/4 * diff(xPoints)


postscript("Plots/Application.eps", width = 6.5, height = 8, horizontal = FALSE, font = "Times")
layout(matrix(c(1,1,2,3,4,4,5,6), nrow = 4, byrow = TRUE), heights = c(1,5,1,5))

par(mar = rep(0.01, 4))
plot(NULL, xlim = c(0,1), ylim = c(0,1), xaxt="n", yaxt="n", bty="n")
text(x = 0.5, y = 0.25, labels = "Recurrent stroke", pos = 3, cex = 2)

par(mar = c(4,4,2,1), family = "Times")
plot(NULL, xlim = c(RoundChoose(min(xPoints * sdBp + meanBp),10,0), RoundChoose(max(xPoints * sdBp + meanBp),10,1)), ylim = c(-0.3,0.3), 
     type = "l", main = "Stratified analysis" #"Recurring stroke"
     , xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
for(i in 1:(length(xPoints)-1)){
  segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Mvb[[i]][["mu"]][1], y1 = delta.Mvb[[i]][["mu"]][1], lty = 1)
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Emp[[i]][["mu"]][1], y1 = delta.Emp[[i]][["mu"]][1], lty = 2)
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Ana[[i]][["mu"]][1], y1 = delta.Ana[[i]][["mu"]][1], lty = 3)
  arrows(x0=xError1[i] * sdBp + meanBp, y0=delta.Mvb[[i]][["mu"]][1] - delta.Mvb[[i]][["se"]][1], x1=xError1[i] * sdBp + meanBp, y1=delta.Mvb[[i]][["mu"]][1] + delta.Mvb[[i]][["se"]][1], code=3, angle=90, length=0.05)
  #arrows(x0=xError2[i] * sdBp + meanBp, y0=delta.Emp[[i]][["mu"]][1] - delta.Emp[[i]][["se"]][1], x1=xError2[i] * sdBp + meanBp, y1=delta.Emp[[i]][["mu"]][1] + delta.Emp[[i]][["se"]][1], code=3, angle=90, length=0.05, lty=2)
  #arrows(x0=xError3[i] * sdBp + meanBp, y0=delta.Ana[[i]][["mu"]][1] - delta.Ana[[i]][["se"]][1], x1=xError3[i] * sdBp + meanBp, y1=delta.Ana[[i]][["mu"]][1] + delta.Ana[[i]][["se"]][1], code=3, angle=90, length=0.05,lty=3)
  
}
abline(h=0,lty=3)

plot(NULL, xlim = c(RoundChoose(min(xPoints * sdBp + meanBp),10,0), RoundChoose(max(xPoints * sdBp + meanBp),10,1)), ylim = c(-0.3,0.3), 
     type = "l", main = "Logistic regression" #"Recurring stroke"
     , xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
for(i in 1:(length(xPoints)-1)){
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Mvb[[i]][["mu"]][1], y1 = delta.Mvb[[i]][["mu"]][1], lty = 1)
  segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Emp[[i]][["mu"]][1], y1 = delta.Emp[[i]][["mu"]][1], lty = 1)
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Ana[[i]][["mu"]][1], y1 = delta.Ana[[i]][["mu"]][1], lty = 3)
  #arrows(x0=xError1[i] * sdBp + meanBp, y0=delta.Mvb[[i]][["mu"]][1] - delta.Mvb[[i]][["se"]][1], x1=xError1[i] * sdBp + meanBp, y1=delta.Mvb[[i]][["mu"]][1] + delta.Mvb[[i]][["se"]][1], code=3, angle=90, length=0.05)
  arrows(x0=xError2[i] * sdBp + meanBp, y0=delta.Emp[[i]][["mu"]][1] - delta.Emp[[i]][["se"]][1], x1=xError2[i] * sdBp + meanBp, y1=delta.Emp[[i]][["mu"]][1] + delta.Emp[[i]][["se"]][1], code=3, angle=90, length=0.05, lty=1)
  #arrows(x0=xError3[i] * sdBp + meanBp, y0=delta.Ana[[i]][["mu"]][1] - delta.Ana[[i]][["se"]][1], x1=xError3[i] * sdBp + meanBp, y1=delta.Ana[[i]][["mu"]][1] + delta.Ana[[i]][["se"]][1], code=3, angle=90, length=0.05,lty=3)
  
}
abline(h=0,lty=3)

par(mar = rep(0.01, 4))
plot(NULL, xlim = c(0,1), ylim = c(0,1), xaxt="n", yaxt="n", bty="n")
text(x = 0.5, y = 0.25, labels = "Dependency", pos = 3, cex = 2)

par(mar = c(4,4,2,1), family = "Times")     
plot(NULL, xlim = c(RoundChoose(min(xPoints * sdBp + meanBp),10,0), RoundChoose(max(xPoints * sdBp + meanBp),10,1)), ylim = c(-0.3,0.3), 
     type = "l", main = "Stratified analysis" #"Dependency"
     ,  xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
for(i in 1:(length(xPoints)-1)){
  segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Mvb[[i]][["mu"]][2], y1 = delta.Mvb[[i]][["mu"]][2], lty = 1)
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Emp[[i]][["mu"]][2], y1 = delta.Emp[[i]][["mu"]][2], lty = 2)
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Ana[[i]][["mu"]][2], y1 = delta.Ana[[i]][["mu"]][2], lty = 3)
  arrows(x0=xError1[i] * sdBp + meanBp, y0=delta.Mvb[[i]][["mu"]][2] - delta.Mvb[[i]][["se"]][2], x1=xError1[i] * sdBp + meanBp, y1=delta.Mvb[[i]][["mu"]][2] + delta.Mvb[[i]][["se"]][2], code=3, angle=90, length=0.05)
  #arrows(x0=xError2[i] * sdBp + meanBp, y0=delta.Emp[[i]][["mu"]][2] - delta.Emp[[i]][["se"]][2], x1=xError2[i] * sdBp + meanBp, y1=delta.Emp[[i]][["mu"]][2] + delta.Emp[[i]][["se"]][2], code=3, angle=90, length=0.05, lty=2)
  #arrows(x0=xError3[i] * sdBp + meanBp, y0=delta.Ana[[i]][["mu"]][2] - delta.Ana[[i]][["se"]][2], x1=xError3[i] * sdBp + meanBp, y1=delta.Ana[[i]][["mu"]][2] + delta.Ana[[i]][["se"]][2], code=3, angle=90, length=0.05,lty=3)
}
abline(h=0, lty=3)

plot(NULL, xlim = c(RoundChoose(min(xPoints * sdBp + meanBp),10,0), RoundChoose(max(xPoints * sdBp + meanBp),10,1)), ylim = c(-0.3,0.3), 
     type = "l", main = "Logistic regression" #"Dependency"
      , xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
for(i in 1:(length(xPoints)-1)){
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Mvb[[i]][["mu"]][2], y1 = delta.Mvb[[i]][["mu"]][2], lty = 1)
  segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Emp[[i]][["mu"]][2], y1 = delta.Emp[[i]][["mu"]][2], lty = 1)
  #segments(x0 = xPoints[i] * sdBp + meanBp, x1= xPoints[i+1] * sdBp + meanBp, y0 = delta.Ana[[i]][["mu"]][2], y1 = delta.Ana[[i]][["mu"]][2], lty = 3)
  #arrows(x0=xError1[i] * sdBp + meanBp, y0=delta.Mvb[[i]][["mu"]][2] - delta.Mvb[[i]][["se"]][2], x1=xError1[i] * sdBp + meanBp, y1=delta.Mvb[[i]][["mu"]][2] + delta.Mvb[[i]][["se"]][2], code=3, angle=90, length=0.05)
  arrows(x0=xError2[i] * sdBp + meanBp, y0=delta.Emp[[i]][["mu"]][2] - delta.Emp[[i]][["se"]][2], x1=xError2[i] * sdBp + meanBp, y1=delta.Emp[[i]][["mu"]][2] + delta.Emp[[i]][["se"]][2], code=3, angle=90, length=0.05, lty=1)
  #arrows(x0=xError3[i] * sdBp + meanBp, y0=delta.Ana[[i]][["mu"]][2] - delta.Ana[[i]][["se"]][2], x1=xError3[i] * sdBp + meanBp, y1=delta.Ana[[i]][["mu"]][2] + delta.Ana[[i]][["se"]][2], code=3, angle=90, length=0.05,lty=3)
}
abline(h=0, lty=3)

#par(mar = rep(0.01, 4))
#plot.new()
#legend(x = "center", legend = c("Reference", "Empirical", "Analytical"), lty = 1:3, ncol = 3) #col = c("black", "red", "darkgreen"))
dev.off()


#### Traceplot ####
postscript("Plots/Traceplot.eps", width = 6.5, height = 8, horizontal = FALSE, font = "Times", family = "sans")
layout(matrix(1:20,nrow=5,ncol=4,byrow=FALSE), heights = c(1,6,6,6,6), widths = c(1,6,6,6))
par(mar=rep(0.01,4))
plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")

plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text(expression(beta[1]^q), x=0.5, y=0.5)

plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text(expression(beta[2]^q), x=0.5, y=0.5)

plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text(expression(beta[3]^q), x=0.5, y=0.5)

plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text(expression(beta[4]^q), x=0.5, y=0.5)

plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text("q=1", x=0.5, y=0.5)

par(mar=c(3,3,1,1))
for(i in 1:4){
  plot(as.vector(Pars.3[["ChainsPG"]][[1]][,i]), type = "l", xlab = "Iterations", ylab = " ", col = "#333333", las=1)
  lines(as.vector(Pars.3[["ChainsPG"]][[2]][,i]), col = "#C0C0C0")  
}

par(mar=rep(0.01,4))
plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text("q=2", x=0.5, y=0.5)


par(mar=c(3,3,1,1))
for(i in 5:8){
  plot(as.vector(Pars.3[["ChainsPG"]][[1]][,i]), type = "l", xlab = "Iterations", ylab = " ", col = "#333333", las=1)
  lines(as.vector(Pars.3[["ChainsPG"]][[2]][,i]), col = "#C0C0C0")  
}

par(mar=rep(0.01,4))
plot(NULL, xlim=c(0,1), ylim = c(0,1), bty="n", xaxt="n", yaxt="n")
text("q=3", x=0.5, y=0.5)

par(mar=c(3,3,1,1))
for(i in 9:12){
  plot(as.vector(Pars.3[["ChainsPG"]][[1]][,i]), type = "l", xlab = "Iterations", ylab = " ", col = "#333333", las=1)
  lines(as.vector(Pars.3[["ChainsPG"]][[2]][,i]), col = "#C0C0C0")  
}

dev.off()
#postscript("Plots/Application.eps", width = 6.5, height = 4, horizontal = FALSE, font = "Times")
#layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE), heights = c(5,2))
#par(mar = c(4,4,2,1), family = "Times")
#plot(x = xPoints[-1] * sdBp + meanBp, y = delta.yPoints[,1], xlim = c(RoundChoose(min(xPoints[-1] * sdBp + meanBp),10,0), RoundChoose(max(xPoints[-1] * sdBp + meanBp),10,1)), ylim = c(-0.5,0.5), 
#     type = "l", main = "Recurring stroke", xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
#lines(x = xPoints[-1] * sdBp + meanBp, y = delta.Emp[,1], lty = 2)
#lines(x = xPoints[-1] * sdBp + meanBp, y = delta.Ana[,1], lty = 3)
#  plot(x = xPoints[-1] * sdBp + meanBp, y = delta.yPoints[,2], xlim = c(RoundChoose(min(xPoints[-1] * sdBp + meanBp),10,0), RoundChoose(max(xPoints[-1] * sdBp + meanBp),10,1)), ylim = c(-0.5,0.5), 
#       type = "l", main = "Dependency",  xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
#  lines(x = xPoints[-1] * sdBp + meanBp, y = delta.Emp[,2], lty = 2)
#  lines(x = xPoints[-1] * sdBp + meanBp, y = delta.Ana[,2], lty = 3)
#  
#  par(mar = rep(0.01, 4))
#  plot.new()
#  legend(x = "center", legend = c("Reference", "Empirical", "Analytical"), lty = 1:3, ncol = 3) #col = c("black", "red", "darkgreen"))
#  dev.off()
  
  #### Descriptives ####
T0_freq <- prop.table(table(Data[Data$TRT == 0, c("STRK14", "DEP6")]))
T1_freq <- prop.table(table(Data[Data$TRT == 1, c("STRK14", "DEP6")]))
T0.tab <- addmargins(T0_freq)
T1.tab <- addmargins(T1_freq)
Prop.tab <- cbind(T0.tab, NA, T1.tab) 

Rho_T0 <- ComputeRho(rev(c(T0_freq)))
Rho_T1 <- ComputeRho(rev(c(T1_freq)))


T0_freq_Lo <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C < -1 * sdBp & Data$RSBP_C > -5 * sdBp, c("STRK14", "DEP6")]))
T1_freq_Lo <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C < -1 * sdBp & Data$RSBP_C > -5 * sdBp, c("STRK14", "DEP6")]))
T0.tab_Lo <- addmargins(T0_freq_Lo)
T1.tab_Lo <- addmargins(T1_freq_Lo)
Prop.tab_Lo <- cbind(T0.tab_Lo, NA, T1.tab_Lo) 

Rho_T0_Lo <- ComputeRho(rev(c(T0_freq_Lo)))
Rho_T1_Lo <- ComputeRho(rev(c(T1_freq_Lo)))


T0_freq_Hi <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C > 1 * sdBp & Data$RSBP_C < 5 * sdBp, c("STRK14", "DEP6")]))
T1_freq_Hi <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C > 1 * sdBp & Data$RSBP_C < 5 * sdBp, c("STRK14", "DEP6")]))
T0.tab_Hi <- addmargins(T0_freq_Hi)
T1.tab_Hi <- addmargins(T1_freq_Hi)
Prop.tab_Hi <- cbind(T0.tab_Hi, NA, T1.tab_Hi) 

Rho_T0_Hi <- ComputeRho(rev(c(T0_freq_Hi)))
Rho_T1_Hi <- ComputeRho(rev(c(T1_freq_Hi)))


T0_freq_Lo1 <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C < 0, c("STRK14", "DEP6")]))
T1_freq_Lo1 <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C < 0, c("STRK14", "DEP6")]))
T0.tab_Lo1 <- addmargins(T0_freq_Lo1)
T1.tab_Lo1 <- addmargins(T1_freq_Lo1)
Prop.tab_Lo1 <- cbind(T0.tab_Lo1, NA, T1.tab_Lo1) 

Rho_T0_Lo1 <- ComputeRho(rev(c(T0_freq_Lo1)))
Rho_T1_Lo1 <- ComputeRho(rev(c(T1_freq_Lo1)))


T0_freq_Hi1 <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C > 0, c("STRK14", "DEP6")]))
T1_freq_Hi1 <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C > 0, c("STRK14", "DEP6")]))
T0.tab_Hi1 <- addmargins(T0_freq_Hi1)
T1.tab_Hi1 <- addmargins(T1_freq_Hi1)
Prop.tab_Hi1 <- cbind(T0.tab_Hi1, NA, T1.tab_Hi1) 

Rho_T0_Hi1 <- ComputeRho(rev(c(T0_freq_Hi1)))
Rho_T1_Hi1 <- ComputeRho(rev(c(T1_freq_Hi1)))



T0_freq_Lo2 <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C < -1 * sdBp & Data$RSBP_C > -2 * sdBp, c("STRK14", "DEP6")]))
T1_freq_Lo2 <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C < -1 * sdBp & Data$RSBP_C > -2 * sdBp, c("STRK14", "DEP6")]))
T0.tab_Lo2 <- addmargins(T0_freq_Lo2)
T1.tab_Lo2 <- addmargins(T1_freq_Lo2)
Prop.tab_Lo2 <- cbind(T0.tab_Lo2, NA, T1.tab_Lo2) 

Rho_T0_Lo2 <- ComputeRho(rev(c(T0_freq_Lo2)))
Rho_T1_Lo2 <- ComputeRho(rev(c(T1_freq_Lo2)))


T0_freq_Hi2 <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C > 1 * sdBp & Data$RSBP_C < 2 * sdBp, c("STRK14", "DEP6")]))
T1_freq_Hi2 <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C > 1 * sdBp & Data$RSBP_C < 2 * sdBp, c("STRK14", "DEP6")]))
T0.tab_Hi2 <- addmargins(T0_freq_Hi2)
T1.tab_Hi2 <- addmargins(T1_freq_Hi2)
Prop.tab_Hi2 <- cbind(T0.tab_Hi2, NA, T1.tab_Hi2) 

Rho_T0_Hi2 <- ComputeRho(rev(c(T0_freq_Hi2)))
Rho_T1_Hi2 <- ComputeRho(rev(c(T1_freq_Hi2)))

T0_freq_Lo3 <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C < -2 * sdBp & Data$RSBP_C > -3 * sdBp, c("STRK14", "DEP6")]))
T1_freq_Lo3 <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C < -2 * sdBp & Data$RSBP_C > -3 * sdBp, c("STRK14", "DEP6")]))
T0.tab_Lo3 <- addmargins(T0_freq_Lo3)
T1.tab_Lo3 <- addmargins(T1_freq_Lo3)
Prop.tab_Lo3 <- cbind(T0.tab_Lo3, NA, T1.tab_Lo3) 

Rho_T0_Lo3 <- ComputeRho(rev(c(T0_freq_Lo3)))
Rho_T1_Lo3 <- ComputeRho(rev(c(T1_freq_Lo3)))


T0_freq_Hi3 <- prop.table(table(Data[Data$TRT == 0 & Data$RSBP_C > 2 * sdBp & Data$RSBP_C < 3 * sdBp, c("STRK14", "DEP6")]))
T1_freq_Hi3 <- prop.table(table(Data[Data$TRT == 1 & Data$RSBP_C > 2 * sdBp & Data$RSBP_C < 3 * sdBp, c("STRK14", "DEP6")]))
T0.tab_Hi3 <- addmargins(T0_freq_Hi3)
T1.tab_Hi3 <- addmargins(T1_freq_Hi3)
Prop.tab_Hi3 <- cbind(T0.tab_Hi3, NA, T1.tab_Hi3) 

Rho_T0_Hi3 <- ComputeRho(rev(c(T0_freq_Hi3)))
Rho_T1_Hi3 <- ComputeRho(rev(c(T1_freq_Hi3)))