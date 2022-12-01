# Presented scenarios
if(any(grepl("x", c(Fixed_m6,Random_m6)))){iTypes_m6 <- 1:nrow(Types)
}else{
  iTypes_m6 <- which(Types[,"Methods"] %in% c("Value", "MvB"))
}

if(any(grepl("x", c(Fixed_NH_m6,Random_NH_m6)))){iTypes_NH_m6 <- 1:nrow(Types)
}else{
  iTypes_NH_m6 <- which(Types[,"Methods"] %in% c("Value", "MvB"))
}


#### Make table ####
# Names of decision rules
Rules <- c("Any", "All", "Compensatory")

# Compute treatment differences
Delta_m6 <- lapply(1:length(iTypes_m6), function(i){
  DeltaMV <- Map("-", Theta_m6[["mTheta.E"]][[i]],Theta_m6[["mTheta.C"]][[i]])
  DeltaW <- lapply(DeltaMV, function(i) i %*% Weights)
  list(DeltaMV = DeltaMV, DeltaW = DeltaW) 
  })

Delta_m6_NH <- lapply(1:length(iTypes_NH_m6), function(i){
  DeltaMV <- Map("-", Theta_NH_m6[["mTheta.E"]][[i]],Theta_NH_m6[["mTheta.C"]][[i]])
  DeltaW <- lapply(DeltaMV, function(i) i %*% Weights)
  list(DeltaMV = DeltaMV, DeltaW = DeltaW)
  })

# Extract point estimates
PointEsts_m6 <- lapply(1:length(iTypes_m6), function(i){
  Delta_H <- c(colMeans(do.call(rbind, Delta_m6[[i]][["DeltaMV"]])),
  colMeans(do.call(rbind, Delta_m6[[i]][["DeltaW"]])))
  names(Delta_H) <- c("delta1", "delta2", "deltaW")
  Delta_NH <- c(colMeans(do.call(rbind, Delta_m6_NH[[i]][["DeltaMV"]])),
               colMeans(do.call(rbind, Delta_m6_NH[[i]][["DeltaW"]])))
  names(Delta_NH) <- c("delta1", "delta2", "deltaW")
  Delta <- rbind(Delta_H, Delta_NH)
  return(Delta)
 # Delta_H
})     

# Extract posterior probabilities
Pops_m6 <- lapply(1:length(iTypes_m6), function(i){
 Pops_H <- c(unlist(Decisions_m6[[i]][["Pop"]][c("PopMultivariate", "PopWeighted")]), unlist( Decisions_m6[[i]][["Decision"]][grep("TS", names(Decisions_m6[[i]][["Decision"]]))]))
 names(Pops_H) <- c("Pop1", "Pop2", "PopW", Rules)
 Pops_NH <- c(unlist(Decisions_NH_m6[[i]][["Pop"]][c("PopMultivariate", "PopWeighted")]), unlist( Decisions_NH_m6[[i]][["Decision"]][grep("TS", names(Decisions_NH_m6[[i]][["Decision"]]))]))
 names(Pops_NH) <- c("Pop1", "Pop2", "PopW", Rules)
 Pops <- rbind(Pops_H, Pops_NH)
 return(Pops)
})  

Res_m6 <-  lapply(1:length(iTypes_m6), function(i){
 Res <- cbind(PointEsts_m6[[i]], Pops_m6[[i]])
 rownames(Res) <- c("H", "NH")
 return(Res)
})

# Indices per (sub)population for later extraction
Ind.ATE <- which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Trial")
names(Ind.ATE) <- Types[which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Trial"),"Methods"]

Ind.CTE.Range_L <- which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_Lo")  
names(Ind.CTE.Range_L) <- Types[which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_Lo"),"Methods"]

Ind.CTE.Range_ML <- which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_MidL")  
names(Ind.CTE.Range_ML) <- Types[which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_MidL"),"Methods"]

Ind.CTE.Range_MH <- which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_MidH")  
names(Ind.CTE.Range_MH) <- Types[which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_MidH"),"Methods"]

Ind.CTE.Range_H <- which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_Hi")  
names(Ind.CTE.Range_H) <- Types[which(Types[,"Methods"] %in% c("Analytical","Empirical", "MvB") & Types[,"Populations"] == "Intra_Hi"),"Methods"]

Ind.CTE.Value_L <- which(Types[,"Methods"] %in% c("Value") & Types[,"Populations"] == "Intra_Lo") 
names(Ind.CTE.Value_L) <- Types[which(Types[,"Methods"] %in% c("Value") & Types[,"Populations"] == "Intra_Lo"),"Methods"]

Ind.CTE.Value_H <- which(Types[,"Methods"] %in% c("Value") & Types[,"Populations"] == "Intra_Hi") 
names(Ind.CTE.Value_H) <- Types[which(Types[,"Methods"] %in% c("Value") & Types[,"Populations"] == "Intra_Hi"),"Methods"]

# Effects per (sub)population
ATE <- lapply(Ind.ATE, function(i){
  y <- Res_m6[[i]]
  return(y)
})
ATE

CTE.Range_L <- lapply(Ind.CTE.Range_L, function(i){
  y <- Res_m6[[i]]
  return(y)
})
CTE.Range_L

CTE.Range_ML <- lapply(Ind.CTE.Range_ML, function(i){
  y <- Res_m6[[i]]
  return(y)
})
CTE.Range_ML

CTE.Range_MH <- lapply(Ind.CTE.Range_MH, function(i){
  y <- Res_m6[[i]]
  return(y)
})
CTE.Range_MH

CTE.Range_H <- lapply(Ind.CTE.Range_H, function(i){
  y <- Res_m6[[i]]
  return(y)
})
CTE.Range_H

CTE.Value_L <- lapply(Ind.CTE.Value_L, function(i){
  y <- Res_m6[[i]]
  return(y)
})
CTE.Value_L

CTE.Value_H <- lapply(Ind.CTE.Value_H, function(i){
  y <- Res_m6[[i]]
   return(y)
})
CTE.Value_H

# Make table
ColNames <- c("delta1", "delta2", "Pop1", "Pop2",  "Any", "All", "deltaW", "PopW","Compensatory")
matApp <- cbind(rbind(ATE[[1]], ATE[["MvB"]][1,],
                      CTE.Range_L[[1]], CTE.Range_L[["MvB"]][1,],
                      CTE.Range_ML[[1]], CTE.Range_ML[["MvB"]][1,],
                      CTE.Range_MH[[1]], CTE.Range_MH[["MvB"]][1,],
                      CTE.Range_H[[1]], CTE.Range_H[["MvB"]][1,],
                      CTE.Value_L[[1]], 
                      CTE.Value_H[[1]]))[,ColNames]
colnames(matApp) <- c(ColNames)
rownames(matApp) <- rep(c("BMMLR", "BMLR", "BMB"), length(unique(Types[,"Populations"])) + 2)[-c(length(unique(Types[,"Populations"]))*3+c(3,6))]

# Format table
tabApp <- cbind.data.frame(rownames(matApp), matrix(NA, nrow = nrow(matApp), ncol = ncol(matApp)))
colnames(tabApp) <- c("Method", "s1", "DeltaMV", "PopMV", "Any", "All", "s2", "DeltaW", "PopW", "Compensatory")
for(i in 1:nrow(matApp)){
  tabApp[i,"DeltaMV"] <- paste0("($", paste0(sprintf("%6.3f", round(matApp[i,c("delta1", "delta2")],3)), collapse = ", "), "$)", collapse = "")
  tabApp[i, "PopMV"] <- paste0("($", paste0(sprintf("%5.3f", round(matApp[i,c("Pop1", "Pop2")],3)), collapse = ", "), "$)", collapse = "")
  
  if(matApp[i,c("Pop1")] > (1 - Alpha / 4) 
     & matApp[i,c("Pop2")] < (Alpha / 4)){tabApp[i,"Any"] <- "$\\bm{>} \\& \\bm{<}$"
  } else if(matApp[i,c("Pop1")] < (Alpha / 4) 
            & matApp[i,c("Pop2")] > (1 - Alpha / 4)){tabApp[i,"Any"] <- "$\\bm{<} \\& \\bm{>}$"
  } else if(max(matApp[i,c("Pop1","Pop2")]) > (1 - Alpha / 4)){tabApp[i,"Any"] <- "$\\bm{>}$"
  } else if(min(matApp[i,c("Pop1","Pop2")]) < (Alpha / 4)){tabApp[i,"Any"] <- "$\\bm{<}$"
  }  else {tabApp[i,"Any"] <- "-"}
  
  if(min(matApp[i,c("Pop1","Pop2")]) > (1 - Alpha / 2)){tabApp[i,"All"] <- "$\\bm{>}$"
  } else if(max(matApp[i,c("Pop1","Pop2")]) < (Alpha / 2)){tabApp[i,"All"] <- "$\\bm{<}$"
  } else {tabApp[i,"All"] <- "-"}
  
  tabApp[i,"DeltaW"] <- paste0("$", sprintf("%6.3f", round(matApp[i,"deltaW"],3)),"$", collapse = "")
  tabApp[i,"PopW"] <- paste0("$", sprintf("%5.3f", round(matApp[i,"PopW"],3)), "$", collapse = "")
  if(matApp[i,c("PopW")] > (1 - Alpha / 2)){tabApp[i,"Compensatory"] <- "$\\bm{>}$"
  } else if(matApp[i,c("PopW")] < (Alpha / 2)){tabApp[i,"Compensatory"] <- "$\\bm{<}$"
  } else {tabApp[i,"Compensatory"] <- "-"}
}

AddToRow.App <- list()
AddToRow.App$pos <- as.list(c(0, seq(0,length(unique(Types[,"Populations"]))*3,3),length(unique(Types[,"Populations"]))*3+2))

# Subsample sizes
DataC.df <- do.call(rbind,DataC)
nT <- lapply(1:length(RangesApp[[1]]), function(pop){
  nE <- nrow(DataC.df[DataC.df[,"trt"] == 1 & DataC.df[,"x"] >= min(RangesApp[["Continuous"]][[pop]]) & DataC.df[,"x"] <= max(RangesApp[["Continuous"]][[pop]]),])
  nC <- nrow(DataC.df[DataC.df[,"trt"] == 0 & DataC.df[,"x"] >= min(RangesApp[["Continuous"]][[pop]]) & DataC.df[,"x"] <= max(RangesApp[["Continuous"]][[pop]]),])
  return(c(nE=nE,nC=nC))})

nE <- paste0("$n_{A} = ", sapply(nT,"[[",1),"$")
nC <- paste0("$n_{C} = ", sapply(nT,"[[",2),"$")
Props <- sapply(1:length(RangesApp[[1]]), function(pop){
  paste0(c(nE[pop], nC[pop]), 
           c(", ", " "), collapse="")})



Command_te <- c(
  " & & ($\\delta^{Strk7}, \\delta^{Indep6}$) & Pop & Any & All & & $\\delta (\\bm{w})$ & Pop & Comp \\\\\n \\midrule \n ",
                 paste0(c(" ", rep(" \\midrule \n ", length(unique(Types[,"Populations"]))+1)),
                       "\\multicolumn{", 2,
                       "}{l}{", c("ATE", 
                                  "CATE - Low range ",
                                  "CATE - Mid-Low range ",
                                  "CATE - Mid-High range ",
                                  "CATE - High range ",
                                  "CATE - Low value ",
                                  "CATE - High value "),
                       "} & \\multicolumn{", ncol(tabApp)-2,"}{l}{", c(Props, " ", " "), "}" 
                       ,
                       rep(" \\\\ \n ", length(unique(Types[,"Populations"]))),
                       c(" ", rep(" \\midrule \n ", length(unique(Types[,"Populations"]))+1))))
                    



AddToRow.App$command <- c(Command_te)
                          
Align.App <- paste0("llp{0.02cm}rrrrp{0.02cm}rrr")

print(xtable(tabApp, 
             caption= paste0("Average (ATE) and conditional average (CATE) treatment effects of the ", 
                             as.character(length(unique(Types[,"Populations"]))+2),
             " specified (sub)populations of the IST3 , including posterior probabilities (Pop) and  superiority ($>$) and inferiority ($<$) conclusions for each decision rule."),
               label="tab:App", align=Align.App,
             digits=c(1,1,rep(3,ncol(tabApp)-1))),
      add.to.row=AddToRow.App, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)


#### Plot predicted probabilities application ####
# Compute data 
xPointsU <- seq(0,42,1)
xPointsEmpU <- seq(xPointsU[1],xPointsU[length(xPointsU)],1)

EstRC_m6 <- lapply(seq(1,nIt, nThin), function(i){
  rbind(if(length(Fixed_m6) > 0){Pars_m6[["Pars"]][[1]][["bDrawPG"]][[i]]}, 
        if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["gDrawPG"]][[i]]})})

thetaAna.E <- lapply(1:(length(xPointsEmpU)-1), function(x) {
  theta <- EstimateThetaAnalytical(EstPars = EstRC_m6, X = do.call(rbind, xDataC_m6),
                                   Trt = 1, RangeX = xPointsEmpU[x+0:1],
                                   Fixed = Fixed_m6, Random = Random_m6)
  return(do.call(rbind, theta))
})

thetaAna.C <- lapply(1:(length(xPointsEmpU)-1), function(x) {
  theta <- EstimateThetaAnalytical(EstPars = EstRC_m6, X = do.call(rbind, xDataC_m6),
                                   Trt = 0, RangeX = xPointsEmpU[x+0:1],
                                   Fixed = Fixed_m6, Random = Random_m6)
  return(do.call(rbind, theta))
})

delta.Ana <- lapply(1:(length(xPointsEmpU)-1), function(x) {
  delta <- thetaAna.E[[x]] - thetaAna.C[[x]]
  mu <- colMeans(delta)
  se <- apply(delta, 2, function(x) sd(x))
  return(list(mu=mu, se=se))
})

# Save and draw plot 
postscript("Application/Plots/Application10_1.eps", width = 6.5, height = 4, horizontal = FALSE, font = "Times")
layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE), heights = c(4,0))

par(mar = c(4,4,2,1), family = "Times")
plot(NULL, xlim = c(0,42), ylim = c(-0.50,0.50), 
     type = "l", main = "Recurring stroke"
     , xlab = "Stroke severity score (NIHSS)", ylab = bquote(theta[A] - theta[C]), las = 1,
     yaxt = "n", xaxs="i", yaxs  ="i")
polygon(c(rev(xPointsEmpU[-1]), xPointsEmpU[-1]), 
        c(rev(sapply(delta.Ana, function(x) x[["mu"]]+x[["se"]])[1,]), sapply(delta.Ana, function(x) x[["mu"]]-x[["se"]])[1,]), col = 'grey', border = NA)
lines(x = xPointsEmpU[-1], y = sapply(delta.Ana, "[[", "mu")[1,], lty = 1, col = "black")

axis(2, at = round(seq(-0.50,0.50,0.10), 2), las=1)
abline(h=0,lty=3)
abline(v=c(6,15,25), lty=2)

par(mar = c(4,4,2,1), family = "Times")     
plot(NULL, xlim = c(0,42), 
     ylim = c(-0.50,0.50), 
     type = "l", main = "Independent living"
     ,  xlab = "Stroke severity scroe (NIHSS)", ylab = bquote(theta[A] - theta[C]), las = 1,
     yaxt = "n", xaxs="i", yaxs  ="i")
polygon(c(rev(xPointsEmpU[-1]), xPointsEmpU[-1]), 
        c(rev(sapply(delta.Ana, function(x) x[["mu"]]+x[["se"]])[2,]), sapply(delta.Ana, function(x) x[["mu"]]-x[["se"]])[2,]), col = 'grey', border = NA)
lines(x = xPointsEmpU[-1], y = sapply(delta.Ana, "[[", "mu")[2,], lty = 1, col = "black")
abline(h=0, lty=3)
abline(v=c(6,15,25), lty=2)
axis(2, at = round(seq(-0.50,0.50,0.10), 2), las=1)

dev.off()

