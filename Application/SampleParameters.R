#### Gibbs sampling ####
set.seed(3)
X3 <- cbind(rep(1,nrow(Data)), Data[,c("TRT", "RSBP_C")], Data[,"TRT"] * Data[,"RSBP_C"])
X3 <- as.matrix(X3*1)

P3 <- ncol(X3)

m0.3 <- rep(b0, P3);
V0.3 <- matrix(B0, P3 , P3);
Pars.3 <- EstimateParametersPG(X = X3, Y = Y12, nBurn = nBurn, nIt = nIt, Start = start, bMu0 = b0, bSigma0 = B0)
Theta.3 <- Transform2Theta(BetaDrawPG = Pars.3[["bDrawPG"]], X = X3, Y = Y12,  Types = Types,
                           Values = Values[["RSBP_C"]], Range = Ranges[["RSBP_C"]],# Ranges = Ranges[["RSBP_C"]], Values = Values[["RSBP_C"]],
                           PriorAlpha = Alpha0, MeasurementLevel = "Continuous")

Decisions.3 <- vector("list", nrow(Types))
for(type in 1:nrow(Types)){
  Decisions.3[[type]] <- EvaluateData(Data = do.call(rbind, Map("-", Theta.3[["mTheta.E"]][[type]],  
                          Theta.3[["mTheta.C"]][[type]])), 
               Weights = weights, Alpha = Alpha, Alternative = c("greater.than", "smaller.than"), 
               Rule = c("Any", "All", "Compensatory"), Truth = Truths, Types = Types[type,])
  }
  

save(Pars.3, Theta.3, Decisions.3,
     file = "Workspaces/Application3.RData")

