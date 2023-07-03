#### 1. Multilevel logistic regression ####
Fixed_m6 <- c("x", "Trt_x")             # Names of fixed predictors
Random_m6 <- c("Intercept", "Trt")      # Names of random predictors

pF_m6 <- length(Fixed_m6)               # Number of fixed predictors
pR_m6 <- length(Random_m6)              # Number of random predictors
bMu0_m6 <- rep(0, pF_m6)                # Prior mean of fixed predictors
gMu0_m6 <- rep(0, pR_m6)                # Prior mean of random predictors

bSigma0_m6 <- diag(1e-1, pF_m6)         # Prior precision matrix of fixed predictors
gSigma0_m6 <- diag(1e-1, pR_m6)         # Prior precision matrix of random predictors

nu0_m6 <- -pR_m6 - 1                      # Prior degrees of freedom covariance matrix random effects
Tau0_m6 <- diag(0, pR_m6)               # Prior covariance matrix of covariance matrix random effects

# Covariate data
xDataC_m6 <- lapply(which(nJ > 0), function(j){
  y <- as.matrix(cbind(1, DataC[[j]][,c(xVars)]))
  colnames(y) <- c("Intercept", "Trt", "x", "Trt_x")
  z <- y[,c(Fixed_m6, Random_m6),drop=FALSE]
  return(z)})
set.seed(2909)

# Sample parameters
Pars_m6 <- EstimateParameters(X = xDataC_m6, Y = yDataC, Fixed = Fixed_m6, Random = Random_m6, 
                              nBurn = nBurn, nIt = nIt, Start = StartVals, 
                              bMu0 = bMu0_m6, bSigma0 = bSigma0_m6, 
                              gMu0 = gMu0_m6, gSigma0 = gSigma0_m6,
                              nu0 = nu0_m6, Tau0 = Tau0_m6, nChain = nChain, ReturnThinned = TRUE, nThin = nThin)

# Diagnose MCMC samples - fixed parameters (if applicable)
if(length(Fixed_m6) > 0){
  Chains.bDraw_m6 <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars_m6[["Pars"]][[chain]][["bDrawPG"]]), nrow = nIt/nThin, ncol = length(Fixed_m6) * Q, byrow = TRUE)[,-((length(Fixed_m6) * (Q-1)) + 1:(length(Fixed_m6)))]))) 
  ESS.bDraw_m6 <- effectiveSize(Chains.bDraw_m6[[1]])
  AC.bDraw_m6 <- autocorr.diag(Chains.bDraw_m6[[1]], lags = 1:100)
  Fit.bDraw_m6 <- summary(Chains.bDraw_m6)$statistics
  Convergence.bDraw_m6 <- gelman.diag(Chains.bDraw_m6)$mpsrf
  mESS.bDraw_m6 <- lapply(1:nChain, function(chain) multiESS(Chains.bDraw_m6[[chain]]))
  mcESS.bDraw_m6 <- lapply(1:nChain, function(chain) ess(Chains.bDraw_m6[[chain]]))
  mcse.bDraw_m6 <- lapply(1:nChain, function(chain) mcse.mat(Chains.bDraw_m6[[chain]]))
}

# Diagnose MCMC samples - random parameters (if applicable)
if(length(Random_m6) > 0){
  Chains.tDraw_m6 <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(lapply(1:(nIt/nThin), function(i) lapply(1:(Q-1), function(q) {x <- Pars_m6[["Pars"]][[chain]][["tauDrawPG"]][[i]][[q]]; x[lower.tri(x, diag = TRUE)]}))), nrow = nIt / nThin, byrow=TRUE))))
  Chains.gDraw_m6 <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars_m6[["Pars"]][[chain]][["gDrawPG"]]), nrow = nIt/nThin, ncol = length(Random_m6) * Q, byrow = TRUE)[,-((length(Random_m6) * (Q-1)) + 1:(length(Random_m6)))])))
  ESS.gDraw_m6 <- effectiveSize(Chains.gDraw_m6[[1]])
  ESS.tDraw_m6 <- effectiveSize(Chains.tDraw_m6[[1]])
  AC.gDraw_m6 <- autocorr.diag(Chains.gDraw_m6[[1]], lags = c(1:10))
  AC.tDraw_m6 <- autocorr.diag(Chains.tDraw_m6[[1]], lags = c(1:10))
  Fit.gDraw_m6 <- summary(Chains.gDraw_m6[[1]])$statistics
  Fit.tDraw_m6 <- summary(Chains.tDraw_m6[[1]])$statistics
  Convergence.gDraw_m6 <- gelman.diag(Chains.gDraw_m6)$mpsrf
  Convergence.tDraw_m6 <- gelman.diag(Chains.tDraw_m6)$mpsrf
  mESS.gDraw_m6 <- lapply(1:nChain, function(chain) multiESS(Chains.gDraw_m6[[chain]]))
  mESS.tDraw_m6 <- lapply(1:nChain, function(chain) multiESS(Chains.tDraw_m6[[chain]]))
  mcESS.gDraw_m6 <- lapply(1:nChain, function(chain) ess(Chains.gDraw_m6[[chain]]))
  mcESS.tDraw_m6 <- lapply(1:nChain, function(chain) ess(Chains.tDraw_m6[[chain]]))
  mcse.gDraw_m6 <- lapply(1:nChain, function(chain) mcse.mat(Chains.gDraw_m6[[chain]]))
  mcse.tDraw_m6 <- lapply(1:nChain, function(chain) mcse.mat(Chains.tDraw_m6[[chain]]))
  
}
ESS.out_m6 <- c(if(length(Fixed_m6)>0){list(bDraw = ESS.bDraw_m6, mc.bDraw = mcESS.bDraw_m6, mESS.bDraw = mESS.bDraw_m6)}, 
                if(length(Random_m6)>0){list(gDraw = ESS.gDraw_m6, mc.gDraw = mcESS.gDraw_m6, mESS.gDraw = mESS.gDraw_m6,
                                             tDraw = ESS.tDraw_m6, mc.tDraw = mcESS.tDraw_m6, mESS.tDraw = mESS.tDraw_m6)})
AC.out_m6 <- c(if(length(Fixed_m6)>0){list(bDraw = AC.bDraw_m6)}, if(length(Random_m6)>0){list(gDraw = AC.gDraw_m6, tDraw = AC.tDraw_m6)})
Fit.out_m6 <- c(if(length(Fixed_m6)>0){list(bDraw = Fit.bDraw_m6)}, if(length(Random_m6)>0){list(gdraw = Fit.gDraw_m6, tDraw = Fit.tDraw_m6)})
Convergence.out_m6 <- c(if(length(Fixed_m6)>0){list(bDraw = Convergence.bDraw_m6)}, if(length(Random_m6)>0){list(gDraw = Convergence.gDraw_m6, tDraw = Convergence.tDraw_m6)})
Se.out_m6 <- c(if(length(Fixed_m6)>0){list(
  mc.bDraw = mcse.bDraw_m6)}, 
  if(length(Random_m6)>0){list( 
    mc.gDraw = mcse.gDraw_m6,
    mc.tDraw = mcse.tDraw_m6)})

Diags_m6 <- list(ESS = ESS.out_m6, AC = AC.out_m6, Fit = Fit.out_m6, Convergence = Convergence.out_m6, Se.out = Se.out_m6)

# Restructure and thin regression parameters
EstRC_m6 <- lapply(1:(nIt/nThin), function(i){
  lapply(1:J, function(j) {rbind(if(length(Fixed_m6) > 0){Pars_m6[["Pars"]][[1]][["bDrawPG"]][[i]]}, 
        if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["gjDrawPG"]][[i]][,,j]})})})

if(any(grepl("x", c(Fixed_m6,Random_m6)))){iTypes_m6 <- 1:nrow(Types)
}else{
  iTypes_m6 <- which(Types[,"Methods"] %in% c("Value", "MvB"))
}

Tau_m6 <- lapply(1:(nIt/nThin), function(i){
        if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["tauDrawPG"]][[i]]}})

Average_Tau_m6 <- lapply(1:3, function(i){
    Reduce("+", lapply(Tau_m6, "[[", i)) / length(Tau_m6)
  })

bDraw_m6 <- lapply(1:(nIt/nThin), function(i){
  if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["bDrawPG"]][[i]]}})

Average_bDraw_m6 <- Reduce("+", bDraw_m6) / length(bDraw_m6)


gDraw_m6 <- lapply(1:(nIt/nThin), function(i){
  if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["gDrawPG"]][[i]]}})

Average_gDraw_m6 <- 
  Reduce("+", gDraw_m6) / length(gDraw_m6)


# Transform to success probabilities
Theta_m6 <- Transform2Theta(EstPars = EstRC_m6, X = xDataC_m6, Y = yDataC, nIt = nIt/nThin,
                             Types = Types[iTypes_m6,,drop=FALSE], 
                            Range = RangesApp[["Continuous"]], Value = ValuesApp[["Continuous"]],
                            Fixed = Fixed_m6, Random = Random_m6, PriorAlpha = rep(0.01,4))

# Compute treatment differences
Decisions_m6 <- lapply(1:length(iTypes_m6), function(i){
  #if(Types[iTypes_m6[i],"Methods"] != "MvB"){
  DeltaData <- Map("-", Theta_m6[["mTheta.E"]][[i]], Theta_m6[["mTheta.C"]][[i]])
  #}else{
   # DeltaData <- Map("-", Theta_m6[["mTheta.E"]][[i]], Theta_m6[["mTheta.C"]][[i]])
  #}
  EvaluateData(Data = DeltaData, 
               Weights = Weights, Alpha = Alpha, Alternative = c("greater.than", "two.sided"), 
               Rule = c("Any", "All", "Compensatory"), Truth = NULL, Types = Types[iTypes_m6[i],], nIt=nIt/nThin)
})

save(xDataC_m6, Pars_m6, Random_m6, Fixed_m6, Diags_m6, Theta_m6, Decisions_m6, 
     file = "Application/Workspaces/Application_m6_p10_H.RData")

#### 2. Single-level (non-hierarchical) logistic regression ####
Fixed_NH_m6 <- c("Intercept", "Trt", "x", "Trt_x")  # Names of fixed predictors
Random_NH_m6 <- c(NULL)                             # Names of random predictors
pF_NH_m6 <- length(Fixed_NH_m6)# Number of fixed predictors
pR_NH_m6 <- length(Random_NH_m6)# Number of random predictors
bMu0_NH_m6 <- rep(0, pF_NH_m6)# Prior mean of fixed predictors
gMu0_NH_m6 <- rep(0, pR_NH_m6)# Prior mean of random predictors

bSigma0_NH_m6 <- diag(1e-1, pF_NH_m6) # Prior covariance matrix of fixed predictors
gSigma0_NH_m6 <- diag(1e-1, pR_NH_m6) # Prior covariance matrix of random predictors

nu0_NH_m6 <- pR_NH_m6           # Prior degrees of freedom covariance matrix random effects
Tau0_NH_m6 <- diag(1e-1, pR_NH_m6)# Prior covariance matrix of covariance matrix random effects

# Covariate data
xDataC_NH_m6 <- lapply(which(nJ > 0), function(j){
  y <- as.matrix(cbind(1, DataC[[j]][,c(xVars)]))
  colnames(y) <- c("Intercept", "Trt", "x", "Trt_x")
  z <- y[,c(Fixed_NH_m6, Random_NH_m6),drop=FALSE]
  return(z)})
set.seed(2909)

# Sample parameters
Pars_NH_m6 <- EstimateParameters(X = xDataC_NH_m6, Y = yDataC, Fixed = Fixed_NH_m6, Random = Random_NH_m6, 
                                 nBurn = nBurn, nIt = nIt, Start = StartVals, 
                                 bMu0 = bMu0_NH_m6, bSigma0 = bSigma0_NH_m6, 
                                 gMu0 = gMu0_NH_m6, gSigma0 = gSigma0_NH_m6,
                                 nu0 = nu0_NH_m6, Tau0 = Tau0_NH_m6,nChain = nChain, ReturnThinned = TRUE, nThin = nThin)


# Diagnose MCMC samples - fixed parameters (if applicable)
if(length(Fixed_NH_m6) > 0){
  Chains.bDraw_NH_m6 <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars_NH_m6[["Pars"]][[chain]][["bDrawPG"]]), nrow = nIt/nThin, ncol = length(Fixed_NH_m6) * Q, byrow = TRUE)[,-((length(Fixed_NH_m6) * (Q-1)) + 1:(length(Fixed_NH_m6)))]))) #[,-((Q-1)*length(Fixed_NH_m6)+1:length(Fixed_NH_m6))])) 
  ESS.bDraw_NH_m6 <- effectiveSize(Chains.bDraw_NH_m6[[1]])
  AC.bDraw_NH_m6 <- autocorr.diag(Chains.bDraw_NH_m6[[1]], lags = 1:100)
  Fit.bDraw_NH_m6 <- summary(Chains.bDraw_NH_m6)$statistics
  Convergence.bDraw_NH_m6 <- gelman.diag(Chains.bDraw_NH_m6)$mpsrf
  mESS.bDraw_NH_m6 <- lapply(1:nChain, function(chain) multiESS(Chains.bDraw_NH_m6[[chain]]))
  mcESS.bDraw_NH_m6 <- lapply(1:nChain, function(chain) ess(Chains.bDraw_NH_m6[[chain]]))
  mcse.bDraw_NH_m6 <- lapply(1:nChain, function(chain) mcse.mat(Chains.bDraw_NH_m6[[chain]]))
}

# Diagnose MCMC samples - random parameters (if applicable)
if(length(Random_NH_m6) > 0){
  Chains.tDraw_NH_m6 <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(lapply(1:(nIt/nThin), function(i) lapply(1:(Q-1), function(q) {x <- Pars_NH_m6[["Pars"]][[chain]][["tauDrawPG"]][[i]][[q]]; x[lower.tri(x, diag = TRUE)]}))), nrow = nIt / nThin, byrow=TRUE))))
  Chains.gDraw_NH_m6 <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars_NH_m6[["Pars"]][[chain]][["gDrawPG"]]), nrow = nIt/nThin, ncol = length(Random_NH_m6) * Q, byrow = TRUE)[,-((length(Random_NH_m6) * (Q-1)) + 1:(length(Random_NH_m6)))])))
  ESS.gDraw_NH_m6 <- effectiveSize(Chains.gDraw_NH_m6[[1]])
  ESS.tDraw_NH_m6 <- effectiveSize(Chains.tDraw_NH_m6[[1]])
  AC.gDraw_NH_m6 <- autocorr.diag(Chains.gDraw_NH_m6[[1]], lags = c(1:10))
  AC.tDraw_NH_m6 <- autocorr.diag(Chains.tDraw_NH_m6[[1]], lags = c(1:10))
  Fit.gDraw_NH_m6 <- summary(Chains.gDraw_NH_m6[[1]])$statistics
  Fit.tDraw_NH_m6 <- summary(Chains.tDraw_NH_m6[[1]])$statistics
  Convergence.gDraw_NH_m6 <- gelman.diag(Chains.gDraw_NH_m6)$mpsrf
  Convergence.tDraw_NH_m6 <- gelman.diag(Chains.tDraw_NH_m6)$mpsrf
  mESS.gDraw_NH_m6 <- lapply(1:nChain, function(chain) multiESS(Chains.gDraw_NH_m6[[chain]]))
  mESS.tDraw_NH_m6 <- lapply(1:nChain, function(chain) multiESS(Chains.tDraw_NH_m6[[chain]]))
  mcESS.gDraw_NH_m6 <- lapply(1:nChain, function(chain) ess(Chains.gDraw_NH_m6[[chain]]))
  mcESS.tDraw_NH_m6 <- lapply(1:nChain, function(chain) ess(Chains.tDraw_NH_m6[[chain]]))
  mcse.gDraw_NH_m6 <- lapply(1:nChain, function(chain) mcse.mat(Chains.gDraw_NH_m6[[chain]]))
  mcse.tDraw_NH_m6 <- lapply(1:nChain, function(chain) mcse.mat(Chains.tDraw_NH_m6[[chain]]))
  
}
ESS.out_NH_m6 <- c(if(length(Fixed_NH_m6)>0){list(bDraw = ESS.bDraw_NH_m6, mc.bDraw = mcESS.bDraw_NH_m6, mESS.bDraw = mESS.bDraw_NH_m6)}, 
                   if(length(Random_NH_m6)>0){list(gDraw = ESS.gDraw_NH_m6, mc.gDraw = mcESS.gDraw_NH_m6, mESS.gDraw = mESS.gDraw_NH_m6,
                                                   tDraw = ESS.tDraw_NH_m6, mc.tDraw = mcESS.tDraw_NH_m6, mESS.tDraw = mESS.tDraw_NH_m6)})
AC.out_NH_m6 <- c(if(length(Fixed_NH_m6)>0){list(bDraw = AC.bDraw_NH_m6)}, if(length(Random_NH_m6)>0){list(gDraw = AC.gDraw_NH_m6, tDraw = AC.tDraw_NH_m6)})
Fit.out_NH_m6 <- c(if(length(Fixed_NH_m6)>0){list(bDraw = Fit.bDraw_NH_m6)}, if(length(Random_NH_m6)>0){list(gdraw = Fit.gDraw_NH_m6, tDraw = Fit.tDraw_NH_m6)})
Convergence.out_NH_m6 <- c(if(length(Fixed_NH_m6)>0){list(bDraw = Convergence.bDraw_NH_m6)}, if(length(Random_NH_m6)>0){list(gDraw = Convergence.gDraw_NH_m6, tDraw = Convergence.tDraw_NH_m6)})
Se.out_NH_m6 <- c(if(length(Fixed_NH_m6)>0){list(
  mc.bDraw = mcse.bDraw_NH_m6)}, 
  if(length(Random_NH_m6)>0){list(
    mc.gDraw = mcse.gDraw_NH_m6,
     mc.tDraw = mcse.tDraw_NH_m6)})

Diags_NH_m6 <- list(ESS = ESS.out_NH_m6, AC = AC.out_NH_m6, Fit = Fit.out_NH_m6, Convergence = Convergence.out_NH_m6, Se.out = Se.out_NH_m6)

# Restructure and thin regression parameters
EstRC_NH_m6 <- lapply(1:(nIt/nThin), function(i){ lapply(1:J, function(j) {
  rbind(if(length(Fixed_NH_m6) > 0){Pars_NH_m6[["Pars"]][[1]][["bDrawPG"]][[i]]}, 
        if(length(Random_NH_m6) > 0){Pars_NH_m6[["Pars"]][[1]][["gDrawPG"]][[i]]})})})

if(any(grepl("x", c(Fixed_NH_m6,Random_NH_m6)))){iTypes_NH_m6 <- 1:nrow(Types)
}else{
  iTypes_NH_m6 <- which(Types[,"Methods"] %in% c("Value", "MvB"))
}

# Compute success probabilities
Theta_NH_m6 <- Transform2Theta(EstPars = EstRC_NH_m6, X = xDataC_NH_m6, Y = yDataC, nIt = nIt/nThin,
                               Types = Types[iTypes_NH_m6,], 
                               Range = RangesApp[["Continuous"]], Value = ValuesApp[["Continuous"]],
                               Fixed = Fixed_NH_m6, Random = Random_NH_m6, PriorAlpha = rep(0.01, 4))

# Compute treatment differences
Decisions_NH_m6 <- lapply(1:length(iTypes_NH_m6), function(i){
  EvaluateData(Data = Map("-", Theta_NH_m6[["mTheta.E"]][[i]],Theta_NH_m6[["mTheta.C"]][[i]]), 
               Weights = Weights, Alpha = Alpha, Alternative = c("greater.than", "two.sided"), 
               Rule = c("Any", "All", "Compensatory"), Truth = NULL, Types = Types[iTypes_NH_m6[i],], nIt=nIt/nThin)
})


save(xDataC_NH_m6, Pars_NH_m6, Random_NH_m6, Fixed_NH_m6, Diags_NH_m6, Theta_NH_m6, Decisions_NH_m6, 
     file = "Application/Workspaces/Application_m6_p10_NH.RData")

