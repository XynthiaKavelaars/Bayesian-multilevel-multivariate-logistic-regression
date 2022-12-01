#### 1. Run simulation ####
# Select conditions - multilevel model
if(any(grepl("x", c(Fixed,Random)))){iTypes <- which(!(Types[,"Methods"] %in% c("MvB")))
}else{
  iTypes <- which(Types[,"Methods"] %in% c("Value", "MvB"))
}

# Select conditions - single-level model
if(any(grepl("x", c(Fixed_NH,Random_NH)))){iTypes_NH <- which(!(Types[,"Methods"] %in% c("MvB")))
}else{
  iTypes_NH <- which(Types[,"Methods"] %in% c("Value"))
}

# Select condition - multivariate Bernoulli model
iTypes_MB <- which(Types[,"Methods"] %in% c("MvB"))
#### 1.1 Multilevel (BMMLR) & Single-level (BMLR) models ####

for(dgm in 1:length(TrueBeta)){
  for(s in 1:length(Sigma)){
    
    # True parameters 
    TruePars <- list(b = TrueBeta[[dgm]][Fixed,,drop=FALSE], g = TrueBeta[[dgm]][Random,,drop=FALSE], Sigma = Sigma[s])
    TruePars_NH <- list(b = TrueBeta[[dgm]][Fixed_NH,,drop=FALSE], g = TrueBeta[[dgm]][Random_NH,,drop=FALSE], Sigma = Sigma[s])
    
    # True values
    Tau0 <- diag(TruePars$Sigma, length(Random))
    Truths <- lapply(1:nrow(Types), function(i){
       if(Types[i,"Methods"] == "Value" | !any(grepl("x", c(Fixed,Random)))){
         y <- sTrueValues(TrueBeta = TrueBeta[[dgm]], ValueX = Values[["Continuous"]][[Types[i,"Populations"]]], Weights = Weights, Fixed = Fixed, Random = Random)
       }else if(Types[i,"Methods"] %in% c("Analytical", "Empirical", "MvB") & any(grepl("x", c(Fixed,Random)))){
         y <- eTrueValues(TrueBeta = TrueBeta[[dgm]], 
                          pXD = pX, MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[["Continuous"]][[Types[i,"Populations"]]], Continuous = TRUE, Weights = Weights,
                          Fixed = Fixed, Random = Random)
      }else{y <- NULL}
      return(y)})

    
    for(ss in 1:(nrow(SampleSizes))){
      set.seed(Seeds[dgm,s,ss])
     foreach(sim = 1:nSim, .packages = packages, .verbose = TRUE)%dopar%{
       set.seed(Seeds[dgm,s,ss]*sim)
       
        # Generate data
        Data <- GenerateData(TruePars, Q = Q, pX = pX, Mux = MuX, SigmaX = SigmaX, J = SampleSizes[ss,"nCluster"], n = SampleSizes[ss,"nObs"], Continuous = Continuous, Fixed = Fixed, Random = Random)
        
        # Sample parameters - multilevel (BMMLR)
        set.seed(sim^2)
        Pars <- EstimateParameters(X = Data[["X"]], Y = Data[["Y"]], Fixed = Fixed, Random = Random,
                                   nBurn = nBurn, nIt = nIt, Start = StartVals, 
                                   bMu0 = bMu0, bSigma0 = bSigma0, 
                                   gMu0 = gMu0, gSigma0 = gSigma0,
                                   nu0 = nu0, Tau0 = Tau0, nChain = nChain)
        
        # Compute bias in regression parameters - BMMLR
        BiasPars <- EstimateBiasPars(Pars = Pars[["Pars"]][[1]], TruePar = TruePars, Fixed = Fixed, Random = Random, nPars = nIt)
        
        Pars.Fixed <- Pars[["Pars"]][[1]]["bDrawPG"]
        Pars.Random <- Pars[["Pars"]][[1]]["gDrawPG"]
        
        # Restructure and thin regression parameters - BMMLR
        EstRC <- lapply(seq(1,nIt, nThin), function(i){
            rbind(if(length(Fixed) > 0){Pars.Fixed[["bDrawPG"]][[i]]}, 
                  if(length(Random) > 0){Pars.Random[["gDrawPG"]][[i]]})})
        
        # Transform parameters - BMMLR
        Theta <- Transform2Theta(EstPars = EstRC, X = Data[["X"]], Y = Data[["Y"]], nIt = nIt,
                                 Types = Types[iTypes,,drop=FALSE], 
                                 Range = Ranges[["Continuous"]], Value = Values[["Continuous"]], 
                                 Fixed = Fixed, Random = Random, PriorAlpha = PriorAlpha)
        
        # Make decision - BMMLR
        Decisions <- lapply(1:length(iTypes), function(i){
          EvaluateData(Data = Map("-", Theta[["mTheta.E"]][[i]],  
                                  Theta[["mTheta.C"]][[i]]),
                       Weights = Weights, Alpha = Alpha, Alternative = c("greater.than"), 
                       Rule = c("Any", "All", "Compensatory"), Truth = Truths[[i]], Types = Types[iTypes[i],])
        })
        
        # Diagnose MCMC samples - fixed parameters (if applicable) - BMMLR
        if(length(Fixed) > 0){
          Chains.bDraw <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars[["Pars"]][[chain]][["bDrawPG"]][seq(1,nIt,nThin)]), nrow = nIt/nThin, ncol = length(Fixed) * Q, byrow = TRUE)[,-((length(Fixed) * (Q-1)) + 1:(length(Fixed)))]))) #[,-((Q-1)*length(Fixed)+1:length(Fixed))])) 
          ESS.bDraw <- effectiveSize(Chains.bDraw)
          AC.bDraw <- autocorr(Chains.bDraw, lags = 1)[[1]][1,,]
          Fit.bDraw <- summary(Chains.bDraw)$statistics
          Convergence.bDraw <- gelman.diag(Chains.bDraw)$mpsrf
          mESS.bDraw <- lapply(1:nChain, function(chain) multiESS(Chains.bDraw[[chain]]))
          mcESS.bDraw <- lapply(1:nChain, function(chain) ess(Chains.bDraw[[chain]]))
          mcse.bDraw <- lapply(1:nChain, function(chain) mcse.mat(Chains.bDraw[[chain]]))
                  }
        
        # Diagnose MCMC samples - random parameters (if applicable) - BMMLR
        if(length(Random) > 0){
          Chains.tDraw <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(lapply(seq(1,nIt,nThin), function(i) lapply(1:(Q-1), function(q) {x <- Pars[["Pars"]][[chain]][["tauDrawPG"]][[i]][[q]]; x[lower.tri(x, diag = TRUE)]}))), nrow = nIt / nThin, byrow=TRUE))))
          Chains.gDraw <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars[["Pars"]][[chain]][["gDrawPG"]][seq(1,nIt,nThin)]), nrow = nIt/nThin, ncol = length(Random) * Q, byrow = TRUE)[,-((length(Random) * (Q-1)) + 1:(length(Random)))])))
          ESS.gDraw <- effectiveSize(Chains.gDraw)
          ESS.tDraw <- effectiveSize(Chains.tDraw)
          AC.gDraw <- autocorr.diag(Chains.gDraw, lags = c(1:10))
          AC.tDraw <- autocorr.diag(Chains.tDraw, lags = c(1:10))
          Fit.gDraw <- summary(Chains.gDraw)$statistics
          Fit.tDraw <- summary(Chains.tDraw)$statistics
          Convergence.gDraw <- gelman.diag(Chains.gDraw)$mpsrf
          Convergence.tDraw <- gelman.diag(Chains.tDraw)$mpsrf
          mESS.gDraw <- lapply(1:nChain, function(chain) multiESS(Chains.gDraw[[chain]]))
          mESS.tDraw <- lapply(1:nChain, function(chain) multiESS(Chains.tDraw[[chain]]))
          mcESS.gDraw <- lapply(1:nChain, function(chain) ess(Chains.gDraw[[chain]]))
          mcESS.tDraw <- lapply(1:nChain, function(chain) ess(Chains.tDraw[[chain]]))
          mcse.gDraw <- lapply(1:nChain, function(chain) mcse.mat(Chains.gDraw[[chain]]))
          mcse.tDraw <- lapply(1:nChain, function(chain) mcse.mat(Chains.tDraw[[chain]]))
          
        }
        
        if(all(c(length(Random), length(Fixed)) > 0)){
         mESS <- lapply(1:nChain, function(chain) multiESS(cbind(Chains.bDraw[[chain]], Chains.gDraw[[chain]], Chains.tDraw[[chain]])))
        }else if(length(Random) > 0){
         mESS <-  lapply(1:nChain, function(chain) multiESS(cbind(Chains.gDraw[[chain]], Chains.tDraw[[chain]]))) 
        }else if(length(Fixed) > 0){
        mESS <- mESS.bDraw
           }
        
        ESS.out <- c(list(mESS = mESS),
                     if(length(Fixed)>0){list(bDraw = ESS.bDraw, mc.bDraw = mcESS.bDraw, mESS.bDraw = mESS.bDraw)}, 
                     if(length(Random)>0){list(gDraw = ESS.gDraw, mc.gDraw = mcESS.gDraw, mESS.gDraw = mESS.gDraw,
                                               tDraw = ESS.tDraw, mc.tDraw = mcESS.tDraw, mESS.tDraw = mESS.tDraw)})
        AC.out <- c(if(length(Fixed)>0){list(bDraw = AC.bDraw)}, if(length(Random)>0){list(gDraw = AC.gDraw, tDraw = AC.tDraw)})
        Fit.out <- c(if(length(Fixed)>0){list(bDraw = Fit.bDraw)}, if(length(Random)>0){list(gdraw = Fit.gDraw, tDraw = Fit.tDraw)})
        Convergence.out <- c(if(length(Fixed)>0){list(bDraw = Convergence.bDraw)}, if(length(Random)>0){list(gDraw = Convergence.gDraw, tDraw = Convergence.tDraw)})
        Se.out <- c(if(length(Fixed)>0){list(bDraw = mcse.bDraw)},
                    if(length(Random)>0){list(gDraw = mcse.gDraw, 
                                              tDraw = mcse.tDraw
)})
        
        Diags <- list(ESS = ESS.out, AC = AC.out, Fit = Fit.out, Convergence = Convergence.out, Se.out = Se.out)
       
        # Save data - BMMLR
        if(SaveRDS){
          saveRDS(Data, file = paste0(wd, "/Workspaces/Workspaces_H/Data/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Data", sim, ".RDS", collapse = ""))
          saveRDS(BiasPars, file = paste0(wd, "/Workspaces/Workspaces_H/Bias/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_BiasPars", sim, ".RDS", collapse = ""))
          saveRDS(Pars, file = paste0(wd, "/Workspaces/Workspaces_H/Parameters/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Pars", sim, ".RDS", collapse = ""))
         saveRDS(Decisions, file = paste0(wd, "/Workspaces/Workspaces_H/Decisions/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Decisions", sim, ".RDS", collapse = ""))
          saveRDS(Diags, file = paste0(wd, "/Workspaces/Workspaces_H/Diags/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Diags", sim, ".RDS", collapse = ""))
        }
        
        set.seed(sim^2)
        
        # Sample parameters - single-level (BMLR)
        Pars_NH <- EstimateParameters(X = Data[["X"]], Y = Data[["Y"]], Fixed = Fixed_NH, Random = Random_NH,
                                      nBurn = nBurn, nIt = nIt, Start = StartVals, 
                                      bMu0 = rep(b0, length(Fixed_NH)), bSigma0 = chol2inv(chol(diag(B0, length(Fixed_NH)))), 
                                      gMu0 = rep(g0, length(Random_NH)), gSigma0 = chol2inv(chol(diag(G0, length(Random_NH)))),
                                      nu0 = nu0, Tau0 = Tau0, nChain = nChain)
        
        # Compute bias in regression parameters - BMLR
        BiasPars_NH <- EstimateBiasPars(Pars = Pars_NH[["Pars"]][[1]], TruePar = TruePars_NH, Fixed = Fixed_NH, Random = Random_NH, nPars = nIt)
        
        Pars.Fixed_NH <- Pars_NH[["Pars"]][[1]]["bDrawPG"]
        Pars.Random_NH <- Pars_NH[["Pars"]][[1]]["gDrawPG"]
        
        # Restructure and thin regression parameters - BMLR
        EstRC_NH <- lapply(seq(1,nIt, nThin), function(i){
          rbind(if(length(Fixed_NH) > 0){Pars.Fixed_NH[["bDrawPG"]][[i]]}, 
                if(length(Random_NH) > 0){Pars.Random_NH[["gDrawPG"]][[i]]})})
        
        # Transform to success probabilities - BMLR
        Theta_NH <- Transform2Theta(EstPars = EstRC_NH, X = Data[["X"]], Y = Data[["Y"]], nIt = nIt,
                                    Types = Types[iTypes_NH,,drop=FALSE], 
                                    Range = Ranges[["Continuous"]], Value = Values[["Continuous"]], 
                                    Fixed = Fixed_NH, Random = Random_NH, PriorAlpha = PriorAlpha)
        
        # Make decision - BMLR
        Decisions_NH <- lapply(1:length(iTypes_NH), function(i){
          EvaluateData(Data = Map("-", Theta_NH[["mTheta.E"]][[i]],  
                                  Theta_NH[["mTheta.C"]][[i]]),
                       Weights = Weights, Alpha = Alpha, Alternative = c("greater.than"), 
                       Rule = c("Any", "All", "Compensatory"), Truth = Truths[[i]], Types = Types[iTypes_NH[i],])
        })
        
        # Diagnose MCMC samples - random parameters (if applicable) - BMLR
        if(length(Fixed_NH) > 0){
          Chains.bDraw_NH <- mcmc.list(lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Pars_NH[["Pars"]][[chain]][["bDrawPG"]][seq(1,nIt,nThin)]), nrow = nIt/nThin, ncol = length(Fixed_NH) * Q, byrow = TRUE)[,-((length(Fixed_NH) * (Q-1)) + 1:(length(Fixed_NH)))]))) #[,-((Q-1)*length(Fixed_NH)+1:length(Fixed_NH))])) 
          ESS.bDraw_NH <- effectiveSize(Chains.bDraw_NH)
          AC.bDraw_NH <- autocorr(Chains.bDraw_NH, lags = 1)[[1]][1,,]
          Fit.bDraw_NH <- summary(Chains.bDraw_NH)$statistics
          Convergence.bDraw_NH <- gelman.diag(Chains.bDraw_NH)$mpsrf
          mESS.bDraw_NH <- lapply(1:nChain, function(chain) multiESS(Chains.bDraw_NH[[chain]]))
          mcESS.bDraw_NH <- lapply(1:nChain, function(chain) ess(Chains.bDraw_NH[[chain]]))
          mcse.bDraw_NH <- lapply(1:nChain, function(chain) mcse.mat(Chains.bDraw_NH[[chain]]))
        }
        
        mESS_NH <- mESS.bDraw_NH
        ESS.out_NH <- list(mESS = mESS_NH)
        AC.out_NH <- list(bDraw = AC.bDraw_NH)
        Fit.out_NH <- list(bDraw = Fit.bDraw_NH)
        Convergence.out_NH <-list(bDraw = Convergence.bDraw_NH)
        Se.out_NH <- list(bDraw = mcse.bDraw_NH)
        
        Diags_NH <- list(ESS = ESS.out_NH, AC = AC.out_NH, Fit = Fit.out_NH, Convergence = Convergence.out_NH, Se.out = Se.out_NH)
        
        # Save data - BMLR
        if(SaveRDS){
          saveRDS(BiasPars_NH, file = paste0(wd, "/Workspaces/Workspaces_NH/Bias/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_BiasPars_NH", sim, ".RDS", collapse = ""))
          saveRDS(Pars_NH, file = paste0(wd, "/Workspaces/Workspaces_NH/Parameters/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Pars_NH", sim, ".RDS", collapse = ""))
          saveRDS(Decisions_NH, file = paste0(wd, "/Workspaces/Workspaces_NH/Decisions/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Decisions_NH", sim, ".RDS", collapse = ""))
          saveRDS(Diags_NH, file = paste0(wd, "/Workspaces/Workspaces_NH/Diags/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Diags_NH", sim, ".RDS", collapse = ""))
        }
       }
    }
    }
}

#### 1.2 Multivariate Bernoulli (BMB) model ####
iTypes_MB <- which(Types[,"Methods"] %in% c("MvB"))

for(dgm in 1:length(TrueBeta)){
  for(s in 1:length(Sigma)){
    
    # True parameters 
    TruePars <- list(b = TrueBeta[[dgm]][Fixed,,drop=FALSE], g = TrueBeta[[dgm]][Random,,drop=FALSE], Sigma = Sigma[s])
    TruePars_NH <- list(b = TrueBeta[[dgm]][Fixed_NH,,drop=FALSE], g = TrueBeta[[dgm]][Random_NH,,drop=FALSE], Sigma = Sigma[s])
    
    # True values
    Tau0 <- diag(TruePars$Sigma, length(Random))
    Truths <- lapply(1:nrow(Types), function(i){
      if(Types[i,"Methods"] == "Value" | !any(grepl("x", c(Fixed,Random)))){
        y <- sTrueValues(TrueBeta = TrueBeta[[dgm]], ValueX = Values[["Continuous"]][[Types[i,"Populations"]]], Weights = Weights, Fixed = Fixed, Random = Random)
      }else if(Types[i,"Methods"] %in% c("Analytical", "Empirical", "MvB") & any(grepl("x", c(Fixed,Random)))){
        y <- eTrueValues(TrueBeta = TrueBeta[[dgm]], 
                         pXD = pX, MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[["Continuous"]][[Types[i,"Populations"]]], Continuous = TRUE, Weights = Weights,
                         Fixed = Fixed, Random = Random)
      }else{y <- NULL}
      return(y)})
  
    for(ss in 1:(nrow(SampleSizes))){
      set.seed(Seeds[dgm,s,ss])
        foreach(sim = 1:nSim, .packages = packages, .verbose = TRUE)%dopar%{
        set.seed(Seeds[dgm,s,ss]*sim)
        # Generate data
        Data <- GenerateData(TruePars, Q = Q, pX = pX, Mux = MuX, SigmaX = SigmaX, J = SampleSizes[ss,"nCluster"], n = SampleSizes[ss,"nObs"], Continuous = Continuous, Fixed = Fixed, Random = Random)
       
        # Transform to success probabilities - BMB
        Theta <- Transform2Theta(EstPars = NULL, X = Data[["X"]], Y = Data[["Y"]], nIt = nIt,
                                 Types = Types[iTypes_MB,,drop=FALSE], 
                                 Range = Ranges[["Continuous"]], Value = Values[["Continuous"]], 
                                 Fixed = Fixed, Random = Random, PriorAlpha = PriorAlpha)
        
        # Make decisions - BMB
        Decisions <- lapply(1:length(iTypes_MB), function(i){
          EvaluateData(Data = Map("-", Theta[["mTheta.E"]][[i]],  
                                  Theta[["mTheta.C"]][[i]]),
                       Weights = Weights, Alpha = Alpha, Alternative = c("greater.than"), 
                       Rule = c("Any", "All", "Compensatory"), Truth = Truths[[i]], Types = Types[iTypes_MB[i],])  
              
        })
        
        # Save data - BMB
        saveRDS(Theta, file = paste0(wd, "/Workspaces/Workspaces_MB/Theta/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Theta", sim, ".RDS", collapse = ""))
        saveRDS(Decisions, file = paste0(wd, "/Workspaces/Workspaces_MB/Decisions/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Decisions", sim, ".RDS", collapse = ""))
      }}
  }}



