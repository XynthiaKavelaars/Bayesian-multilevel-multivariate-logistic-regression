#### Data Generating Mechanisms ####
Rho <- Rho_Lo <- Rho_Hi <- -0.20

#### 1. Delta = 0 - heterogeneous delta ####
# True success probabilities for low (Lo) and high (Hi) subpopulations of experimental (E) and control (C) treatments.
Theta1E_Lo <- c(0.625,0.575)
Theta1C_Lo <- c(0.375,0.425)
Theta1E_Hi <- c(0.375,0.425)
Theta1C_Hi <- c(0.625,0.575)

# Compute associated true regression parameters
TrueBeta1.C <- FindTrueBeta(Theta1C_Lo, Theta1E_Lo, Theta1C_Hi, Theta1E_Hi, Rho, Rho, xLo = Values[["Continuous"]][["Intra_Lo"]], xHi = Values[["Continuous"]][["Intra_Hi"]], Fixed = Fixed, Random = Random)


#### TrueBeta: list of true regression coefficients for all DGMs ####
TrueBeta <- list( 
  TrueBeta1.C = TrueBeta1.C)


# Number of data generating mechisms
nDgm <- length(TrueBeta)

# True transformed parameters for values (sTrueVal) and ranges (eTrueVal) of trial populations (ATE) and subpopulations (CTE)
Truth <- lapply(1:nDgm, function(dgm) vector("list", length(Sigma)))
for(dgm in 1:nDgm){
  for(s in 1:length(Sigma)){
    TruePars <- list(b = TrueBeta[[dgm]][Fixed,,drop=FALSE], g = TrueBeta[[dgm]][Random,,drop=FALSE], Sigma = Sigma[s])
    TruePars_NH <- list(b = TrueBeta[[dgm]][Fixed_NH,,drop=FALSE], g = TrueBeta[[dgm]][Random_NH,,drop=FALSE], Sigma = Sigma[s])
    Tau0 <- diag(TruePars$Sigma, length(Random))
    sTrueVal.ATE <- sTrueValues(TrueBeta = TrueBeta[[dgm]], ValueX = Values[["Continuous"]][["Trial"]], Weights = Weights, Fixed = Fixed, Random = Random)
    sTrueVal.CTE <- sTrueValues(TrueBeta = TrueBeta[[dgm]], ValueX = Values[["Continuous"]][["Intra_Lo"]], Weights = Weights, Fixed = Fixed, Random = Random)
    eTrueVal.ATE <- eTrueValues(TrueBeta = TrueBeta[[dgm]], 
                                pXD = pX, MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[["Continuous"]][["Trial"]], Continuous = TRUE, Weights = Weights,
                                Fixed = Fixed, Random = Random)
    eTrueVal.CTE <- eTrueValues(TrueBeta = TrueBeta[[dgm]], 
                                pXD = pX, MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[["Continuous"]][["Intra_Lo"]], Continuous = TRUE, Weights = Weights,
                                Fixed = Fixed, Random = Random)
    
    Truth[[dgm]][[s]] <- list(sATE = sTrueVal.ATE, sCTE = sTrueVal.CTE,
                              eATE = eTrueVal.ATE, eCTE = eTrueVal.CTE)
  }}
