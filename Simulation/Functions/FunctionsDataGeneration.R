#### 1. FindTrueBeta: Transform phi to regression coefficients ####
# Function to compute regression coefficients from success probabilities for a model with treatment indicator, covariate, and interaction between treatment and covariate.
# Input:
# - TrueThetaC_Lo: Vector of two success probabilities of treatment C for a population with low value of x
# - TrueThetaE_Lo: Vector of two success probabilities of treatment E for a population with low value of x
# - TrueThetaC_Hi: Vector of two success probabilities of treatment C for a population with high value of x
# - TrueThetaE_Hi: Vector of two success probabilities of treatment E for a population with high value of x
# - TrueRho_Lo: Scalar correlation between outcomes for a population with low value of x
# - TrueRho_Hi: Scalar correlation between outcomes for a population with high value of x
# - xLo: Scalar: low value of x
# - xHi: Scalar: high value of x
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.

# Output: 
# - TrueBeta: (P x Q) matrix with regression coefficients per joint response category. Rows ordered as (Fixed, Random). Currently supported for K=2/Q=4 only.
FindTrueBeta <- function(TrueThetaC_Lo, TrueThetaE_Lo, TrueThetaC_Hi, TrueThetaE_Hi, TrueRho_Lo, TrueRho_Hi, xLo, xHi, Fixed, Random){
  K <- length(TrueThetaC_Lo)
  Q <- 2^K
  if(Q == 4){
    PhiC_Lo <- Theta2Phi(TrueThetaC_Lo,TrueRho_Lo)
    PhiC_Hi <- Theta2Phi(TrueThetaC_Hi,TrueRho_Hi)
    PhiE_Lo <- Theta2Phi(TrueThetaE_Lo,TrueRho_Lo)
    PhiE_Hi <- Theta2Phi(TrueThetaE_Hi,TrueRho_Hi)
  }else if(Q == 2){
    PhiC_Lo <- c(TrueThetaC_Lo, 1 - TrueThetaC_Lo)
    PhiC_Hi <- c(TrueThetaC_Hi, 1 - TrueThetaC_Hi)
    PhiE_Lo <- c(TrueThetaE_Lo, 1 - TrueThetaE_Lo)
    PhiE_Hi <- c(TrueThetaE_Hi, 1 - TrueThetaE_Hi)
  }
  
  PsiC_Lo <- log(PhiC_Lo) - log(PhiC_Lo[Q])
  PsiC_Hi <- log(PhiC_Hi) - log(PhiC_Hi[Q])
  PsiE_Lo <- log(PhiE_Lo) - log(PhiE_Lo[Q])
  PsiE_Hi <- log(PhiE_Hi) - log(PhiE_Hi[Q])
  
  Intercept <- (xHi * PsiC_Lo - xLo * PsiC_Hi) / (xHi - xLo)
  Trt <- (xLo * (PsiC_Hi - PsiE_Hi) + xHi * (PsiE_Lo - PsiC_Lo)) / (xHi - xLo)
  x <- (PsiC_Hi - PsiC_Lo) / (xHi - xLo)
  Trt_x <- ((PsiE_Hi - PsiC_Hi) - (PsiE_Lo - PsiC_Lo)) / (xHi - xLo)
  
  TrueBeta <- do.call(rbind, c(if(length(Fixed>0)){mget(Fixed)}, if(length(Random>0)){mget(Random)}))
  return(TrueBeta)
}

#### 2. eTrueValues: True parameters over range of covariate ####
# Function to perform numerical integration with a single covariate
# Input: 
# - TrueBeta: (P x Q) matrix with true regression coefficients. Rows ordered as (Fixed, Random). Currently supported for Q=2 and Q=4 (i.e. K=1 and K=2).
# - pXD: Scalar (0 =< pXD =< 1) with proportion of discrete x
# - MuX: Scalar with mean of continuous x
# - SigmaX: Scalar with standard deviation of continuous x
# - RangeX: Vector with lower and upper bound of integration
# - Continuous: Logical. TRUE if x is continuous; FALSE if x is discrete.
# - Weights: Vector of two positively valued weights of weighted linear combination. Weights sum to 1.
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.

# Output: 
# - ThetaE: Vector of two success probabilities of treatment E
# - ThetaC: Vector of two success probabilities of treatment C
# - RhoE: Scalar correlation of treatment E
# - RhoC: Scalar correlation of treatment C
# - PhiE: Vector of four joint response probabilities of treatment E
# - PhiC: Vector of four joint response probabilities of treatment C
# - Delta: Vector of two treatment differences (E-C)
# - DeltaW: Scalar weighted treatment difference (E-C) (if K= 2)
eTrueValues <- function(TrueBeta, pXD = NULL, MuX = NULL, SigmaX = NULL, RangeX, Continuous, Weights, Fixed, Random){
  Q <- dim(TrueBeta)[2] 
  PhiE <- PhiC <- rep(NA, Q)
  for(q in 1:Q){
    PhiE[q] <- integrate(Integrand, lower = RangeX[1], upper = RangeX[2], 
                         EstRC = TrueBeta, Trt = 1, q = q,
                         MuX = MuX, SigmaX = SigmaX, RangeX = RangeX,
                         Fixed = Fixed, Random = Random)$value
    PhiC[q] <- integrate(Integrand, lower = RangeX[1], upper = RangeX[2], 
                         EstRC = TrueBeta, Trt = 0, q = q,
                         MuX = MuX, SigmaX = SigmaX, RangeX = RangeX,
                         Fixed = Fixed, Random = Random)$value
  }
  
  if(Q == 4){
    ThetaE <- c(sum(PhiE[c(1,2)]), sum(PhiE[c(1,3)]))
    ThetaC <- c(sum(PhiC[c(1,2)]), sum(PhiC[c(1,3)]))
    RhoE <- ComputeRho(PhiE)
    RhoC <- ComputeRho(PhiC)
  }else if(Q ==2){
    ThetaE <- PhiE
    ThetaC <- PhiC
  }
  
  Delta <- ThetaE - ThetaC
  
  if(Q == 4){
    DeltaW <- Delta %*% Weights
    Out <- list(ThetaE = ThetaE, ThetaC = ThetaC, RhoE = RhoE, RhoC = RhoC, PhiE = PhiE, PhiC = PhiC, Delta = Delta, DeltaW = DeltaW)
  }else if(Q == 2){
    Out <- list(ThetaE = ThetaE, ThetaC = ThetaC, PhiE = PhiE, PhiC = PhiC, Delta = Delta)
  }
  return(Out)
}

#### 3. sTrueValues: True parameters for value of covariate ####  
# Function to compute true parameters for a value of a single covariate
# Input: 
# - TrueBeta: P x Q matrix with true regression coefficients. Rows ordered as (Fixed, Random). 
# - ValueX: Scalar value of x
# - Weights: Vector of two positively valued weights of weighted linear combination. Weights sum to 1.
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.

# Output: 
# - ThetaE: Vector of two success probabilities of treatment E
# - ThetaC: Vector of two success probabilities of treatment C
# - RhoE: Scalar correlation of treatment E
# - RhoC: Scalar correlation of treatment C
# - PhiE: Vector of four joint response probabilities of treatment E
# - PhiC: Vector of four joint response probabilities of treatment C
# - Delta: Vector of two treatment differences (E-C)
# - DeltaW: Scalar weighted treatment difference (E-C)
sTrueValues <- function(TrueBeta, ValueX, Weights, Fixed, Random){
  Intercept <- 1
  x <- ValueX
  
  Trt <- 1
  Trt_x <- Trt * ValueX
  xValsE <- c(Fixed, Random)
  XE <- t(do.call(rbind, mget(xValsE)))
  
  Trt <- 0
  Trt_x <- Trt * ValueX
  xValsC <- c(Fixed, Random)
  XC <- t(do.call(rbind, mget(xValsC)))
  
  PsiE <- XE %*% TrueBeta
  PsiC <- XC %*% TrueBeta
  
  PhiE <- exp(PsiE) / rowSums(exp(PsiE))
  PhiC <- exp(PsiC) / rowSums(exp(PsiC))
  
  if(ncol(PhiE) == 4){
    ThetaE <- c(rowSums(PhiE[,c(1,2), drop=FALSE]), rowSums(PhiE[,c(1,3), drop=FALSE]))
    ThetaC <- c(rowSums(PhiC[,c(1,2), drop=FALSE]), rowSums(PhiC[,c(1,3), drop=FALSE]))
    
    RhoE <- ComputeRho(PhiE)
    RhoC <- ComputeRho(PhiC)
  }else if(ncol(PhiE) == 2){
    ThetaE <- PhiE
    ThetaC <- PhiC
  }
  
  Delta <- ThetaE - ThetaC
  
  
  if(ncol(PhiE) == 4){
    DeltaW <- Delta %*% Weights
    Out <- list(ThetaE = ThetaE, ThetaC = ThetaC, RhoE = RhoE, RhoC = RhoC, PhiE = PhiE, PhiC = PhiC, Delta = Delta, DeltaW = DeltaW)
  }else if(ncol(PhiE) == 2){
    Out <- list(ThetaE = ThetaE, ThetaC = ThetaC,PhiE = PhiE, PhiC = PhiC, Delta = Delta)
  }
  
  return(Out)
}


