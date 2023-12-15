#### 1. EstimateParametersPG: Wrapper for MCMC procedure to estimate regression coefficients - Polya-Gamma ####
# Function to estimate regression coefficients using Polya-Gamma gibbs sampling. 
# Input:
# - X: n x P design matrix
# - Y: n x Q matrix of multinomial response data
# - nBurn: Scalar. Number of burnin iterations
# - nIt: Scalar. Number of iterations
# - Start: Vector of two starting values, each used for a different chain
# - bMu0: Scalar. Prior mean.
# - bSigma0: Scalar. Prior variance of regression coefficients

# Output: 
# A list with:
# - ConvergencePG: Multivariate Gelman-Rubin statistic to asses convergence
# - bDrawPG: List of length nIt with P x Q matrices with estimated regression coefficients
EstimateParametersPG <- function(X, Y, nBurn, nIt, Start, bMu0, bSigma0){
  tryCatch({
    P <- ncol(X)
    Q <- ncol(Y)
    
    Chain1PG <- SampleBetaPG(X=X, Y=Y, nBurn=nBurn, nIt=nIt, Start=Start[1], bMu0=bMu0, bSigma0=bSigma0)
    Chain2PG <- SampleBetaPG(X=X, Y=Y, nBurn=nBurn, nIt=nIt, Start=Start[2], bMu0=bMu0, bSigma0=bSigma0)
    
    ChainsPG <- mcmc.list(mcmc(do.call(rbind, lapply(Chain1PG[["bDrawPG"]], function(i) as.vector(i[,-seq(Q, Q * P, Q)])))),
                          mcmc(do.call(rbind, lapply(Chain2PG[["bDrawPG"]], function(i) as.vector(i[,-seq(Q, Q * P, Q)])))))
    
    ConvergencePG <- gelman.diag(ChainsPG)$mpsrf
    bDrawPG <- Chain1PG[["bDrawPG"]]  
    res <- list(ConvergencePG = ConvergencePG, bDrawPG = bDrawPG, ChainsPG = ChainsPG)
  },
  warning = function(warn) {
    print(paste("Warning: ", warn))
    res <- list(NULL, NULL)
  },
  error = function(err) {
    print(paste("Error: ", err))
    res <- list(NULL,NULL)
  },
  finally = function(f) {
    return(res)
  })
  
}

#### 2. SampleBetaPG: Sample regression coefficients via Polya-Gamma logistic regression ####
# Function to sample regression coefficients via Polya-Gamma multinomial logistic regression
# Input:
# - X: n x P design matrix
# - Y: n x Q matrix of multinomial response data
# - nBurn: Scalar. Number of burnin iterations
# - nIt: Scalar. Number of iterations
# - Start: Vector of two starting values, each used for a different chain
# - bMu0: Scalar. Prior mean.
# - bSigma0: Scalar. Prior variance of regression coefficients

# Output:
# A list with:
# - bDrawPG: List of length nIt with P x Q matrices with estimated regression coefficients

SampleBetaPG <- function(X, Y, nBurn, nIt, Start, bMu0, bSigma0){
  
  P <- ncol(X);
  n <- nrow(Y);
  Q <- ncol(Y);
  
  bDrawPG <- list("vector", nIt);
  betaPG <- array(Start, dim = c(P,Q));
  betaPG[,Q] <- rep(0,P);
  kappa <- Y - 1/2 # Note: As opposed to the paper: kappa here is not divided by omega to improve programming efficiency.
  
  b0 <- rep(bMu0,P)
  B0 <- diag(bSigma0,P,P) # Prior precision matrix
  P0 <- B0 %*% b0;
  
  # Start Gibbs sampler
  for(i in 2:(nBurn + nIt)){
    for(q in 1:(Q-1)){
      # Compute linear predictor
      C <- log(rowSums(exp(X %*% betaPG[,-q])));
      eta <- X %*% betaPG[,q] - C;
      
      # Draw Polya-Gamma variable
      omegaDraw <- pgdraw(1,eta)
      
      # Draw regression coefficients
      bSigma <- chol2inv(chol(t(X) %*% (X * omegaDraw) + B0))
      bMu <- bSigma %*% (t(X) %*% (kappa[,q] + omegaDraw * C) + P0)
      betaPG[,q] <- bMu + t(chol(bSigma)) %*% rnorm(P) 
    } 
    if(i > nBurn){bDrawPG[[i-nBurn]] <- betaPG}
  }
  
  return(list(bDrawPG = bDrawPG))
}

#### 3. Sample population ####
# Input: 
# - X: n x P design matrix with covariate data
# - Method: "Empirical" for empirical marginalization, "Value" for vector of fixed values, "MvB" for unconditional multivariate Bernoulli (reference approach)
# - Values: Scalar. Value of x representing subpopulation.
# - Range: Vector of lower and upper bound of covariate that represents subpopulation.

# Output: 
# - xE: nE x P design matrix with covariate data for treatment E
# - xC: nC x P design matrix with covariate data for treatment C
SamplePopulation <- function(X, Method, Values, Ranges){
  
  if("Value" %in% Method){
      xE <- cbind(1, 1, Values, Values)
      xC <- cbind(1, 0, Values, 0)
    
  } else if(any(c("MvB", "Empirical") %in% Method)){
        xE <- X[X[,2] == 1 & X[,3] >= min(Ranges) & X[,3] <= max(Ranges),]
        xC <- X[X[,2] == 0 & X[,3] >= min(Ranges) & X[,3] <= max(Ranges),]
      }
  
  
  return(list(xE = xE, xC = xC))
}


#### 4. EstimateThetaMvB: Function to estimate theta with Multivariate Bernoulli distribution ####
# Function to estimate success probailities with a Multivariate Bernoulli distribution (reference approach)
# Input: 
# - Y: n x Q matrix with multinomial responses.
# - PriorAlpha: Vector of length Q with prior hyperparameters.
# - nIt: Scalar. Number of draws from posterior distributions.

# Output:
# - mTheta: nIt x K matrix with bivariate Bernoulli probabilities. Currently supported for K=2 only.

EstimateThetaMvB <- function(Y, PriorAlpha, nIt){
  mPhi <- MCMCpack::rdirichlet(nIt, colSums(Y) + PriorAlpha)
  mTheta <- cbind(rowSums(mPhi[,c(1,2)]), rowSums(mPhi[,c(1,3)]))
  return(list(MVB = mTheta))
}

#### 5. EstimateThetaEmpirical: Function to estimate theta via empirical marginalization ####
# Function to estimate success probabilities via empirical marginalization
# Input:
# - BetaPG: P x Q x nIt array of posterior draws of regression coefficients
# - X: n x P matrix with covariate data

# Output:
# - mTheta: nIt x K matrix with bivariate probabilities. Currently supported for K=2 only.

EstimateThetaEmpirical <- function(BetaPG, X){
  nIt <- length(BetaPG)
  
  mPhi <- lapply(1:nIt, function(i) colMeans(exp(X %*% BetaPG[[i]]) / rowSums(exp(X %*% BetaPG[[i]]))))
  mThetaPG <- lapply(1:nIt, function(i) c(sum(mPhi[[i]][c(1,2)]), sum(mPhi[[i]][c(1,3)])))
  
  return(mTheta = mThetaPG)
}

#### 6. IntegrandPG: Function to integrate in numerical marginalization ####
# Integrand function for integration over a range of a covariate in numerical marginalization. 
# Input:
# - x: value of covariate x
# - Beta: P x Q matrix of regression coefficients
# - MuX: Mean of distribution of covariate
# - SigmaX: Standard deviation of distribution of covariate
# - RangeX: Vector of lower and upper bound of range to integrate over.
# - Trt: Scalar. Value of treatment indicator
# - q: Scalar. Response category in 1 to Q

# Output:
# - Phi[,q]: Scalar. Joint response of response category q.

IntegrandPG <- function(x, Beta, MuX, SigmaX, RangeX, Trt, q){
  xInt <- cbind(1,Trt,x,Trt*x)
  Psi <- xInt %*% Beta
  pX <- msm::dtnorm(x, mean = MuX, sd = SigmaX, lower = RangeX[1], upper = RangeX[2], log=TRUE)
  Phi <- exp(Psi - log(rowSums(exp(Psi))) + pX)
  return(Phi[,q])
}

#### 7. EstimateThetaAnalytical: Function to integrate over a range of x ####
# Function to integrate over a range of x
# Input: 
# - BetaPG: P x Q x nIt array of posterior regression coefficients
# - X: n x P matrix of covariate data, where the covariate of interest is in the third column
# - Trt: Scalar. Value of treatment indicator.
# - RangeX: Vector of lower and upper bound of range to integrate over.

# Output: 
# - mTheta: nIt x K matrix with bivariate probabilities. Currently supported for K=2 only.

EstimateThetaAnalytical <- function(BetaPG, X, Trt, RangeX){
  nIt <- length(BetaPG)
  
  mThetaPG <- vector("list", nIt)
  mPhi <- array(NA, dim = c(1, ncol(BetaPG[[1]])))
  for(i in 1:nIt){
    for(q in 1:(ncol(BetaPG[[1]]) - 1)){
      mPhi[,q] <- integrate(IntegrandPG, lower = RangeX[1], upper = RangeX[2],
                            Beta = BetaPG[[i]], MuX = mean(X[,3]), SigmaX = sd(X[,3]), RangeX = RangeX, Trt = Trt, q = q)$value
    }
    mThetaPG[[i]] <- c(rowSums(mPhi[,c(1,2),drop=FALSE]), rowSums(mPhi[,c(1,3),drop=FALSE]))
  }
  
  return(mTheta = mThetaPG)
}


#### 8. Transform2Theta: Function to compute theta from regression coefficients via various methods ####
# Function to compute theta from regression coefficients via various methods.
# Input: 
# - BetaDrawPG: P x Q x nIt array with posterior regression coefficients
# - X: n x P matrix with covariate data
# - Y: n x Q matrix with multinomial response data
# - Types: R x 4 matrix with R conditions, including columns with 
#	-- MeasurementLevels ("Discrete", "Continuous"),
#	-- Methods ("Value" for vector of fixed values, "Empirical" for empirical marginalization, "Analytical" for numerical marginalization, "MvB" for Multivariate Bernoulli (unconditional), 
#	-- Populations ("Trial" for study population, "Intra_Lo" for small subpopulation within trial), 
#	--Scopes ("Old" for existing covariate data, "New" for a new sample of covariate data to be drawn within the function.)
# - Ranges: List of S ranges that defines the population of interest by a vector of a lower and an upper bound. 
# - Values: List of S values that defines the population of interest by a scalar value.
# - PriorAlpha: Vector of length Q with prior hyperparameters. Currently supported for Q = 4 only.
# - MeasurementLevel: String. "Concrete" or "Discrete"

# Output:
# - List of:
# -- mTheta.E: nIt x K matrix with multivariate probabilities for experimental treatment (T=1). Currently supported for K=2 only.
# -- mTheta.C: nIt x K matrix with multivariate probabilities for control treatment (T=0). Currently supported for K=2 only.

Transform2Theta <- function(BetaDrawPG, X, Y = NULL, 
                            Types, Ranges, Values, PriorAlpha, MeasurementLevel){
  IndicesML <- which(Types[,"MeasurementLevels"] == MeasurementLevel)

  Empirical <- which(apply(Types, 1, function(x) any(c("Value", "Empirical") %in% x["Methods"])))
  Analytical <- which(apply(Types, 1, function(x) c("Analytical") %in% x["Methods"]))
  MvB <- which(apply(Types, 1, function(x) c("MvB") %in% x["Methods"]))
  
  xEmpirical <- lapply(1:length(IndicesML), function(x) vector("list", 2))
  mTheta.E <- mTheta.C <- vector("list", length(IndicesML))
  
  for(i in IndicesML){
    
    if(i %in% Empirical){
      xEmpirical[[which(IndicesML == i)]] <- SamplePopulation(X = X, 
                                                              Method = Types[i,"Methods"], 
                                                              Values = Values[[Types[i,"Populations"]]], 
                                                              Range = Ranges[[Types[i,"Populations"]]])
      mTheta.E[[which(IndicesML == i)]] <- EstimateThetaEmpirical(BetaPG = BetaDrawPG, X = xEmpirical[[which(IndicesML == i)]][["xE"]])
      mTheta.C[[which(IndicesML == i)]] <- EstimateThetaEmpirical(BetaPG = BetaDrawPG, X = xEmpirical[[which(IndicesML == i)]][["xC"]])
      
    } else if(i %in% Analytical){
      mTheta.E[[which(IndicesML == i)]] <- EstimateThetaAnalytical(BetaPG = BetaDrawPG, X = X, Trt = 1, RangeX = Ranges[[Types[i,"Populations"]]])
      mTheta.C[[which(IndicesML == i)]] <- EstimateThetaAnalytical(BetaPG = BetaDrawPG, X = X, Trt = 0, RangeX = Ranges[[Types[i,"Populations"]]])
      
    } else if(i %in% MvB){
      if(Types[i,"MeasurementLevels"] == "Discrete"){
        Indices <- which(X[,3] %in% Ranges[[Types[i,"Populations"]]])
      }
      else if(Types[i,"MeasurementLevels"] == "Continuous"){
        Indices <- which(X[,3] > min(Ranges[[Types[i,"Populations"]]]) & X[,3] < max(Ranges[[Types[i,"Populations"]]]))
      }
      mTheta.E[[which(IndicesML == i)]] <- EstimateThetaMvB(Y[intersect(Indices, which(X[,2] == 1)),], PriorAlpha, nIt = nIt)
      mTheta.C[[which(IndicesML == i)]] <- EstimateThetaMvB(Y[intersect(Indices, which(X[,2] == 0)),], PriorAlpha, nIt = nIt)
      
    }
  }
  
  return(list(mTheta.E = mTheta.E, mTheta.C = mTheta.C))
}

#### 9. EvaluateData: Function to evaluate data ####
# Function to compute bias and posterior probabilities 
# Input:
# - Data: nIt x K matrix with treatment differences (E - C). Currently supported for K=2 only. 
# - Weights: Vector of length K with weights for linear combination of treatment differences (Compensatory rule). Currently supported for K=2 only.
# - Alpha: Scalar. Desired Type I error rate. 
# - Alternative: String. Type of decision. "greater.than" for a right-sided test, "smaller.than" for a left-sided test, and "two.sided" for a two-sided test. 
# - Rule: String. Options: "All", "Any", "Compensatory" 
# - Truth: List of length 2 with true values: "Delta" is a vector of multivariate treatment differences and "DeltaW" is a scalar of a weighted treatment difference.
# - Types: R x 4 matrix with R conditions, including columns with 
#	-- MeasurementLevels ("Discrete", "Continuous"),
#	-- Methods ("Value" for vector of fixed values, "Empirical" for empirical marginalization, "Analytical" for numerical marginalization, "MvB" for Multivariate Bernoulli (unconditional), 
#	-- Populations ("Trial" for study population, "Intra_Lo" for small subpopulation within trial), 
#	-- Scopes ("Old" for existing covariate data. Currently no other option supported.)


# Output: 
# - Res: List of:
# -- Bias: List of length 2 with one vector of bias on individual treatment differences and one scalar of bias on the weighted linear combination.
# -- Pop: List of length 2 with one vector of posterior probabilities on individual treatment differences and one scalar of a posterior probability on the weighted linear combination. 
# -- Decision: Logical vector of S decisions.

EvaluateData <- function(Data, Weights, Alpha, Alternative, Rule, Truth, Types){
  DecisionTypes <- c("DecisionAny.RS", "DecisionAny.LS", "DecisionAny.TS", "DecisionAll.RS", "DecisionAll.LS","DecisionAll.TS", "DecisionCompensatory.RS", "DecisionCompensatory.LS", "DecisionCompensatory.TS")
  BiasTypes <- c("BiasMultivariate", "BiasWeighted")
  PopTypes <- c("PopMultivariate", "PopWeighted")
  
  if(any(c("Any", "All") %in% Rule)){
    BiasMultivariate <- t(apply(Data, 1, function(x) x - Truth[["Delta"]]))
    PopMultivariate <- colMeans(Data > 0)
    
    if("Any" %in% Rule){
      if("greater.than" %in% Alternative){DecisionAny.RS <- max(PopMultivariate) > (1 - Alpha / length(PopMultivariate))} 
      if("smaller.than" %in% Alternative){DecisionAny.LS <- min(PopMultivariate) < (Alpha / length(PopMultivariate))}
      if("two.sided" %in% Alternative){DecisionAny.TS <- max(PopMultivariate) > (1 - Alpha / (length(PopMultivariate) * 2)) | min(PopMultivariate) < (Alpha / (length(PopMultivariate) * 2))}
    }
    
    if("All" %in% Rule){
      if("greater.than" %in% Alternative){DecisionAll.RS <- min(PopMultivariate) > (1 - Alpha)} 
      if("smaller.than" %in% Alternative){DecisionAll.LS <- max(PopMultivariate) < (Alpha)} 
      if("two.sided" %in% Alternative){DecisionAll.TS <- min(PopMultivariate) > (1 - Alpha / 2) | max(PopMultivariate) < (Alpha / 2)}
    }
  }
  
  if("Compensatory" %in% Rule){
    DataWeighted <- Data %*% Weights
    BiasWeighted <- as.matrix(apply(DataWeighted, 1, function(x) x - Truth[["DeltaW"]]))
    PopWeighted <- colMeans(as.matrix(DataWeighted) > 0)
    
    if("greater.than" %in% Alternative){DecisionCompensatory.RS <- PopWeighted > (1 - Alpha)} 
    if("smaller.than" %in% Alternative){DecisionCompensatory.LS <- PopWeighted < (Alpha)} 
    if("two.sided" %in% Alternative){DecisionCompensatory.TS <- PopWeighted > (1 - Alpha / 2)| PopWeighted < (Alpha / 2)}
  }
  
  Bias <- sapply(BiasTypes, function(x){ y <- get0(x); if(is.array(y)){colMeans(y)}}, USE.NAMES = TRUE)
  Pop <- sapply(PopTypes, function(x) get0(x), USE.NAMES = TRUE)
  Decision <- sapply(DecisionTypes, function(x) get0(x), USE.NAMES = TRUE)
  
  Res <- list(Bias = Bias[!sapply(Bias,is.null)], Pop = Pop[!sapply(Pop,is.null)], Decision = Decision[!sapply(Decision,is.null)])
  
  return(Res)
}



#### 10. RoundChoose ####
# Function to round number up or down to a chosen interval
# Input:
# x: Scalar. Number to be rounded.
# roundTo: Scalar. Interval to be rounded to. E.g. 5, to round to the next 5th number
# dir: "1" for rounding up; "0" for rounding down. Defaults to 1.

# Output:
# roundedX: Scalar. Rounded number.
RoundChoose <- function(x, roundTo, dir = 1) {
  if(dir == 1) {  ##ROUND UP
    roundedX <- x + (roundTo - x %% roundTo)
  } else {
    if(dir == 0) {  ##ROUND DOWN
      roundedX <- x - (x %% roundTo)
    }
  }
  return(roundedX)
}