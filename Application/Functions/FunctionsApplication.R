#### 1. Generate data ####
# Function to generate covariate and response data 
# Input: 
# - TruePar: List with:
# -- b = (PF x Q) array with true parameters of fixed effects(regression coefficients)
# -- g = (PF x Q) array with true parameters of random effects(regression coefficients)
# -- sigma = scalar with true variance parameter parameters of random effects
# - Q: Scalar number of joint response categories (=2^K). 
# - pX: Scalar probability that discrete covariate = 1
# - MuX: Scalar mean of continuous covariate
# - SigmaX: Scalar standard deviation of continuous covariate 
# - n: Scalar sample size
# - J: Scalar number of clusters
# - RangeX: Vector with values of interval of discrete variables or (2) vector with limits of low interval of continuous covariate 
# - Continuous: Logical. TRUE if covariate is continuous; FALSE if covariate is discrete.
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.
# Output:
# - Data: 
# -- List with:
# -- X: List of J (n_j x Q) design matrices (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# -- Y: List of J (n_j x Q) responses. 
# -- Phi: List of J (n_j x Q) response probabilities.
# -- RC: List of J (p x Q) regression coefficients per cluster j.

GenerateData <- function(TruePar, Q, pX = NULL, Mux = NULL, SigmaX = NULL, J, n, Continuous, Fixed, Random){

  Intercept <- lapply(1:J, function(j) rep(1, n * 2))
  Trt <- lapply(1:J, function(j) rep(c(0,1), n))
  if(Continuous){x <- lapply(1:J, function(j) rnorm(n * 2, MuX, SigmaX))
  }else{x <- lapply(1:J, function(j) rbinom(n * 2, 1, pX))}
  Trt_x <- Map("*", Trt, x)
  
  if(length(Fixed) > 0){
  XF <- mget(Fixed)
  Xf <- lapply(1:J, function(j) do.call(cbind, lapply(XF, "[[", j)))
  }else{Xf <- vector("list", J)}
  
  if(length(Random) > 0){
   XR <- mget(Random)
  Xr <- lapply(1:J, function(j) do.call(cbind, lapply(XR, "[[", j)))
   gj <- array(0, dim = c(J,Q,length(Random)))

   for(q in 1:(Q-1)){
      gj[,q,] <- mvrnorm(n = J, mu = TruePar[["g"]][,q,drop=FALSE], Sigma = diag(TruePar[["Sigma"]], length(Random)))}
    dimnames(gj)[[3]] <- Random
  }else{Xr <- vector("list", J)}
  
  X <- Map("cbind", Xf, Xr)
  lapply(X, function(x) colnames(x) <- c(Fixed, Random))
  
  RC <- lapply(1:J, function(j){x <- rbind(if(exists("b", where = TruePar)){TruePar[["b"]]}, 
                                           if(exists("gj")){matrix(gj[j,,], ncol = Q, byrow = TRUE)})
     dimnames(x) <- list(c(Fixed, Random), paste0("q", 1:Q))
  return(x)
  })
  
  psi <- lapply(1:J, function(j) lapply(1:Q, function(q) exp(X[[j]] %*% RC[[j]][,q])))
  phi <- lapply(1:J, function(j) sapply(1:Q, function(q) psi[[j]][[q]] / Reduce("+", psi[[j]])))
  
  Y <- lapply(1:J, function(j) t(sapply(1:nrow(phi[[j]]), function(i) rmultinom(1, 1, phi[[j]][i,]))))
  
      Data <- list(X = X, Y = Y, Phi = phi, RC = RC)
  return(Data)
}



#### 2. SampleBetaPG_ML ####
# Function to sample regression coefficients via Polya-Gamma multinomial logistic regression
# Input:
# - X: List of J (n_j x Q) design matrices (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# - Y: List of J (n_j x Q) responses. 
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.
# - nBurn: Scalar. Number of burnin iterations
# - nIt: Scalar. Number of iterations
# - Start: Vector of two starting values, each used for a different chain
# - bMu0: Vector of prior means of fixed regression coefficients, where PF = no. of fixed covariates.
# - bSigma0: (PF x PF) Prior variance matrix of fixed regression coefficients, where PF = no. of fixed covariates.
# - gMu0 Vector of prior means of fixed regression coefficients, where PR = no. of random covariates.
# - gSigma0: (PR x PR) Prior covariance matrix of fixed regression coefficients, where PR = no. of random covariates.
# - nu0: Scalar. Degrees of freedom of prior inverse-Wishart distribution.
# - Tau0: (PR x PR) prior matrix of inverse-Wishart distribution.
# - ReturnThinned: Logical. Should thinned chains be returned? Default = FALSE.
# - nThin: Thinning rate. Default is 1. Adjustable to integers with a value large than one, only when ReturnThinned = TRUE. 

# Output:
# A list with:
# - bDrawPG: List of length nIt/nThin with PF x Q matrices with estimated fixed regression coefficients
# - gDrawPG: List of length nIt/nThin with PR x Q matrices with estimated random regression coefficients
# - gjDrawPG: List of length nIt/nThin with J lists of PR x Q matrices with estimated random regression coefficients (only when ReturnThinned = TRUE)
# - tauDrawPG: List of length nIt/nThin with PR x PR covariance matrices of random regression coefficients.

SampleBetaPG_ML <- function(X, Y, Fixed, Random, nBurn, nIt, Start, bMu0 = NULL, bSigma0 = NULL, gMu0 = NULL, gSigma0 = NULL, nu0 = NULL, Tau0 = NULL, ReturnThinned = FALSE, nThin = 1){
  out <- tryCatch(
    {
    J <- length(X)    
    z <- lapply(1:J, function(j) Y[[j]] - 1/2)
   
       pF <- length(Fixed)
    if(pF > 0){IndF <- 1:pF
    XF <- lapply(1:J, function(j) matrix(X[[j]], ncol = length(c(Fixed, Random)))[,IndF,drop=FALSE])
     pB0 <- bSigma0 %*% bMu0
     }
    pR <- length(Random)
    if(pR > 0){IndR <- 1:pR + pF
    XR <-  lapply(1:J, function(j) matrix(X[[j]], ncol = length(c(Fixed, Random)))[,IndR,drop=FALSE])
    pG0 <- gSigma0 %*% gMu0
    }
 
    n <- sapply(Y, nrow);
    Q <- ncol(Y[[1]]);
    

    if(pF > 0){
    bDrawPG <- vector("list", nIt)
    }
    if(pR > 0){
      gDrawPG <- gjDrawPG <- tauDrawPG <- omegaDrawPG <- w.gDrawPG <- list("vector", nIt)
    }

    
 omegaDraw <- lapply(1:(Q-1), function(q) lapply(1:J, function(j) rep(Start, n[j])));
  
   if(pF > 0){
 bDraw <- cbind(matrix(Start, nrow = pF, ncol = Q-1), 0)}
 
  if(pR > 0){
     gDraw <- cbind(matrix(Start, nrow = pR, ncol = Q-1), 0);
    gjDraw <- w.gDraw <- abind::abind(array(Start, dim = c(pR,Q-1,J)), matrix(0,pR,J), along = 2);
    varMat <- lapply(1:(Q-1), function(q) diag(Start, pR));
    tauMat <- lapply(varMat, function(x) chol2inv(chol(x)));
  }
 
 
  
   # create progress bar
  pb <- tkProgressBar(title = "progress bar", min = 1,
                      max = nBurn + nIt, width = 300)
  
  
  # Start Gibbs sampler
  for(i in 2:(nBurn + nIt)){
    for(q in 1:(Q-1)){
      
      # Draw auxiliary variable (i.e. variance of kappa)
      RC <- lapply(1:J, function(j) rbind(if(pF > 0){bDraw}, if(pR > 0){matrix(gjDraw[,,j], nrow = pR, ncol = Q)}))
      C <- lapply(1:J, function(j) log(rowSums(exp(X[[j]] %*% RC[[j]][,-q,drop = FALSE]))))
      omegaDraw[[q]] <- lapply(1:J, function(j){pgdraw(1,X[[j]] %*% RC[[j]][,q,drop=FALSE] - C[[j]])})
      

       # Compute linear predictor
      if(pF > 0){
        if(pR > 0){
          Cb <- lapply(1:J, function(j)  - XR[[j]] %*% matrix(gjDraw[,q,j], nrow = pR) + 
                             log(rowSums(exp(X[[j]] %*% RC[[j]][,-q]))))
        }else{
         Cb <- C
         }
       
      # Draw regression coefficients
      bSigma <- chol2inv(chol(Reduce('+', lapply(1:J, function(j) t(XF[[j]]) %*% (XF[[j]] * omegaDraw[[q]][[j]]))) + bSigma0))
      bMu <- bSigma %*% (Reduce('+', lapply(1:J, function(j) t(XF[[j]]) %*% (z[[j]][,q,drop=FALSE] + omegaDraw[[q]][[j]] * Cb[[j]]))) + pB0) 
      bDraw[,q] <- bMu + t(chol(bSigma)) %*% rnorm(pF) 
      }
      
      if(pR > 0){
        if(pF > 0){
        Cg <- lapply(1:J, function(j) - XF[[j]] %*% bDraw[,q,drop=FALSE] + 
                         log(rowSums(exp(X[[j]] %*% rbind(bDraw[,-q,drop=FALSE], matrix(gjDraw[,-q,j], nrow = pR, ncol = Q-1))))))
        }else{
          Cg <- C
        }
      
         for(j in 1:J){
        gjSigma <- chol2inv(chol(t(XR[[j]]) %*%  (XR[[j]] * omegaDraw[[q]][[j]]) + tauMat[[q]])) 
        gjMu <- gjSigma %*% (t(XR[[j]]) %*% (z[[j]][,q,drop=FALSE] + omegaDraw[[q]][[j]] * Cg[[j]]) + tauMat[[q]] %*% gDraw[,q,drop=FALSE]) 
        gjDraw[,q,j] <- c(gjMu + t(chol(gjSigma)) %*% rnorm(pR))
                }
        
        gSigma <- chol2inv(chol(J * tauMat[[q]] + gSigma0)) 
        gMu <- gSigma %*% (tauMat[[q]] %*% rowSums(matrix(gjDraw[,q,], nrow = pR, ncol = J)) + pG0) 
        gDraw[,q] <- gMu + t(chol(gSigma)) %*% rnorm(pR) 
        
        tauShape <- nu0 + J 
        e <- matrix(NA, pR, J)
        for(j in 1:J){e[,j] <- gjDraw[,q,j] - gDraw[,q]}
        tauScale <- Tau0 + e %*% t(e)
        varMat[[q]] <- riwish(tauShape, tauScale)
        tauMat[[q]] <- chol2inv(chol(varMat[[q]])) # Precision 
        
        
     }
        
    } 
    
      if(i > nBurn){

  if(pF > 0){
      bDrawPG[[i-nBurn]] <- bDraw}
    if(pR > 0){
      gDrawPG[[i-nBurn]] <- gDraw
      gjDrawPG[[i-nBurn]] <- gjDraw
      tauDrawPG[[i-nBurn]] <- varMat
      }
   } 
    setTkProgressBar(pb, i, label=paste("Iteration", i, " of ", nBurn + nIt)) 
  }
  close(pb)

  if(ReturnThinned){
  if(pF > 0){
    Thinned.bDrawPG <- lapply(seq(1,nIt,nThin), function(i) {bDrawPG[[i]]})
    Out.Fixed <- list(bDrawPG = Thinned.bDrawPG)}
  if(pR > 0){
    Thinned.gDrawPG <- lapply(seq(1,nIt,nThin), function(i) {gDrawPG[[i]]})
  Thinned.gjDrawPG <- lapply(seq(1,nIt,nThin), function(i) {gjDrawPG[[i]]})
  Thinned.tauDrawPG <- lapply(seq(1,nIt,nThin), function(i) {tauDrawPG[[i]]})
  Out.Random <- list(gDrawPG = Thinned.gDrawPG, gjDrawPG = Thinned.gjDrawPG, tauDrawPG = Thinned.tauDrawPG)}
  } else {
    if(pF > 0){Out.Fixed <- list(bDrawPG = bDrawPG)} 
    if(pR > 0){Out.Random <- list(gDrawPG = gDrawPG, tauDrawPG = tauDrawPG)}
    
  }
return(c(if(pF > 0){Out.Fixed}, if(pR > 0){Out.Random}))
    },
error=function(e) {
  message(paste("Error"))
  #if(exists(i)){paste("in iteration", i)}
  message("Here's the original error message:")
  message(e)
  # Choose a return value in case of error
  return(NA)
})
  return(out)
  
}

#### 3. EstimateParameters ####
# Function to estimate regression coefficients using Polya-Gamma gibbs sampling. 
# Input:
# - X: List of J (n_j x Q) design matrices (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# - Y: List of J (n_j x Q) responses. 
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.
# - nBurn: Scalar. Number of burnin iterations
# - nIt: Scalar. Number of iterations
# - Start: Vector of two starting values, each used for a different chain
# - bMu0: Vector of prior means of fixed regression coefficients, where PF = no. of fixed covariates.
# - bSigma0: (PF x PF) Prior variance matrix of fixed regression coefficients, where PF = no. of fixed covariates.
# - gMu0 Vector of prior means of fixed regression coefficients, where PR = no. of random covariates.
# - gSigma0: (PR x PR) Prior covariance matrix of fixed regression coefficients, where PR = no. of random covariates.
# - nu0: Scalar. Degrees of freedom of prior inverse-Wishart distribution.
# - Tau0: (PR x PR) prior matrix of inverse-Wishart distribution.
# - nChain: Scalar. Number of chains. 
# - ReturnThinned: Logical. Should thinned chains be returned? Default = FALSE.
# - nThin: Thinning rate. Default is 1. Adjustable to integers with a value large than one, only when ReturnThinned = TRUE. 

# Output: 
# A list with:
# - Pars: List of length nChain sublists of:
# -- bDrawPG: List of length nIt/nThin with PF x Q matrices with estimated fixed regression coefficients
# -- gDrawPG: List of length nIt/nThin with PR x Q matrices with estimated random regression coefficients
# -- gjDrawPG: List of length nIt/nThin with J lists of PR x Q matrices with estimated random regression coefficients (only when ReturnThinned = TRUE)
# -- tauDrawPG: List of length nIt/nThin with PR x PR covariance matrices of random regression coefficients.
# - Convergence: Multivariate Gelman-Rubin statistic to asses convergence

EstimateParameters <- function(X, Y, Fixed, Random, nBurn, nIt, Start, bMu0 = NULL, bSigma0 = NULL, gMu0 = NULL, gSigma0 = NULL, nu0 = NULL, Tau0 = NULL, nChain, ReturnThinned = FALSE, nThin = 1){

  Chain <- vector("list", nChain)
  for(chain in 1:nChain){
   Chain[[chain]] <- SampleBetaPG_ML(X = X, Y = Y, Fixed = Fixed, Random = Random, 
                                    nBurn = nBurn, nIt = nIt, Start = Start[chain], 
                                    bMu0 = bMu0, bSigma0 = bSigma0, 
                                    gMu0 = gMu0, gSigma0 = gSigma0,
                                    nu0 = nu0, Tau0 = Tau0, ReturnThinned = ReturnThinned, nThin = nThin)
  
  
  }
 
  if(length(Fixed) > 0){
  Chains.bDraw <- lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Chain[[chain]][["bDrawPG"]]), nrow = nIt/nThin, ncol = length(Fixed) * Q, byrow = TRUE)[,-((length(Fixed) * (Q-1)) + 1:(length(Fixed)))])) #[,-((Q-1)*length(Fixed)+1:length(Fixed))])) 
  Convergence.bDraw <- gelman.diag(Chains.bDraw)$mpsrf
  Convergence.Fixed <- Convergence.bDraw
  }
 
  if(length(Random) > 0){
    Chains.tDraw <- lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(lapply(1:(nIt/nThin), function(i) lapply(1:(Q-1), function(q) {x <- Chain[[chain]][["tauDrawPG"]][[i]][[q]]; x[lower.tri(x, diag = TRUE)]}))), nrow = nIt/nThin, byrow=TRUE)))
  Chains.gDraw <- lapply(1:nChain, function(chain) as.mcmc(matrix(unlist(Chain[[chain]][["gDrawPG"]]), nrow = nIt/nThin, ncol = length(Random) * Q, byrow = TRUE)[,-((length(Random) * (Q-1)) + 1:(length(Random)))]))
    Convergence.gDraw <- gelman.diag(Chains.gDraw)$mpsrf
  Convergence.tDraw <- gelman.diag(Chains.tDraw)$mpsrf
  Convergence.Random <- list(gDraw = Convergence.gDraw, tDraw = Convergence.tDraw)}
  
  Convergence <- c(if(length(Fixed) > 0){Convergence.Fixed}, if(length(Random) > 0){Convergence.Random})

    return(list(Pars = Chain, Convergence = Convergence))
}



#### 4. EstimateBiasPars ####
# Input:
# - Pars: List of length nChain sublists of:
# -- bDrawPG: List of length nIt with PF x Q matrices with estimated fixed regression coefficients
# -- gDrawPG: List of length nIt with PR x Q matrices with estimated random regression coefficients
# -- tauDrawPG: List of length nIt with PR x PR covariance matrices of random regression coefficients.
# - Convergence: Multivariate Gelman-Rubin statistic to asses convergence
# - TruePar: List with:
# -- b = (PF x Q) array with true parameters of fixed effects(regression coefficients)
# -- g = (PF x Q) array with true parameters of random effects(regression coefficients)
# -- sigma = scalar with true variance parameter parameters of random effects
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.
# - nPars: Scalar. No of iterations used in computations. Usually nIt, but can be thinned. 
# Output:
# - Bias: A list with:
# -- bDraw: A (PF x Q) matrix with bias of fixed regression coefficients. 
# -- gDraw: A (PR x Q) matrix with bias of random regression coefficients. 
# -- tDraw: A (PR x PR) matrix with bias of the covariance matrix of random regression coefficients.
EstimateBiasPars <- function(Pars, TruePar, Fixed, Random, nPars){
  if(length(Fixed) > 0){
  Bias.bDrawTrue <- Reduce("+", lapply(1:nPars, function(i) Pars[["bDrawPG"]][[i]] - TruePar[["b"]])) / nPars
    Bias.Fixed <- list(bdrawTrue = Bias.bDrawTrue)}
  if(length(Random) > 0){
  Bias.gDraw <- Reduce("+", lapply(1:nPars, function(i) Pars[["gDrawPG"]][[i]] - TruePar[["g"]])) / nPars
  
  Bias.tDraw <- Reduce("+", lapply(1:nPars, function(i) sapply(1:(Q-1), function(q) Pars[["tauDrawPG"]][[i]][[q]] - diag(TruePar[["Sigma"]], nrow(TruePar[["g"]]))))) / nPars
  Bias.Random <- list(gDraw = Bias.gDraw, tDraw = Bias.tDraw)}
  
  
  Bias <- c(if(length(Fixed) > 0){Bias.Fixed}, if(length(Random) > 0){Bias.Random})
  
  return(Bias)
}


#### 5. SamplePopulation: Function to (sub)sample covariate data ####
# Input: 
# - X: (n x Q) design matrix (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# - Method: "Empirical" for empirical marginalization, "Value" for vector of fixed values, "MvB" for unconditional multivariate Bernoulli (reference approach)
# - Values: Scalar. Value of x representing subpopulation.
# - Range: Vector of lower and upper bound of covariate that represents subpopulation.
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.

# Output: 
# List of:
# -- xE: nE x P design matrix with covariate data for treatment E
# -- xC: nC x P design matrix with covariate data for treatment C
SamplePopulation <- function(X, Method, Values, Range, Fixed, Random){
  
  if("Value" %in% Method){
    x <- Values
    Intercept <- 1
    Trt <- 1
    Trt_x <- Trt * Values
    xValsE <- c(Fixed, Random)
    xE <- t(do.call(rbind, mget(xValsE)))
    
    Trt <- 0
    Trt_x <- Trt * Values
    xValsC <- c(Fixed, Random)
    xC <- t(do.call(rbind, mget(xValsC)))
    
  } else if(any(c("MvB", "Empirical") %in% Method)){
    if(ncol(X) > 2){
    xE <- X[X[,"Trt"] == 1 & X[,"x"] >= min(Range) & X[,"x"] <= max(Range),]
    xC <- X[X[,"Trt"] == 0 & X[,"x"] >= min(Range) & X[,"x"] <= max(Range),]
  }else{
    xE <- X[X[,"Trt"] == 1,]
    xC <- X[X[,"Trt"] == 0,] 
  }
  }
  
  return(list(xE = xE, xC = xC))
}

#### 6. EstimateThetaMvB: Function to estimate theta with Multivariate Bernoulli distribution ####
# Function to estimate success probailities with a Multivariate Bernoulli distribution (reference approach)
# Input: 
# - Y: n x Q matrix with multinomial responses.
# - PriorAlpha: Vector of length Q with prior hyperparameters.
# - nPars: Scalar. No of iterations used in computations. Usually nIt, but can be thinned. 

# Output:
# - mTheta: nIt x K matrix with bivariate Bernoulli probabilities. Currently supported for K=2 only.

EstimateThetaMvB <- function(Y, PriorAlpha, nPars){
  mPhi <- MCMCpack::rdirichlet(nIt, colSums(Y) + PriorAlpha)
  
  if(ncol(Y) == 4){
  mTheta <- cbind(rowSums(mPhi[,c(1,2)]), rowSums(mPhi[,c(1,3)]))
  }else if(ncol(Y) == 2){
    mTheta <- mPhi
  }
  return(mTheta)
}

#### 7. EstimateThetaEmpirical: Function to estimate theta via empirical marginalization ####
# Function to estimate success probabilities via empirical marginalization
# Input:
# - EstPars: List of nIt (P x Q) arrays of posterior draws of regression coefficients.
# - X: (n x P) matrix with covariate data

# Output:
# - mTheta: List of nIt vectors of length K matrix with multivariate probabilities. Currently supported for K=1 and K=2 only.

EstimateThetaEmpirical <- function(EstPars, X){
  nIt <- length(EstPars)
  mPhi <- lapply(1:nIt, function(i) exp(X %*% EstPars[[i]]) / rowSums(exp(X %*% EstPars[[i]])))
  if(length(mPhi[[1]]) == 4 | ncol(mPhi[[1]]) == 4){
  mTheta <- lapply(1:nIt, function(i) colMeans(cbind(rowSums(mPhi[[i]][,c(1,2),drop=FALSE]), rowSums(mPhi[[i]][,c(1,3),drop=FALSE]))))
  }else if(length(mPhi[[1]]) == 2 | ncol(mPhi[[1]]) == 2){
    mTheta <- colMeans(mPhi)
  }
  
  return(mTheta)
}

#### 8. IntegrandPG: Function to integrate in numerical marginalization ####
# Integrand function for integration over a range of a covariate in numerical marginalization. 
# Input:
# - x: value of covariate x
# - Beta: P x Q matrix of regression coefficients
# - MuX: Mean of distribution of covariate
# - SigmaX: Standard deviation of distribution of covariate
# - RangeX: Vector of lower and upper bound of range to integrate over.
# - Trt: Scalar. Value of treatment indicator
# - q: Scalar. Response category in 1 to Q
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.

# Output:
# - Phi[,q]: Scalar. Joint response of response category q.

Integrand <- function(x, EstRC, MuX, SigmaX, RangeX, Trt, q, Fixed, Random){
  Intercept <- 1
  Trt_x <- Trt * x
  xVals <- c(Fixed, Random)
  
  xInt <- t(do.call(rbind, mget(xVals)))
  Psi <- xInt %*% EstRC
  pX <- msm::dtnorm(x, mean = MuX, sd = SigmaX, lower = RangeX[1], upper = RangeX[2], log=TRUE)
  Phi <- exp(Psi - log(rowSums(exp(Psi))) + pX)
  return(Phi[,q])
}

#### 9. EstimateThetaAnalytical: Function to integrate over a range of x ####
# Function to integrate over a range of x
# Input: 
# - EstPars: List of nIt (P x Q) arrays of posterior regression coefficients, ordered by rows
# - X: (n x P) matrix of covariate data, where the covariate of interest is in the third column named "x".
# - Trt: Scalar. Value of treatment indicator.
# - RangeX: Vector of lower and upper bound of range to integrate over.
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.

# Output: 
# - mTheta: List of nIt vectors of length K with multivariate probabilities. Currently supported for K=2 only.

EstimateThetaAnalytical <- function(EstPars, X, Trt, RangeX, Fixed, Random){
  Q <- ncol(EstPars[[1]])
  nIt <- length(EstPars)
  
  mTheta <- vector("list", nIt)
  mPhi <- array(NA, dim = c(1, Q))
  for(i in 1:nIt){
    for(q in 1:(Q - 1)){
      mPhi[,q] <- integrate(Integrand, lower = RangeX[1], upper = RangeX[2],
                            EstRC = EstPars[[i]], MuX = mean(X[,"x"]), SigmaX = sd(X[,"x"]), RangeX = RangeX, Trt = Trt, q = q,
                            Fixed = Fixed, Random = Random)$value
    }
    if(Q == 4){
      mTheta[[i]] <- cbind(rowSums(mPhi[,c(1,2),drop=FALSE]), rowSums(mPhi[,c(1,3),drop=FALSE]))
    }else if(Q == 2){
      mTheta[[i]] <- mPhi
    }

  }
  
  return(mTheta)
}

#### 10. Transform2Theta: Function to compute theta from regression coefficients via various methods ####
# Function to compute theta from regression coefficients via various methods.
# Input: 
# - EstPars: List of nIt (P x Q) array with posterior regression coefficients
# - X: List of J (n_j x Q) design matrices (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# - Y: List of J (n_j x Q) responses.
# - nIt: Scalar. Number of iterations.
# - Types: (R x 4) matrix with R conditions, including columns with 
#	-- MeasurementLevels ("Discrete", "Continuous"),
#	-- Methods ("Value" for vector of fixed values, "Empirical" for empirical marginalization, "Analytical" for numerical marginalization, "MvB" for Multivariate Bernoulli (unconditional), 
#	-- Populations ("Trial" for study population, "Intra_Lo" for small subpopulation within trial), 
# - Range: List of S ranges that defines the population of interest by a vector of a lower and an upper bound. 
# - Value: List of S values that defines the population of interest by a scalar value.
# - Fixed: Character vector with names of fixed variables in covariate vector.
# - Random: Character vector with names of random variables in covariate vector.
# - PriorAlpha: Vector of length Q with prior hyperparameters. Currently supported for Q = 4 only.

# Output:
# - List of length R, storing for each condition r in 1,...,R:
# -- mTheta.E: nIt x K matrix with multivariate probabilities for experimental treatment (T=1). Currently supported for K=2 only.
# -- mTheta.C: nIt x K matrix with multivariate probabilities for control treatment (T=0). Currently supported for K=2 only.

Transform2Theta <- function(EstPars = NULL, X, Y = NULL, nIt, Types, Range, Value, Fixed, Random, PriorAlpha){
   nPars <- length(EstPars)
  J <- length(X)
  Empirical <- which(apply(Types, 1, function(x) any(c("Empirical", "Value") %in% x["Methods"])))
  Analytical <- which(apply(Types, 1, function(x) c("Analytical") %in% x["Methods"]))
  MvB <- which(apply(Types, 1, function(x) c("MvB") %in% x["Methods"]))
  
  xEmpirical <- lapply(1:nrow(Types), function(x) lapply(1:J, function(j) vector("list", 2)))
  Theta.E <- Theta.C <- lapply(1:nrow(Types), function(i) vector("list", J))
  mTheta.E <- mTheta.C <- vector("list", nrow(Types))#lapply(1:nrow(Types), function(type) vector("list", J))
  
  nj.E <- nj.C <- rep(NA, J)
  Indices <- vector("list", J)
   for(i in 1:nrow(Types)){
     for(j in 1:J){
     if(ncol(X[[j]]) <= 2){
       Indices[[j]] <- 1:nrow(X[[j]])    
     }
     else if(Types[i,"MeasurementLevels"] == "Discrete" & ncol(X[[j]]) > 2){
       Indices[[j]] <- which(X[[j]][,"x"] %in% Range[[Types[i,"Populations"]]])
     }
     else if(Types[i,"MeasurementLevels"] == "Continuous" & ncol(X[[j]]) > 2){
       Indices[[j]] <- which(X[[j]][,"x"] > min(Range[[Types[i,"Populations"]]]) & X[[j]][,"x"] < max(Range[[Types[i,"Populations"]]]))
     }
     nj.E[j] <- length(intersect(Indices[[j]], which(X[[j]][,"Trt"] == 1)))
     nj.C[j] <- length(intersect(Indices[[j]], which(X[[j]][,"Trt"] == 0)))
     }
    if(i %in% Empirical){
      for(j in 1:J){
       xEmpirical[[i]] <- SamplePopulation(X = X[[j]], Method = Types[i,"Methods"], Values = Value[[Types[i,"Populations"]]], Range = Range[[Types[i,"Populations"]]],
                                          Fixed = Fixed, Random = Random)
      Theta.E[[i]][[j]] <- EstimateThetaEmpirical(EstPars = lapply(1:(nIt), function(i) EstPars[[i]][[j]]), X = xEmpirical[[i]][["xE"]])
      Theta.C[[i]][[j]] <- EstimateThetaEmpirical(EstPars = lapply(1:(nIt), function(i) EstPars[[i]][[j]]), X = xEmpirical[[i]][["xC"]])
          }
      mTheta.E[[i]] <- lapply(1:(nIt), function(k) {colSums(do.call(rbind, lapply(which(nj.E > 0), function(j) Theta.E[[i]][[j]][[k]] * nj.E[j] / sum(nj.E))))})
      mTheta.C[[i]] <- lapply(1:(nIt), function(k) {colSums(do.call(rbind, lapply(which(nj.C > 0), function(j) Theta.C[[i]][[j]][[k]] * nj.C[j] / sum(nj.C))))})
        } else if(i %in% Analytical){
      if(ncol(X[[1]]) > 2){
         mTheta.E[[i]] <- EstimateThetaAnalytical(EstPars = EstPars[[j]], X = do.call(rbind, X), Trt = 1, RangeX = Range[[Types[i,"Populations"]]], Fixed = Fixed, Random = Random)
        mTheta.C[[i]] <- EstimateThetaAnalytical(EstPars = EstPars[[j]], X = do.call(rbind, X), Trt = 0, RangeX = Range[[Types[i,"Populations"]]], Fixed = Fixed, Random = Random)
      }else{
         mTheta.E[[i]] <- mTheta.C[[i]] <- NULL
       }
    } else if(i %in% MvB){
        yE <- do.call(rbind, lapply(which(nj.E > 0), function(j){
        Y[[j]][intersect(Indices[[j]], which(X[[j]][,"Trt"] == 1)),,drop=FALSE]}))
      yC <- do.call(rbind, lapply(which(nj.C > 0), function(j){
        Y[[j]][intersect(Indices[[j]], which(X[[j]][,"Trt"] == 0)),,drop=FALSE]}))
      
      theta.E <- EstimateThetaMvB(yE, PriorAlpha, nPars = nPars)
      theta.C <- EstimateThetaMvB(yC, PriorAlpha, nPars = nPars)
      
      mTheta.E[[i]] <- lapply(seq_len(nrow(theta.E)),function(x) theta.E[x,])
      mTheta.C[[i]] <- lapply(seq_len(nrow(theta.C)),function(x) theta.C[x,])
  
    }
  }
  return(list(mTheta.E = mTheta.E, mTheta.C = mTheta.C))
}


#### 11. EvaluateData: Function to evaluate data ####
# Function to compute bias and posterior probabilities 
# Input:
# - Data: List of nIt vectors of length K with treatment differences (E - C). Currently supported for K=2 only. 
# - Weights: Vector of length K with weights for linear combination of treatment differences (Compensatory rule). Currently supported for K=2 only.
# - Alpha: Scalar (0 =< Alpha =< 1). Desired Type I error rate. 
# - Alternative: String. Type of decision. "greater.than" (default) for a right-sided test and "two.sided" for a two-sided test. 
# - Rule: String. Options: "All", "Any", "Compensatory" 
# - Truth: List of length 2 with true values: "Delta" is a vector of multivariate treatment differences and "DeltaW" is a scalar of a weighted treatment difference.
# - Types: R x 4 matrix with R conditions, including columns with 
#	-- MeasurementLevels ("Discrete", "Continuous"),
#	-- Methods ("Value" for vector of fixed values, "Empirical" for empirical marginalization, "Analytical" for numerical marginalization, "MvB" for Multivariate Bernoulli (unconditional), 
#	-- Populations ("Trial" for study population, "Intra_Lo" for small subpopulation within trial), 
# - nIt = Scalar. Numbver of iterations

# Output: 
# - Res: List of:
# -- Bias: List of length 2 with one vector of bias on individual treatment differences and one scalar of bias on the weighted linear combination.
# -- Pop: List of length 2 with one vector of posterior probabilities on individual treatment differences and one scalar of a posterior probability on the weighted linear combination. 
# -- Decision: Logical vector of S decisions.

EvaluateData <- function(Data, Weights, Alpha, Alternative = "greater.than", Rule, Truth = NULL, Types, nIt){
  if("greater.than" %in% Alternative){DecisionTypes.RS <- c("DecisionAny.RS", "DecisionAll.RS", "DecisionCompensatory.RS")}
  if("two.sided" %in% Alternative){DecisionTypes.TS <- c("DecisionAny.TS", "DecisionAll.TS", "DecisionCompensatory.TS")}
  DecisionTypes <- unlist(mget(c("DecisionTypes.RS","DecisionTypes.TS"), ifnotfound = list(c(), c())), use.names = FALSE)
  #names(DecisionTypes) <- unlist(mget(c("names(DecisionTypes.RS)", "names(DecisionTypes.TS")))
  BiasTypes <- c("BiasMultivariate", "BiasWeighted")
  PopTypes <- c("PopMultivariate", "PopWeighted")
  
 # if(Types["Methods"]== "MvB"){
#    DeltaEst <- lapply(1:nIt, function(i) {lapply(1:length(nJ), function(j) Data[[j]][[i]])})
 # }else{
#  DeltaEst <- lapply(1:nIt, function(i) colSums(do.call(rbind, lapply(1:length(Data), function(j) nJ[[j]]/sum(nJ) * Data[[j]][[i]]))))
#}
  DeltaEst <- Data
  if(any(c("Any", "All") %in% Rule)){
    if(!is.null(Truth)){BiasMultivariate <- do.call(rbind, lapply(1:nIt, function(i) DeltaEst[[i]] - Truth[["Delta"]]))}
    PopMultivariate <- colMeans(do.call(rbind, lapply(1:nIt, function(i) DeltaEst[[i]] > 0)))
    
    if("Any" %in% Rule){
      if("greater.than" %in% Alternative){DecisionAny.RS <- max(PopMultivariate) > (1 - Alpha / length(PopMultivariate))} 
      if("two.sided" %in% Alternative){DecisionAny.TS <- max(PopMultivariate) > (1 - Alpha / (length(PopMultivariate) * 2)) | min(PopMultivariate) < (Alpha / (length(PopMultivariate) * 2))}
    }
    
    if("All" %in% Rule){
      if("greater.than" %in% Alternative){DecisionAll.RS <- min(PopMultivariate) > (1 - Alpha)} 
      if("two.sided" %in% Alternative){DecisionAll.TS <- min(PopMultivariate) > (1 - Alpha / 2) | max(PopMultivariate) < (Alpha / 2)}
    }
  }
  
  if("Compensatory" %in% Rule){
    DataWeighted <- do.call(rbind, lapply(1:nIt, function(i) DeltaEst[[i]] %*% Weights))
    if(!is.null(Truth)){BiasWeighted <- as.matrix(apply(DataWeighted, 1, function(x) x - Truth[["DeltaW"]]))}
    PopWeighted <- colMeans(as.matrix(DataWeighted) > 0)
    
    if("greater.than" %in% Alternative){DecisionCompensatory.RS <- PopWeighted > (1 - Alpha)} 
    if("two.sided" %in% Alternative){DecisionCompensatory.TS <- PopWeighted > (1 - Alpha / 2)| PopWeighted < (Alpha / 2)}
  }
  
  Bias <- sapply(BiasTypes, function(x){ y <- get0(x); if(is.array(y)){colMeans(y)}}, USE.NAMES = TRUE)
  Pop <- sapply(PopTypes, function(x) get0(x), USE.NAMES = TRUE)
  Decision <- sapply(DecisionTypes, function(x){y <- get0(x)}, USE.NAMES = TRUE)
  
  Res <- list(Bias = Bias[!sapply(Bias,is.null)], Pop = Pop[!sapply(Pop,is.null)], Decision = Decision[!sapply(Decision,is.null)])
  
  return(Res)
}


