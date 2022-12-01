#### Functions data ####
#### 1. ComputeRho: Compute pairwise correlation ####
# From Dai, Ding, & Wahba (2013). Multivariate Bernoulli distribution.
# Input: 
# - phi: Vector of length 4 with joint response probabilities.

# Output:
# - rho: Scalar. Correlation between outcome variables.
ComputeRho <- function(phi){
  rho <- (phi[1] * phi[4] - phi[2] * phi[3]) / 
    sqrt(sum(phi[c(1,2)]) * sum(phi[c(1,3)]) * sum(phi[c(3,4)]) * sum(phi[c(2,4)]))
  return(rho)
}

#### 2. Theta2Phi: Transform theta to phi ####
# Currently supported for K=2/Q=4 only.
# Input: 
# - theta: Vector of bivariate success probabilities 
# - rho: Scalar pairwise correlation between success probabilities in theta.

# Ouput:
# - phi: Vector of 4 joint response probabilities, summing to 1 and ordered as {11},{10},{01},{00}.
Theta2Phi <- function(theta,rho){
  phi11 <- rho * sqrt(prod(theta)*prod(1-theta)) + prod(theta)
  phi <- c(phi11, theta[1]-phi11, theta[2]-phi11,1-theta[1]-theta[2]+phi11)
  if(!all(phi>= 0)) {
    warning("One or more negative probabilities. Rho adjusted.")
    while(!all(phi >= 0)){
      if(rho < 0) rho <- rho + 0.01
      if(rho > 0) rho <- rho - 0.01
      phi11 <- rho * sqrt(prod(theta)*prod(1-theta)) + prod(theta)
      phi <- c(phi11, theta[1]-phi11, theta[2]-phi11,1-theta[1]-theta[2]+phi11)
    }
  }
  return(phi)
}

#### 3. Phi2Theta: Transform phi to theta ####
# Currently supported forK=2/Q=4 only.
# Input: 
# - phi: Vector of 4 joint response probabilities, summing to 1 and ordered as {11},{10},{01},{00}.

# Output:
# - theta: Vector of bivariate success probabilities.
Phi2Theta <- function(phi){
  c(sum(phi[c(1,2)]), sum(phi[c(1,3)]))
}

#### 4. Theta2DeltaW: Transform bivaraite success probabilities to weighted treatment difference.
# Currently supported for K=2 only.
# Input:
# - thetaE: Vector of bivariate success probabilities for experimental treatment.
# - thetaC: Vector of bivariate success probabilities for control treatment.
# - Weights: Vector of length K with weights for linear combination of treatment differences (Compensatory rule). 

# Output: 
# - deltaW: Scalar. Weighted treatment difference.

Theta2DeltaW <- function(thetaE, thetaC, Weights){
  deltaW <- (thetaE - thetaC) %*% Weights
  return(deltaW)}

#### 4. Phi2DeltaW: Transform bivariate success probabilities to weighted treatment difference. ####
# Currently supported for Q=4/K=2 only.
# Input:
# - phiE: Vector of Q joint response probabilities for experimental treatment, ordered as {11},{10},{01},{00}.
# - phiC: Vector of Q joint response probabilities for control treatment, ordered as {11},{10},{01},{00}.
# - Weights: Vector of length K with weights for linear combination of treatment differences (Compensatory rule). 

# Output: 
# - deltaW: Scalar. Weighted treatment difference.

Phi2DeltaW <- function(phiE, phiC, Weights){
  delta <- phiE - phiC
  deltaTheta <- cbind(rowSums(delta[,c(1,2)]), rowSums(delta[,c(1,3)]))
  deltaW <- deltaTheta %*% Weights
  return(deltaW)
}


#### 5. Multivariate2Multinomial: Transform bivariate binomial response to multinomial response ####
# Input: 
# yBV: n x 2 matrix with bivariate binomial responses

# Output: 
# yMult: n x 4 matrix with multinomial responses, ordered as {11},{10},{01},{00}
Multivariate2Multinomial <- function(yBV){
  Answers <- rev(expand.grid(rev(list(c(0,1), c(0,1)))))
  Answers <- Answers[nrow(Answers):1,]
  yMultVec <- apply(yBV, 1, function(x) which(apply(Answers, 1, function(y) all(y == x))))
  yMult <- matrix(0, length(yMultVec), 2^K)
  for(i in 1:length(yMultVec)){yMult[i,yMultVec[i]] <- 1}
  return(yMult) 
}

#### 6. RoundChoose ####
# Adapted from: https://stackoverflow.com/a/32508105
# Function to round number up or down to a chosen interval
# Input:
# - x: Scalar. Number to be rounded.
# - roundTo: Scalar. Interval to be rounded to. E.g. 5, to round to the next 5th number
# - dir: "1" for rounding up; "0" for rounding down. Defaults to 1.

# Output:
# - roundedX: Scalar. Rounded number.
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
