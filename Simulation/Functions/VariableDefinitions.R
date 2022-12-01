# Variable definitions
Fixed <- c("x", "Trt_x") # Fixed covariates hierarchical model
Random <- c("Intercept", "Trt") # Random covariates hierarchical model
Fixed_NH <- c("Intercept", "Trt", "x", "Trt_x") # Fixed covariates non-hierarchical model
Random_NH <- c(NULL) # Random covariates non-hierarchical model (non-existing, but needed for function)

Q <- 4 # No. of joint response categories
Continuous <- TRUE # Continuous covariate?

pF <- length(Fixed) # No. of fixed regression coefficients
pR <- length(Random) # No. of random regression coefficients

# Covariate data - Continuous 
MuX <- 0			# Mean
SigmaX <- 1		# Standard deviation

# Define subpopulations 
# Range defining subpopulations
Ranges <- list()
Ranges[["Continuous"]] <- list(Trial = c(-Inf, Inf),
                               Intra_Lo = c(-1, 0),
                               Intra_Hi = c(0, 1))

# Values defining subpopulations
Values <- list()
Values[["Continuous"]] <- list(Trial = c(0),
                               Intra_Lo = c(-1),
                               Intra_Hi = c( 1))
                            
# Probability of discrete covariate
Sigma <- c(1e-1) # True variance of random effects (identical for all random effects)
nCluster <- c(10,100) # No. of clusters
nObs <- c(10,100) # No. of observations per cluster
SampleSizes <- cbind(expand.grid(nCluster, nObs), c(outer(nCluster, nObs))) # Matrix of sample sizes
colnames(SampleSizes) <- c("nCluster", "nObs", "Total")

#### Prior parameters ####
# Variance parameter
B0 <- 10  # Fixed effects
G0 <- 10  # Random effects

# Precision matrix
if(pF>0)bSigma0 <- chol2inv(chol(diag(B0, pF))) # Prior covariance matrix fixed effects
if(pR>0)gSigma0 <- chol2inv(chol(diag(G0, pR))) # Prior covariance matrix random effects

# Mean parameters
b0 <- 0   # Fixed effects
g0 <- 0   # Random effects
bMu0 <- rep(b0, pF)   # Vector of mean parameters fixed effects
gMu0 <- rep(g0, pR)   # Vector of mean parameters random effects
nu0 <- pR # Degrees of freedom covariance matrix

PriorAlpha <- rep(1e-2, Q)  # Prior alpha Dirichlet - multivariate Bernoulli model

# MCMC parameters 
StartVals <- c(1,2)
nBurn <- 1e4  # Burnin iterations
nIt <- 5e4    # Posterior draws
nChain <- 2   # Number of chains
nSim <- 1e3   # Number of simulated datasets
nThin <- 10   # Thinning rate

nDgm <- 1
MeasurementLevels <- c("Continuous") # Measurementlevels
Methods <- c("Value", "Empirical", "MvB") # Estimation/transformation methods
Populations <- c("Trial", "Intra_Lo") # (sub)populations
TypesSpace <- as.data.frame(t(expand.grid(MeasurementLevels, Methods, Populations, stringsAsFactors = FALSE)))

Types <- t(TypesSpace)
colnames(Types) <- rownames(TypesSpace) <- c("MeasurementLevels", "Methods", "Populations")

# Seeds for simulation
set.seed(2022)
Seeds <- array(sample(1:1e4, (nDgm + 1) * (length(Sigma) + 2) * (nrow(SampleSizes) + 2), replace = FALSE),
                dim = c(nDgm + 1, length(Sigma) + 2, nrow(SampleSizes) + 1))

Rules <- c("Any", "All", "Compensatory")		# Name rules
TestSide <- c("Comp_rs", "All_rs", "Any_rs")		# Names of performed tests
DecisionSelect_RS <- c("DecisionAny.RS", "DecisionAll.RS", "DecisionCompensatory.RS") # Names of evaluated decisions


Weights <- c(0.50,0.50)		# Weights compensatory rule
Alpha <- 0.05			# Type I error rate

SaveRDS <- TRUE # Save simulated data?
