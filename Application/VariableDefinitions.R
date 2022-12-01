#### Variable definitions ####
# Subset of variables 
SelectedVars <- c("randhosp_id", "country", "treatment", "gender",
                  "age", "weight", "nihss", "sbprand",
                  "sbpstart", "dbpstart",
                  "adjudicated", "ohs6", "ohs18", "alivefav18")



#  MCMC parameters 
nChain <- 2 # Number of chains
nIt <- 50e4 # Number of iterations
nBurn <- 1e4 # Number of burnin iterations
nThin <- 1e1 # Thinning rate
StartVals <- c(0.5,1) # Starting values
x <- "nihss" # Name of predictor variable
idVars <- c("country_id") # Identification of country
xVars <- c("trt", "x", "trt_x") # names of cvariates
yVars <- c("nostrk7", "indep6") # Names of response variables
K <- length(yVars) # Number of reponse variables
Q <- 2^K # Number of (multinomial) joint response categories


MeasurementLevels <- c("Continuous") # Measurementlevel of predictor variable
Methods <- c("Value", "Empirical", "MvB") # Methods to compute treatment effects
Populations <- c("Trial", "Intra_Lo", "Intra_MidL", "Intra_MidH", "Intra_Hi") # Names of subpopulations
TypesSpace <- as.data.frame(t(expand.grid(MeasurementLevels, Methods, Populations, stringsAsFactors = FALSE))) # All possible combinations of effects and (sub)populations
Exclude <- lapply(TypesSpace, function(x) # Excluded combinations
  all(c("Intra_MidL", "Value") %in% x) |
    all(c("Intra_MidH", "Value") %in% x) |
    all(c("Discrete", "Analytical", "Trial") %in% x) |
    all(c("Discrete", "Analytical", "Intra_Lo") %in% x) | 
    all(c("Discrete", "Empirical", "Intra_Lo") %in% x) |
    all(c("Discrete", "Analytical", "Intra_Hi") %in% x) | 
    all(c("Discrete", "Empirical", "Intra_Hi") %in% x))


# Conditions
Types <- t(TypesSpace[,which(do.call(rbind, Exclude) == FALSE)])
colnames(Types) <- rownames(TypesSpace) <- c("MeasurementLevels", "Methods", "Populations")

# Function to define ranges and values of subpopulations
RangeFun <- function(x){
  list(Trial = c(-Inf, Inf),
       Intra_Lo = c(0,5),
       Intra_MidL = c(6,14),
       Intra_MidH = c(15,24),
       Intra_Hi = c(25,Inf))}
# Values defining subpopulations
 ValueFun <- function(x){
  list(Trial = mean(x, na.rm=TRUE),
       Intra_Lo = mean(x, na.rm=TRUE) - 1 * sd(x, na.rm = TRUE),
       Intra_Hi = mean(x, na.rm=TRUE) + 1 * sd(x, na.rm = TRUE))}




# Variable definitions
Continuous <- TRUE

PriorAlpha <- rep(0.1,4) # Prior parameters multivariate Bernoulli model

Rules <- c("Any", "All", "Compensatory")		# Name rules
TestSide <- c("Comp_rs", "All_rs", "Any_rs")		# Names of performed tests
DecisionSelect_RS <- c("DecisionAny.RS", "DecisionAll.RS", "DecisionCompensatory.RS") # Name of decision rules and test side


Weights <- c(0.2,0.8)		# Weights compensatory rule
Alpha <- 0.05			# Type I error rate

SaveRDS <- TRUE # Save simulated data
