#### Variable definitions Application ####
Q <- ncol(Y12)


# Sample size
n <- nrow(Data)
nE <- sum(Data$TRT == 1)
nC <- n - nE

# Prior
Alpha0 <- rep(1, ncol(Y12))

# Decision
Alpha <- 0.05

# Gibbs sampling parameters
nIt <- 2e4
nBurn <- 1e4
start <- c(0,0.5,1)

# Prior parameters
b0 <- 0
B0 <- 1e-2

# Decision parameters
weights <- c(0.25,0.75)


MeasurementLevels <- c("Discrete", "Continuous")
Methods <- c("MvB", "Value", "Empirical", "Analytical")
Populations <- c("Trial", "Lo", "Lo3", "Lo2", "Lo1", "Hi1", "Hi2", "Hi3", "Hi")
Scopes <- c("Old")
MethodFit <- c("PG")
TypesSpace <- as.data.frame(t(expand.grid(MeasurementLevels, Methods, Populations, Scopes, stringsAsFactors = FALSE)))
Exclude <- lapply(TypesSpace, function(x) {
    any(c("Discrete") %in% x) |
    all(c("Trial", "Value") %in% x) | 
    all(c("Lo", "Value") %in% x) |
    all(c("Hi", "Value") %in% x) 
    })

Types <- t(as.matrix(TypesSpace[,which(do.call(rbind, Exclude) == FALSE)]))
colnames(Types) <- rownames(TypesSpace) <- c("MeasurementLevels", "Methods", "Populations", "Scopes")

Types.Range <- Types[Types[,"Methods"] %in% c("MvB", "Empirical") & Types[,"Populations"] %in% c("Trial", "Lo", "Hi"),]
Types.Value <- Types[Types[,"Methods"] %in% c("Value"),]

Rules <- c("Any", "All", "Compensatory")

Truths <- list(Delta = c(0,0), DeltaW = 0)

# Conditions
Covariates <- list(RSBP_C = sdBp)
Ranges <- lapply(Covariates, function(x)  
  list(Trial = c(-Inf, Inf),
       Lo = c(-Inf, -1 * x),
       Lo3 = c(-Inf, -2 * x),
       Lo2 = c(-2 * x, -1 * x),
       Lo1 = c(-1 * x, 0),
       Hi1 = c(0, 1 * x),
       Hi2 = c(1 * x, 2 * x),
       Hi3 = c(2 * x, Inf),
       Hi = c(1 * x, Inf)))


Values <- lapply(Covariates, function(x)  
  list(Lo3 = c(-3 * x),
       Lo2 = c(-2 * x),
       Lo1 = c(-1 * x),
       Hi1 = c( 1 * x),
       Hi2 = c( 2 * x),
       Hi3 = c( 3 * x)))

nRange <- sapply(names(Covariates), function(cv){
  do.call(rbind, sapply(Populations, function(pop){
    n <- nrow(Data[Data[,cv] > min(Ranges[[cv]][[pop]]) & Data[,cv] < max(Ranges[[cv]][[pop]]),])
    nE <- nrow(Data[Data[,"TRT"] == 1 & Data[,cv] > min(Ranges[[cv]][[pop]]) & Data[,cv] < max(Ranges[[cv]][[pop]]),])
    nC <- nrow(Data[Data[,"TRT"] == 0 & Data[,cv] > min(Ranges[[cv]][[pop]]) & Data[,cv] < max(Ranges[[cv]][[pop]]),])
    unlist(list(n = n, nE = nE, nC = nC))
  }, simplify = FALSE, USE.NAMES = TRUE))
}, simplify = FALSE, USE.NAMES = TRUE)

names.Rules <- c("Any", "All", "Comp")
range.Populations <- lapply(c("Age", "Delay", "Bp"), function(cv) c(
  paste0("$ - \\infty < \\text{", cv, "} < \\infty$", collapse = ""),
  paste0("$ - \\infty < \\text{", cv, "} < -1 \\text{ SD }$", collapse = ""),
  paste0("$ - \\infty < \\text{", cv, "} < -2 \\text{ SD }$", collapse = ""),
  paste0("$-2 \\text{ SD } < \\text{", cv, "} < -1 \\text{ SD }$", collapse = ""),
                                                       paste0("$-1 \\text{ SD } < \\text{", cv, "} < 0 \\text{ SD }$", collapse = ""), 
                                                       paste0("$0 \\text{ SD } < \\text{", cv, "} < +1 \\text{ SD }$", collapse = ""),
                                                       paste0("$+1 \\text{ SD } < \\text{", cv, "} < +2 \\text{ SD }$", collapse = ""),
                                                       paste0("$+2 \\text{ SD } < \\text{", cv, "} < \\infty$", collapse = ""),
  paste0("$+1 \\text{ SD } < \\text{", cv, "} < \\infty$", collapse = "")))

names.Populations <- paste0(c("ATE: ", rep("CTE: ", length(Ranges[["RSBP_C"]]) - 1)))
