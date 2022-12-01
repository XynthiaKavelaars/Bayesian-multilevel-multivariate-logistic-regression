#### Prepare data ####
Dataset <- Data[,SelectedVars]
Dataset$id <- 1:nrow(Dataset)

# Group IDs as factors
Dataset$hosp_id <- as.numeric(Dataset$randhosp_id)
Dataset$country_id <- as.factor(Dataset$country)
levels(Dataset$country_id) <- 1:length(levels(Dataset$country_id))

# Recode treatment and dependency
Dataset$trt <- (Dataset$treatment == "rt-PA") * 1
Dataset$indep6 <- (Dataset$ohs6 < 3)
Dataset$indep18 <- (Dataset$ohs18 < 3)
Dataset$nostrk7 <- (Dataset$adjudicated == 0) * 1
Dataset$x <- Dataset[[x]]
Dataset$trt_x <- Dataset$trt * Dataset$x # Interaction

# Subset UK-data
DataC <- lapply(unlist(unique(Dataset[grep("UK", Dataset$country),"hosp_id"])), function(j) {
  Dataset[grepl("UK", Dataset$country) & 
            Dataset[,"hosp_id"] == j & 
            apply(Dataset[,xVars,drop=FALSE], 1, function(x) all(!is.na(x))) &
            apply(Dataset[,yVars,drop=FALSE], 1, function(x) all(!is.na(x))),]})

# Sample sizes per cluster
nJ <- sapply(DataC, nrow)
J <- sum(sapply(DataC, function(j) nrow(j)> 0))

# Extract response data of subset
yDataC <- lapply(which(sapply(DataC, function(j) nrow(j)> 0)), function(j){
  y <- Multivariate2Multinomial(as.matrix(DataC[[j]][,yVars,drop=FALSE]))
  colnames(y) <- c("11", "10", "01", "00")
  return(y)})

# Define ranges and values of subpopulations
ValuesApp <- RangesApp <- list()
ValuesApp[["Continuous"]] <- ValueFun(Dataset$x)
RangesApp[["Continuous"]] <- RangeFun(Dataset$x)
