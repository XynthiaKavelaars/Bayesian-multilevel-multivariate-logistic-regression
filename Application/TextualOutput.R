# Sample sizes 
nTotal <- sum(nJ)
nTrt <- table(do.call(rbind,DataC)[,"trt"])

# Cluster sizes
clusterSize <- summary(nJ[which(nJ>0)])
clusterSD <- sd(nJ[which(nJ>0)])

# Summary of covariate 
xChars.E <- summary(do.call(rbind, lapply(1:sum(nJ > 0), function(j) DataC[[j]][DataC[[j]][,"trt"]==1,]))[,"x"])
xSD.E <- sd(unlist(do.call(rbind, lapply(1:sum(nJ > 0), function(j) DataC[[j]][DataC[[j]][,"trt"]==1,]))[,"x"]))
xChars.C <- summary(do.call(rbind, lapply(1:sum(nJ > 0), function(j) DataC[[j]][DataC[[j]][,"trt"]==0,]))[,"x"])
xSD.C <- sd(unlist(do.call(rbind, lapply(1:sum(nJ > 0), function(j) DataC[[j]][DataC[[j]][,"trt"]==0,]))[,"x"]))
xMean <- mean(unlist(do.call(rbind, lapply(1:sum(nJ > 0), function(j) DataC[[j]]))[,"x"]))
xSD <- sd(unlist(do.call(rbind, lapply(1:sum(nJ > 0), function(j) DataC[[j]]))[,"x"]))

xRange <- range(do.call(rbind,DataC)$x)
RangesApp
ValuesApp

# Autocorrelation
AC_H <- lapply(Diags_m6[["AC"]], function(x) apply(x[c(1,10),], 1, range))
AC_NH <- lapply(Diags_NH_m6[["AC"]], function(x) apply(x[c(1,10),], 1, range))

# Multivariate scale reduction factor
Convergence_H <- Diags_m6[["Convergence"]]
Convergence_NH <- Diags_NH_m6[["Convergence"]]

# Log Bayes factor BMLR - BMB
LogBF_BMLR_BMB.BFpack <- log(BF_BMLR_BMB.BFpack[["BFtu_confirmatory"]][1])
