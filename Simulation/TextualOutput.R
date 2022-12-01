#### Textual output ####
  Diags <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) vector("list", nrow(SampleSizes))))
  for(dgm in 1:nDgm){
    for(s in 1:length(Sigma)){
      for(ss in 1:nrow(SampleSizes)){
        Diags[[dgm]][[s]][[ss]] <- foreach(sim = 1:nSim, .packages = packages, .verbose = TRUE)%dopar%{
          H <-  readRDS(file = paste0(wd, "/Workspaces/Workspaces_H/Diags/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Diags", sim, ".RDS", collapse = ""))
          NH <- readRDS(file = paste0(wd, "/Workspaces/Workspaces_NH/Diags/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Diags_NH", sim, ".RDS", collapse = ""))
          return(list(H=H, NH=NH))
        }
      }
    }}


#### Autocorrelation ####
AC_H_gDraw_l1 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) range(do.call(rbind, lapply(1:nrow(SampleSizes), function(ss) Reduce("+", lapply(1:nSim, function(sim){
  range(Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["AC"]][["gDraw"]][1,])
}))/nSim)))))
AC_H_gDraw_l5 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) range(do.call(rbind, lapply(1:nrow(SampleSizes), function(ss) Reduce("+", lapply(1:nSim, function(sim){
  range(Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["AC"]][["gDraw"]][5,])
}))/nSim)))))

AC_H_bDraw_l1 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) range(do.call(rbind, lapply(1:nrow(SampleSizes), function(ss) Reduce("+", lapply(1:nSim, function(sim){
  range(Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["AC"]][["bDraw"]][1,])
}))/nSim)))))
AC_H_bDraw_l5 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) range(do.call(rbind, lapply(1:nrow(SampleSizes), function(ss) Reduce("+", lapply(1:nSim, function(sim){
  range(Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["AC"]][["bDraw"]][5,])
}))/nSim)))))

AC_H_tDraw_l1 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) range(do.call(rbind, lapply(1:nrow(SampleSizes), function(ss) Reduce("+", lapply(1:nSim, function(sim){
  range(Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["AC"]][["tDraw"]][1,])
}))/nSim)))))
AC_H_tDraw_l5 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) range(do.call(rbind, lapply(1:nrow(SampleSizes), function(ss) Reduce("+", lapply(1:nSim, function(sim){
  range(Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["AC"]][["tDraw"]][5,])
}))/nSim)))))

rangeAC_l1 <- range(c(unlist(AC_H_bDraw_l1), unlist(AC_H_gDraw_l1), unlist(AC_H_tDraw_l1)))

cat(paste0(c("Minimum", "Maximum"), " autocorrelation lag 1 in regression coefficients: ", round(range(unlist(rangeAC_l1)), 3), "\n"))
cat(paste0("Maximum autocorrelation lag 1 in "), c("Fixed rc", "Random rc", "Random variance matrix")[which(c(max(rangeAC_l1) %in% unlist(AC_H_bDraw_l1), 
        max(rangeAC_l1) %in% unlist(AC_H_gDraw_l1),
        max(rangeAC_l1) %in% unlist(AC_H_tDraw_l1)))])

rangeAC_l5 <- range(c(unlist(AC_H_bDraw_l5), unlist(AC_H_gDraw_l5), unlist(AC_H_tDraw_l5)))
cat(paste0(c("Minimum", "Maximum"), " autocorrelation lag 5 in regression coefficients: ", round(range(unlist(rangeAC_l5)), 3), "\n"))
cat(paste0("Maximum autocorrelation lag 5 in "), c("Fixed rc", "Random rc", "Random variance matrix")[which(c(max(rangeAC_l5) %in% unlist(AC_H_bDraw_l5), 
                                                                                                            max(rangeAC_l5) %in% unlist(AC_H_gDraw_l5),
                                                                                                              max(rangeAC_l5) %in% unlist(AC_H_tDraw_l5)))])
#### Convergence ####
Convergence.bDraw <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) do.call(rbind, lapply(1:nrow(SampleSizes), function(ss){
  max(sapply(1:nSim, function(sim){
  Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["Convergence"]][["bDraw"]]
}))}
))))

Convergence.gDraw <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) do.call(rbind, lapply(1:nrow(SampleSizes), function(ss){
  max(sapply(1:nSim, function(sim){
    Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["Convergence"]][["gDraw"]]
  }))}
))))

Convergence.tDraw <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) do.call(rbind, lapply(1:nrow(SampleSizes), function(ss){
  max(sapply(1:nSim, function(sim){
    Diags[[dgm]][[s]][[ss]][[sim]][["H"]][["Convergence"]][["tDraw"]]
  }))}
))))

nonConverged <- sum(c(unlist(Convergence.bDraw), unlist(Convergence.gDraw), unlist(Convergence.tDraw)) > 1.10)
cat(paste0(nonConverged, " simulations had a multivariate scale reduction factor > 1.10"))


#### Bias treatment effects
BiasWeighted <- lapply(1:nDgm, function(dgm)
  lapply(1:length(Sigma), function(s) 
    lapply(1:nrow(SampleSizes), function(ss){ 
      H <- Reduce("+", lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) do.call(rbind, lapply(1:length(iTypes), function(i) Decisions[[dgm]][[s]][[ss]][[sim]][["H"]][[i]][["Bias"]][["BiasWeighted"]])))) / length(Decisions[[dgm]][[s]][[ss]])
      NH <- Reduce("+", lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) do.call(rbind, lapply(1:length(iTypes_NH), function(i) Decisions[[dgm]][[s]][[ss]][[sim]][["NH"]][[i]][["Bias"]][["BiasWeighted"]])))) / length(Decisions[[dgm]][[s]][[ss]])
      MB <- Reduce("+", lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) do.call(rbind, lapply(1:length(iTypes_MB), function(i) Decisions[[dgm]][[s]][[ss]][[sim]][["MB"]][[i]][["Bias"]][["BiasWeighted"]])))) / length(Decisions[[dgm]][[s]][[ss]])
      list(H=H, NH=NH,MB=MB)
    })))

rangeBiasWeighted <- lapply(1:nDgm, function(dgm)
  lapply(1:length(Sigma), function(s) 
    lapply(1:nrow(SampleSizes), function(ss){
      range(BiasWeighted[[dgm]][[s]][[ss]])})))

cat(paste0(c("Minimum", "Maximum"), " bias in weighted delta: ", round(range(unlist(rangeBiasWeighted)), 3), "\n"))

BiasMV <- lapply(1:nDgm, function(dgm) 
  lapply(1:length(Sigma), function(s) 
    lapply(1:nrow(SampleSizes), function(ss){
      H <- Reduce("+", lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) do.call(rbind, lapply(1:length(iTypes), function(i) Decisions[[dgm]][[s]][[ss]][[sim]][["H"]][[i]][["Bias"]][["BiasMultivariate"]])))) / length(Decisions[[dgm]][[s]][[ss]])
      NH <- Reduce("+", lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) do.call(rbind, lapply(1:length(iTypes_NH), function(i) Decisions[[dgm]][[s]][[ss]][[sim]][["NH"]][[i]][["Bias"]][["BiasMultivariate"]])))) / length(Decisions[[dgm]][[s]][[ss]])
      MB <- Reduce("+", lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) do.call(rbind, lapply(1:length(iTypes_MB), function(i) Decisions[[dgm]][[s]][[ss]][[sim]][["MB"]][[i]][["Bias"]][["BiasMultivariate"]])))) / length(Decisions[[dgm]][[s]][[ss]])
      list(H=H, NH=NH,MB=MB)
    })))

rangeBiasMV <- lapply(1:nDgm, function(dgm)
  lapply(1:length(Sigma), function(s) 
    lapply(1:nrow(SampleSizes), function(ss){
      range(BiasMV[[dgm]][[s]][[ss]])})))

cat(paste0(c("Minimum", "Maximum"), " bias in delta: ", round(range(unlist(rangeBiasMV)), 3), "\n"))

