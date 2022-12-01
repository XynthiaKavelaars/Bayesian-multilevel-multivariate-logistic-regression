#### 2. Evaluate ####
# Select conditions - multilevel model
if(any(grepl("x", c(Fixed,Random)))){iTypes <- which(!(Types[,"Methods"] %in% c("MvB")))
}else{
  iTypes <- which(Types[,"Methods"] %in% c("Value", "MvB"))
}

# Select conditions - single-level model
if(any(grepl("x", c(Fixed_NH,Random_NH)))){iTypes_NH <- which(!(Types[,"Methods"] %in% c("MvB")))
}else{
  iTypes_NH <- which(Types[,"Methods"] %in% c("Value"))
}

# Select condition - multivariate Bernoulli model
iTypes_MB <- which(Types[,"Methods"] %in% c("MvB"))

# Load workspaces and select decisions
  Decisions <- Decisions1 <- Decisions2 <- lapply(1:nDgm, function(dgm) lapply(1:length(Sigma), function(s) vector("list", nrow(SampleSizes))))
  for(dgm in 1:nDgm){
    for(s in 1:length(Sigma)){
      for(ss in 1:nrow(SampleSizes)){
         Decisions1 <- foreach(sim = 1:nSim, .packages = packages, .verbose = TRUE)%dopar%{
          H <- tryCatch({readRDS(file = paste0(wd, "/Workspaces/Workspaces_H/Decisions/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Decisions", sim, ".RDS", collapse = ""))},
        error = function(e){
        return(NA)})
          NH <- tryCatch({readRDS(file = paste0(wd, "/Workspaces/Workspaces_NH/Decisions/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Decisions_NH", sim, ".RDS", collapse = ""))},
        error = function(e){
        return(NA)})
          MB <- tryCatch({readRDS(file = paste0(wd, "/Workspaces/Workspaces_MB/Decisions/Dgm", dgm, "_var", s, "_clus", SampleSizes[ss,"nCluster"], "_obs", SampleSizes[ss,"nObs"], "_Decisions", sim, ".RDS", collapse = ""))},
                         error = function(e){
                           return(NA)})
          
          list(H=H,NH=NH,MB=MB)}
         Decisions2[["H"]] <- lapply(Decisions1[which(unlist(sapply(Decisions1, function(sim) any(!is.na(sim[["H"]])))))], "[[", "H")
         Decisions2[["NH"]] <- lapply(Decisions1[which(unlist(sapply(Decisions1, function(sim) any(!is.na(sim[["NH"]])))))], "[[", "NH")
         Decisions2[["MB"]] <- lapply(Decisions1[which(unlist(sapply(Decisions1, function(sim) any(!is.na(sim[["MB"]])))))], "[[", "MB")
         Decisions[[dgm]][[s]][[ss]] <- lapply(1:min(length(Decisions2[["H"]]), length(Decisions2[["NH"]]), length(Decisions2[["MB"]])), function(sim){
           list(H = Decisions2[["H"]][[sim]], NH = Decisions2[["NH"]][[sim]], MB = Decisions2[["MB"]][[sim]])
           })
      }
      
      }
  }

# Rejection probabilities - average treatment effects
pReject.ATE_RS <- lapply(1:nDgm, function(dgm){ 
  lapply(1:length(Sigma), function(s){
    lapply(1:nrow(SampleSizes), function(ss){
      ar <- aperm(do.call(abind, c(sapply(seq_along(DecisionSelect_RS), function(d){
        
        # Proportion of rejected samples - multilevel model
        pH <- rowMeans(do.call(cbind, lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) {
          tryCatch({ 
            
            do.call(rbind, 
                    sapply(which(Types[iTypes,"Populations"] == "Trial"), function(i){
                      Decisions[[dgm]][[s]][[ss]][[sim]][["H"]][[i]][["Decision"]][[DecisionSelect_RS[d]]]
                    }, simplify = FALSE, USE.NAMES = TRUE))
          },
          error = function(e){
            message(paste0(e, " in sim ", sim, collapse = ""))
          })
        })))
        
        # Proportion of rejected samples - single-level model
        pNH <- tryCatch({
          rowMeans(do.call(cbind, lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) {
                      do.call(rbind, 
                    sapply(which(Types[iTypes_NH,"Populations"] == "Trial"), function(i){
                      Decisions[[dgm]][[s]][[ss]][[sim]][["NH"]][[i]][["Decision"]][[DecisionSelect_RS[d]]]
                    }, simplify = FALSE, USE.NAMES = TRUE))
          }
        )))},
        error = function(e){
          message(paste0(e, " in sim "#, sim
                         , collapse = ""))
        })
        
        # Proportion of rejected samples - multivariate Bernoulli model
        pMB <- rowMeans(do.call(cbind, lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) {
          tryCatch({ 
            
            do.call(rbind, 
                    sapply(which(Types[iTypes_MB,"Populations"] == "Trial"), function(i){
                      Decisions[[dgm]][[s]][[ss]][[sim]][["MB"]][[i]][["Decision"]][[DecisionSelect_RS[d]]]
                    }, simplify = FALSE, USE.NAMES = TRUE))
          },
          error = function(e){
            message(paste0(e, " in sim ", sim, collapse = ""))
          })
        })))
        
        # Standard error pReject
        seH <- sqrt(pH * (1-pH) / length(Decisions[[dgm]][[s]][[ss]]))
        seNH <- sqrt(pNH * (1-pNH) / length(Decisions[[dgm]][[s]][[ss]]))
        seMB <- sqrt(pMB * (1-pMB) / length(Decisions[[dgm]][[s]][[ss]]))
        
        res <- cbind(c(pH, pNH, pMB
                       ), c(seH, seNH, seMB
                         ))
        return(res)
      }, simplify = FALSE, USE.NAMES = TRUE), along = 3)), c(1,3,2))
      dimnames(ar) <- list(c(paste0(Types[intersect(iTypes, which(Types[,"Populations"] == "Trial")),"Methods"], "_H")
                             ,paste0(Types[intersect(iTypes_NH, which(Types[,"Populations"] == "Trial")),"Methods"], "_NH")
                             ,paste0(Types[intersect(iTypes_MB, which(Types[,"Populations"] == "Trial")),"Methods"], "_MB")
                             ),
                           DecisionSelect_RS,
                           c("p", "se"))
      ar
    })
  })
})



# Rejection probabilities - conditional treatment effects
pReject.CTE_RS <- lapply(1:nDgm, function(dgm){ 
  lapply(1:length(Sigma), function(s){
    lapply(1:nrow(SampleSizes), function(ss){
      ar <- aperm(do.call(abind, c(sapply(seq_along(DecisionSelect_RS), function(d){
        
        # Proportion of rejected samples - multilevel model
        pH <- 
          tryCatch({ 
            rowMeans(
            do.call(cbind, lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) {
            do.call(rbind, 
                    sapply(which(Types[iTypes,"Populations"] == "Intra_Lo"), function(i){
                      Decisions[[dgm]][[s]][[ss]][[sim]][["H"]][[i]][["Decision"]][[DecisionSelect_RS[d]]]
                    }, simplify = FALSE, USE.NAMES = TRUE))
          
        }))
            )
             },
        error = function(e){
          message(paste0(e, " in sim ", sim, collapse = ""))
        })
        
        # Proportion of rejected samples - single-level model
        pNH <- 
          rowMeans(do.call(cbind, lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) {
            tryCatch({   do.call(rbind, 
                    sapply(which(Types[iTypes_NH,"Populations"] == "Intra_Lo"), function(i){
                      Decisions[[dgm]][[s]][[ss]][[sim]][["NH"]][[i]][["Decision"]][[DecisionSelect_RS[d]]]
                    }, simplify = FALSE, USE.NAMES = TRUE))
         
          },
          error = function(e){
            message(paste0(e, " in sim ", sim
                           , collapse = ""))
          })
            })))
        
        # Proportion of rejected samples - multivariate Bernoulli model
        pMB <- rowMeans(do.call(cbind, lapply(1:length(Decisions[[dgm]][[s]][[ss]]), function(sim) {
          tryCatch({ 
            
            do.call(rbind, 
                    sapply(which(Types[iTypes_MB,"Populations"] == "Intra_Lo"), function(i){
                      Decisions[[dgm]][[s]][[ss]][[sim]][["MB"]][[i]][["Decision"]][[DecisionSelect_RS[d]]]
                    }, simplify = FALSE, USE.NAMES = TRUE))
          },
          error = function(e){
            message(paste0(e, " in sim ", sim, collapse = ""))
          })
        })))
        
        # Standard error pReject
        seH <- sqrt(pH * (1-pH) / length(Decisions[[dgm]][[s]][[ss]]))
        seNH <- sqrt(pNH * (1-pNH) / length(Decisions[[dgm]][[s]][[ss]]))
        seMB <- sqrt(pMB * (1-pMB) / length(Decisions[[dgm]][[s]][[ss]]))
        
        res <- cbind(c(pH, pNH, pMB
        ), c(seH, seNH, seMB
        ))
        return(res)
      }, simplify = FALSE, USE.NAMES = TRUE), along = 3)), c(1,3,2))
      dimnames(ar) <- list(c(paste0(Types[intersect(iTypes, which(Types[,"Populations"] == "Intra_Lo")),"Methods"], "_H")
                             ,paste0(Types[intersect(iTypes_NH, which(Types[,"Populations"] == "Intra_Lo")),"Methods"], "_NH")
                             ,paste0(Types[intersect(iTypes_MB, which(Types[,"Populations"] == "Intra_Lo")),"Methods"], "_MB")
      ),
      DecisionSelect_RS,
      c("p", "se"))
      ar
    })
  })
})

pReject.CTE_RS


# Tables
Effects <- c("ATEs", "CTEs (range)", "CTEs (value)")
Captions <- c("Proportions of superiority decisions (P) and standard errors (SE)
                    by data-generating mechanism, estimation method, and decision rule.")
names(Captions) <- c("ATE")

#### Table average treatment effect ####
Ind.ATE_H <- which(rownames(pReject.ATE_RS[[1]][[1]][[1]]) %in% c("Empirical_H", "Analytical_H"))
Ind.ATE_NH <- which(rownames(pReject.ATE_RS[[1]][[1]][[1]]) %in% c("Empirical_NH", "Analytical_NH"))
Ind.ATE_MVB <- which(grepl("MB", rownames(pReject.ATE_RS[[1]][[1]][[1]])))
names(Ind.ATE_H) <- paste0(Types[which(Types[,"Methods"] %in% c("Analytical","Empirical") & Types[,"Populations"] == "Trial"),"Methods"], "_H")
names(Ind.ATE_NH) <- paste0(Types[which(Types[,"Methods"] %in% c("Analytical","Empirical") & Types[,"Populations"] == "Trial"),"Methods"], "_NH")
names(Ind.ATE_MVB) <- "MVB"
Ind.ATE <- c(Ind.ATE_H, Ind.ATE_NH, Ind.ATE_MVB)
names(Ind.ATE) <- c(names(Ind.ATE_H), names(Ind.ATE_NH), names(Ind.ATE_MVB))

# Extract proportions per decision rule
pReject.ATE.Any_RS <- lapply(1:nDgm, function(dgm){lapply(1:length(Sigma), function(s){lapply(1:nrow(SampleSizes), function(ss){
  pReject.ATE_RS[[dgm]][[s]][[ss]][,"DecisionAny.RS",]
})})})

pReject.ATE.All_RS <- lapply(1:nDgm, function(dgm){lapply(1:length(Sigma), function(s){lapply(1:nrow(SampleSizes), function(ss){
  pReject.ATE_RS[[dgm]][[s]][[ss]][,"DecisionAll.RS",]
})})})

pReject.ATE.Comp_RS <- lapply(1:nDgm, function(dgm){lapply(1:length(Sigma), function(s){lapply(1:nrow(SampleSizes), function(ss){
  pReject.ATE_RS[[dgm]][[s]][[ss]][,"DecisionCompensatory.RS",]
})})})

pReject.ATE <- lapply(1:nDgm, function(dgm){
  lapply(1:length(Sigma), function(s){
  lapply(1:nrow(SampleSizes), function(ss){
    x <- cbind.data.frame(pReject.ATE.Any_RS[[dgm]][[s]][[ss]], 
                            pReject.ATE.All_RS[[dgm]][[s]][[ss]], 
                            pReject.ATE.Comp_RS[[dgm]][[s]][[ss]])
    y <- x[Ind.ATE,,drop=FALSE]
    rownames(y) <- names(Ind.ATE)
    y
    })})})
pReject.ATE

#### Latex ####
Truth.ATE <- unlist(lapply(1:nDgm, function(dgm) sapply(1:length(Sigma), function(s) {
  MV <- paste0("$\\bm{\\delta} = (", paste0(sprintf("%1.3f", Truth[[dgm]][[s]][["eATE"]][["Delta"]]), collapse = ", "), ")$", collapse = "")
  W <- paste0("$\\bm{\\delta} (\\bm{w}) = ", sprintf("%1.3f", Truth[[dgm]][[s]][["eATE"]][["DeltaW"]]), "$", collapse = "")
  return(paste(MV, W, sep = ", "))
})))

pReject.ATE.P <- do.call(rbind, lapply(1:nDgm, function(dgm){
  do.call(rbind, lapply(1:length(Sigma), function(s){
    do.call(rbind, pReject.ATE[[dgm]][[s]])}))
}))

pReject.ATE.P1 <- matrix(NA, nrow = nrow(pReject.ATE.P), ncol = length(Rules) + ncol(pReject.ATE.P))
pReject.ATE.P1[,-seq(1,ncol(pReject.ATE.P1),3)] <- as.matrix(pReject.ATE.P)

pReject.ATE.Print <- pReject.ATE.P1
pReject.ATE.Print[,seq(2,ncol(pReject.ATE.P1),3)] <- apply(pReject.ATE.P1[,seq(2,ncol(pReject.ATE.P1),3)], 2, function(x) paste0("$", sprintf("%.3f", x), "$"))
pReject.ATE.Print[,seq(3,ncol(pReject.ATE.P1),3)] <- apply(pReject.ATE.P1[,seq(3,ncol(pReject.ATE.P1),3)], 2, function(x) paste0("(", sprintf("%.3f", x), ")"))

tabPRejectATE <- cbind.data.frame(rep(c("BMMLR", "BMLR", "BMB"), times = nrow(SampleSizes)), pReject.ATE.Print)



#### Table conditional treatment effect - range####
Ind.CTE.Range_H <- which(rownames(pReject.CTE_RS[[1]][[1]][[1]]) %in% c("Empirical_H", "Analytical_H"))
Ind.CTE.Range_NH <- which(rownames(pReject.CTE_RS[[1]][[1]][[1]]) %in% c("Empirical_NH", "Analytical_NH"))
Ind.CTE.Range_MVB <- which(grepl("MB", rownames(pReject.CTE_RS[[1]][[1]][[1]])))
names(Ind.CTE.Range_H) <- paste0(Types[which(Types[,"Methods"] %in% c("Analytical","Empirical") & Types[,"Populations"] == "Intra_Lo"),"Methods"], "_H")
names(Ind.CTE.Range_NH) <- paste0(Types[which(Types[,"Methods"] %in% c("Analytical","Empirical") & Types[,"Populations"] == "Intra_Lo"),"Methods"], "_NH")
names(Ind.CTE.Range_MVB) <- "MVB"
Ind.CTE.Range <- c(Ind.CTE.Range_H, Ind.CTE.Range_NH, Ind.CTE.Range_MVB)
names(Ind.CTE.Range) <- c(names(Ind.CTE.Range_H), names(Ind.CTE.Range_NH), names(Ind.CTE.Range_MVB))


# Extract proportions per decision rule
pReject.CTE.Any_RS <- lapply(1:nDgm, function(dgm){lapply(1:length(Sigma), function(s){lapply(1:nrow(SampleSizes), function(ss){
  pReject.CTE_RS[[dgm]][[s]][[ss]][,"DecisionAny.RS",]
})})})

pReject.CTE.All_RS <- lapply(1:nDgm, function(dgm){lapply(1:length(Sigma), function(s){lapply(1:nrow(SampleSizes), function(ss){
  pReject.CTE_RS[[dgm]][[s]][[ss]][,"DecisionAll.RS",]
})})})

pReject.CTE.Comp_RS <- lapply(1:nDgm, function(dgm){lapply(1:length(Sigma), function(s){lapply(1:nrow(SampleSizes), function(ss){
  pReject.CTE_RS[[dgm]][[s]][[ss]][,"DecisionCompensatory.RS",]
})})})

pReject.CTE.Range <- lapply(1:nDgm, function(dgm){
  lapply(1:length(Sigma), function(s){
    lapply(1:nrow(SampleSizes), function(ss){
      x <- cbind.data.frame(pReject.CTE.Any_RS[[dgm]][[s]][[ss]], 
                            pReject.CTE.All_RS[[dgm]][[s]][[ss]], 
                            pReject.CTE.Comp_RS[[dgm]][[s]][[ss]])
      y <- x[Ind.CTE.Range,,drop=FALSE]
      rownames(y) <- names(Ind.CTE.Range)
      y
    })})})
pReject.CTE.Range

Truth.CTE.Range <- unlist(lapply(1:nDgm, function(dgm) sapply(1:length(Sigma), function(s) {
  MV <- paste0("$\\bm{\\delta} = (", paste0(sprintf("%1.3f", Truth[[dgm]][[s]][["eCTE"]][["Delta"]]), collapse = ", "), ")$", collapse = "")
  W <- paste0("$\\bm{\\delta} (\\bm{w}) = ", sprintf("%1.3f", Truth[[dgm]][[s]][["eCTE"]][["DeltaW"]]), "$", collapse = "")
return(paste(MV, W, sep = ", "))
  })))

#### Latex ####
pReject.CTE.Range.P <- do.call(rbind, lapply(1:nDgm, function(dgm){
  do.call(rbind, lapply(1:length(Sigma), function(s){
    do.call(rbind, pReject.CTE.Range[[dgm]][[s]])}))
}))

pReject.CTE.Range.P1 <- matrix(NA, nrow = nrow(pReject.CTE.Range.P), ncol = length(Rules) + ncol(pReject.CTE.Range.P))
pReject.CTE.Range.P1[,-seq(1,ncol(pReject.CTE.Range.P1),3)] <- as.matrix(pReject.CTE.Range.P)

pReject.CTE.Range.Print <- pReject.CTE.Range.P1
pReject.CTE.Range.Print[,seq(2,ncol(pReject.CTE.Range.P1),3)] <- apply(pReject.CTE.Range.P1[,seq(2,ncol(pReject.CTE.Range.P1),3)], 2, function(x) paste0("$", sprintf("%.3f", x), "$"))
pReject.CTE.Range.Print[,seq(3,ncol(pReject.CTE.Range.P1),3)] <- apply(pReject.CTE.Range.P1[,seq(3,ncol(pReject.CTE.Range.P1),3)], 2, function(x) paste0("(", sprintf("%.3f", x), ")"))

tabPRejectCTE.Range <- cbind.data.frame(rep(c("BMMLR", "BMLR", "BMB"), times = nrow(SampleSizes)), pReject.CTE.Range.Print)

#### pReject combined table ATE + CTE ####
tabPReject <- rbind(tabPRejectATE, tabPRejectCTE.Range)

AddToRow.pReject <- list()
AddToRow.pReject$pos <- as.list(seq(0,nDgm * nrow(SampleSizes) * length(c(Ind.ATE, Ind.CTE.Range)) -1, length(Ind.CTE.Range)))

Command_dgm <- paste0(unlist(rbind.data.frame(
  matrix(paste0(c(" ", "\\midrule \n "), "\\multicolumn{", ncol(tabPReject), "}{c}{", 
                c("\\textbf{Average treatment effect: }", "\\textbf{Conditional treatment effect: }"),
                c(Truth.ATE, Truth.CTE.Range), "} \\\\ \n \\midrule \n", 
                paste0("& & \\multicolumn{2}{l}{", Rules, "}", collapse = ""),"\\\\\n"),
               ncol = 2), 
  matrix(rep(" ", nDgm * 2 * (nrow(SampleSizes) - 1)), ncol = 2))))
Command_SampleSize <- unlist(matrix(rep(paste0(c("", rep("\\midrule \n ", nrow(SampleSizes)-1)), 
                                        paste0(" $J$ = ", SampleSizes[,"nCluster"], ", $n_{j}$ = ", SampleSizes[,"nObs"]))
                                 
                                 , 2), ncol = 2 * nrow(SampleSizes)))
Command_pars <- unlist(rbind.data.frame(matrix(paste0(rep(paste0(paste0(" & ", rep(c(" & \\multicolumn{1}{l}{p} & \\multicolumn{1}{l}{SE}"), 
                           times = length(Rules)), collapse = ""), " \\\\\n ", collapse = ""), 2), c("", "\\midrule \n")), ncol = 2),
                           matrix(rep(paste0("& \\multicolumn{", ncol(tabPReject) - 1, "}{l}{ } \\\\ \n \\midrule \n", collapse = ""), 6), ncol = 2)))




AddToRow.pReject$command <- paste0(Command_dgm,
                                         Command_SampleSize, Command_pars)
                              

Align.pReject <- paste0("ll", paste(rep("p{0.02cm}rr", length(Rules)), collapse = ""))

print(xtable(tabPReject, 
             caption= Captions["ATE"],
             label="tab:pReject", align=Align.pReject,
             digits=c(1,1,rep(3,ncol(tabPReject)-1))),
      add.to.row=AddToRow.pReject, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)
