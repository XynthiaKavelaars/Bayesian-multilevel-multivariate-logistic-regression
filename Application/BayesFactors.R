gDraw_m6 <- lapply(1:(nIt/nThin), function(i){
  if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["gDrawPG"]][[i]]}})


gjDraw_m6 <- lapply(1:(nIt/nThin), function(i){
  if(length(Random_m6) > 0){Pars_m6[["Pars"]][[1]][["gjDrawPG"]][[i]]}})

Tau_BMMLR <-  lapply(1:(nIt/nThin), function(i) {
  Pars_m6[["Pars"]][[1]][["tauDrawPG"]][[i]]}) 


Pars_BMLR <- lapply(1:(nIt/nThin), function(i) {rbind(
    Pars_NH_m6[["Pars"]][[1]][["bDrawPG"]][[i]],
    Pars_NH_m6[["Pars"]][[1]][["gDrawPG"]][[i]])}) 


#### BMMLR ####
Average_Tau_BMMLR <- lapply(1:3, function(i){
  Reduce("+", lapply(Tau_BMMLR, "[[", i)) / length(Tau_BMMLR)
})

Average_gDraw_m6 <- 
  Reduce("+", gDraw_m6) / length(gDraw_m6)



#### Bayes factor BMLR vs BMMLR ####
Estimates.BMLR <- colMeans(do.call(rbind, lapply(Pars_BMLR, as.vector))[,1:12])
names(Estimates.BMLR) <- paste0("b", rep(1:3, each = 4), rep(1:4, times=3))
n.BMLR <- sum(nJ)
Sigma.BMLR <- cor(do.call(rbind, lapply(Pars_BMLR, as.vector))[,1:12])
Hyp.BMLR <- c("b13=b14=b23=b24=b33=b34=0")
BF_BMLR_BMB.BFpack <- BFpack::BF(x = Estimates.BMLR, Sigma = Sigma.BMLR, n = n.BMLR, hypothesis = Hyp.BMLR)
LogBF_BMLR_BMB.BFpack <- log(BF_BMLR_BMB.BFpack[["BFtu_confirmatory"]][1])

#### Bayes factor on individual regression coefficients ####
# compute estimates of random effects
post.mean.ranef <- Reduce("+", lapply(1:(nIt/nThin), function(i) {do.call(abind::abind, c(lapply(1:J, function(j) gjDraw_m6[[i]][,,j] - gDraw_m6[[i]]), along=3))})) / length(gjDraw_m6)

# posterior covariance matrix of random effects
post.covm.ranef <- lapply(1:length(Random_m6), function(p){
  lapply(1:(Q-1), function(q){
    cov(do.call(rbind,lapply(1:(nIt/nThin), function(i) {gjDraw_m6[[i]][p,q,]})))
})})

# matrix of contrasts
Con <- do.call(rbind,lapply(1:(J-1),function(c){
  vec1 <- rep(0,J)
  vec1[c+0:1] <- c(1,-1)
  vec1
}))

# Compute posterior pairwise differences and covariance matrix
post.mean.ranef.con <- lapply(1:length(Random_m6), function(p){lapply(1:(Q-1), function(q){Con %*% post.mean.ranef[p,q,]})})
post.covm.ranef.con <- lapply(1:length(Random_m6), function(p){lapply(1:(Q-1), function(q){Con %*% post.covm.ranef[[p]][[q]] %*% t(Con)})})

# Compute prior pairwise differences and covariance matrix
prior.mean.ranef_meanvar <- rep(0,J)
prior.mean.ranef.con_meanvar <- Con %*% prior.mean.ranef_meanvar

# Compute Savage-Dickey density ratio
ranef.var_mean <- vector("list", Q-1)
prior.covm.ranef_meanvar <- prior.covm.ranef.con_meanvar <- logBF01_meanvar <- lapply(1:length(Random_m6), function(p) vector("list", Q-1))
for(q in 1:(Q-1)){
  ranef.var_mean[[q]] <- diag(Average_Tau_BMMLR[[q]])
  for(p in 1:length(Random_m6)){
     prior.covm.ranef_meanvar[[p]][[q]] <- diag(J) * ranef.var_mean[[q]][p]
     prior.covm.ranef.con_meanvar[[p]][[q]] <- Con %*% prior.covm.ranef_meanvar[[p]][[q]] %*% t(Con)
     logBF01_meanvar[[p]][[q]] <- 
       mvtnorm::dmvnorm(rep(0,J-1),mean=post.mean.ranef.con[[p]][[q]],sigma=post.covm.ranef.con[[p]][[q]],log=TRUE) -
         mvtnorm::dmvnorm(rep(0,J-1),mean=prior.mean.ranef.con_meanvar,sigma=prior.covm.ranef.con_meanvar[[p]][[q]],log=TRUE)
  }}

# Print table
BF10_Mean <- do.call(rbind, lapply(logBF01_meanvar,unlist))
BF10_Mean.df <- rbind.data.frame(
  c("", paste0("$q = ", 1:(Q-1), "$")), 
  cbind.data.frame(c("NIHSS", "NIHSS $\\times$ Trt"),
                   do.call(cbind, lapply(as.data.frame(BF10_Mean), sprintf, fmt="%.3f"))))

print(xtable(BF10_Mean.df, 
             label = "tab:BF", 
             caption = "Logarithmic transformations of Bayes factors of BMLR vs. BMMLR",
             align = paste0("ll", paste0(rep("r", Q-1), collapse=""))), 
      sanitize.text.function = identity, 
      include.rownames = FALSE, include.colnames = FALSE, booktabs = TRUE,
      caption.placement = "top", hline.after = c(-1,1,nrow(x)))

      