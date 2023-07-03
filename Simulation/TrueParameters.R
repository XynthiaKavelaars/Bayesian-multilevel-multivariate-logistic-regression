gammas_BMMLR <- round(rbind(TruePars[["b"]], TruePars[["g"]]), 3)
gammas_BMLR <- round(rbind(TruePars_NH[["b"]],TruePars_NH[["g"]]), 3)

header <- c("", paste0("$q_{", 1:Q, "}$"))
gammas <- rbind(header, cbind(paste0("$p_{", 0:3, "}$"), format(gammas_BMLR, nsmall=3)))


print(xtable(gammas, caption = "True regression parameters used for numerical evaluation", label = "tab:EffectSize",
             align = paste0("ll", paste0(rep("r", Q), collapse=""))),
     caption.placement="top", hline.after=c(-1,1,nrow(gammas)), include.rownames=FALSE, include.colnames=FALSE,
     sanitize.text.function=identity,
     math.style.negative=TRUE, booktabs=TRUE)

