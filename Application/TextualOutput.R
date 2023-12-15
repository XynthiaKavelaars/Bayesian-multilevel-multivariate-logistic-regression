load("Workspaces/Application3.RData")

# Convergence
Convergence <- paste0("Multivariate Gelman-Rubin statistic = ", round(Pars.3[["ConvergencePG"]], 3), collapse = "")

# Sample size
SampleSizes <- list(paste0("Sample size trial = ", nrow(Dataset), collapse = ""),
                    paste0("Sample size subset = ", nrow(Data), collapse = ""),
                    paste0("Sample size experimental condition (Hep + Asp)  = ", nE, collapse = ""),
                    paste0("Sample size control condition (Asp)  = ", nC, collapse = ""))

# Means and standard deviations of blood pressure
SummaryData <- list(paste0("Average blood pressure in subset = ", round(meanBp, 2), " (SD = ", round(sdBp, 2), ")", collapse = ""),
                    paste0("Average blood pressure in experimental treatment group (Hep + Asp) = ", round(meanBp_E, 2), " (SD = ", round(sdBp_E, 2), ")", collapse = ""),
                    paste0("Average blood pressure in control treatment group (Asp) = ", round(meanBp_C, 2), " (SD = ", round(sdBp_C, 2), ")", collapse = ""))


Out <- list(SampleSizes = SampleSizes,
            SummaryData = SummaryData,
            Convergence = Convergence)
Out
