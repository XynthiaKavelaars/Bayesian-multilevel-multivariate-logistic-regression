rm(list=ls())
#### Load packages ####
packages <- c('MASS', 'tcltk', 'pgdraw', 'MCMCpack', 'abind', 'foreach', 'doParallel', 'msm', 'mcmcse', 'xtable')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }})

#### Prepare parallel computing ####
nCores <- detectCores()
cl <- makeCluster(nCores -1)
registerDoParallel(cl)

#### Set working directory ####
try(setwd(''))

wd <- getwd()
#### Load functions ####
source("Functions/VariableDefinitions.R")
source("Functions/FunctionsDataGeneration.R")
source("Functions/FunctionsDataTransformation.R")
source("Functions/FunctionsSimulation.R")
source("Functions/DataGeneratingMechanisms.R")

#### Run simulation ####
source("Simulate.R")
source("Evaluate.R")
source("TextualOutput.R")
