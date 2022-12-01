rm(list=ls())

#### Load packages ####
packages <- c("haven", "MCMCpack", "coda", "pgdraw", "tcltk", "foreach", "doParallel", "mcmcse", "xtable") 

package.check <- lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }})

#### Prepare parallel computing ####
nCores <- detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)

#### Set working directory####
setwd('')

#### Load functions ####
source("Simulation/Functions/FunctionsSimulation.R")
source("Simulation/Functions/FunctionsDataTransformation.R")
source("Application/VariableDefinitions.R")
Data <- read_sas("Application/DatasetIST3.sas7bdat")
set.seed(2022) 

wd <- getwd()
#### Run application ####
source("Application/PrepareData.R")
source("Application/Evaluate.R")
source("Application/Tabulate.R")
source("Application/TextualOutput.R")
