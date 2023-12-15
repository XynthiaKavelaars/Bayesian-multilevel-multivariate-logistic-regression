rm(list=ls())
# Table of variables retrieved from:
# https://trialsjournal.biomedcentral.com/articles/10.1186/1745-6215-12-101/tables/2

# Reference retrieved from:
# https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(97)04011-7/fulltext


# Load packages
packages <- c('MCMCpack', 'coda', 'pgdraw', 'TruncatedNormal',  
              'foreach', 'doParallel', 'xtable', 'extrafont')

package.check <- lapply(packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }})

# Set working directory
try(setwd("C:/Users/xynni/surfdrive/Werk/Bayesian Adaptive Clinical Trials/Project 2/Submitted version - Multivariate Behavioral Research/Revision V1/Research archive/Application"))


Dataset <- read.csv("IST_Data.csv")

# Run application
source("PrepareData.R")
source("Functions/VariableDefinitions.R")
source("Functions/FunctionsApplication.R")
source("SampleParameters.R")
source("MakeTableFigure.R")
source("TextualOutput.R")


