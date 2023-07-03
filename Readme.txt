This folder contains all necessary files to replicate the simulation study for the manuscript ”Bayesian multivariate logistic regression for superiority and inferiority decision-making under treatment heterogeneity." by Xynthia Kavelaars, Joris Mulder, and Maurits Kaptein.


########################################################################################################
#					Simulation 						       #
########################################################################################################
Description of files and folders:
\Simulation:				Folder containing subfolders with scripts, functions, simulation conditions.
RunSimulation.R				R-script to source all functions, run simulation, and create output in tables
Simulate.R					R-script to sample parameters via Gibbs sampler, transform them and make decisions
Evaluate.R					R-script to extract results from output and tabulate them
TextualOutput.R				R-script to create in-text results

\\Functions:				Folder containing functions used by scripts in folder \Simulation:
DataGeneratingMechanisms.R		Parameters to define data generating mechanisms
FunctionsSimulation.R			Various functions to estimate parameters 
FunctionsDataGeneration.R		Various functions to generate data
FunctionsDataTransformation.R		Various functions to perform transformations of data and parameters
VariableDefinitions.R			Variable definitions

\\Workspaces:				Folder where folders with workspaces will be stored after running "RunSimulation.R"

\\\Workspaces_H:				Folder where results of Bayesian multilevel multivariate logistic regression (BMMLR) will be saved
\\\\Bias:					Empty folder where bias of regression parameters will be saved per simulated dataset
\\\\Data:					Empty folder where simulated datasets will be saved
\\\\Decisions:				Empty folder where analysis results will be saved per simulated dataset
\\\\Diags:					Empty folder where MCMC diagnostics will be saved per simulated dataset
\\\\Parameters:				Empty folder where regression parameters will be saved per simulated dataset

\\\Workspaces_NH:				Folder where results of single-level Bayesian multivariate logistic regression (BMLR) will be saved
\\\\Bias:					Empty folder where bias of regression parameters will be saved per simulated dataset
\\\\Decisions:				Empty folder where analysis results will be saved per simulated dataset
\\\\Diags:					Empty folder where MCMC diagnostics will be saved per simulated dataset
\\\\Parameters:				Empty folder where regression parameters will be saved per simulated dataset

\\\Workspaces_MB:				Folder where results of multivariate Bernoulli model (BMB) will be saved
\\\\Decisions:				Empty folder where analysis results will be saved per simulated dataset
\\\\Theta:					Empty folder where success probabilities will be saved per simulated dataset

Instructions:
1. Set path to working directory to folder "\Simulation" in call "setwd('')" in "RunSimulation.R"
2. Execute "RunSimulation.R"

########################################################################################################
#					Application 						       #
########################################################################################################
Description of files and folders:
\Application:				Folder containing subfolders with scripts, functions, and dataset of the IST.
RunApplication.R				R-script to run application
PrepareData.R				R-script to prepare data for analysis
Evaluate.R					R-script to sample parameters via Gibbs sampler, transform them and make superiority decisions
Tabulate.R					R-script to create Latex tables and plots presented in manuscript
TextualOutput.R				R-script to create in-text results
VariableDefinitions.R			Various variable definitions
Figures.R					R-script to create figures

\\Functions:
FunctionsApplication.R			Various functions to estimate parameters

\\Workspaces:				Empty folder, where workspaces will be stored after running "RunApplication.R"
\\Plots:					Empty folder, where workspaces will be stored after running "RunApplication.R"

Instructions:
1. Set path to working directory to main directory "Research archive" in call "setwd('')" in "RunApplication.R"
2. Download data from Third International Stroke Trial from https://doi.org/http://doi.org/10.7488/ds/1350 to folder "\Application"
3. Name the dataset "DatasetIST3.sas7bdat".
4. Execute "RunApplication.R"
########################################################################################################
#					Further information					       #
########################################################################################################
These scripts are based on parallel computation. The default number of cores is set at the number of cores of the PC minus one. 
To change the number of cores, change parameter(s) "nCores" in "\Simulation\RunSimulation.R" and "\Application\RunApplication.R".

For any help with the files in this archive, please contact Xynthia Kavelaars (x.m.kavelaars@tilburguniversity.edu). 
