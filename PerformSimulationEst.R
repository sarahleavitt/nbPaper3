#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program runs the full set of simulations by calling SimRun.R
################################################################################

rm(list = ls())
options(scipen=999)

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(caret)
library(gtools)
library(pROC)
library(devtools)
library(roxygen2)
library(stringr)


#### Batch Mode ####

setwd("/project/sv-thesis/nbPaper3/")
#Getting sample size from the arguements
args <- commandArgs(trailingOnly = TRUE)
sampleSize <- as.numeric(args[1])
#Setting the observation date
if(as.numeric(args[2]) == 1){
  observationDate = "infectionDate"
  dateID = "ID"
}else if(as.numeric(args[2]) == 2){
  observationDate = "sampleDate"
  dateID = "SD"
}
#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#The number of simulations per split
nSim <- 50



#### Interactive Mode ####

# setwd("~/Boston University/Dissertation/nbPaper3")
# iTask <- 1
# sampleSize <- 300
# observationDate <- "infectionDate"
# dateID <- "ID"
# nSim <- 5


#### Getting access to functions in other programs ####

load_all("../nbTransmission")
source("../nbSimulation/SimOutbreak.R")
source("../nbSimulation/SimulateOutbreakS.R")
source("../nbSimulation/SimEvaluate.R")
source("SimCovariatesNBEffect.R")
source("SimRunEst.R")



#### Setting the oubreak parameters ####

print(sampleSize)
print(observationDate)

#Effective population size times generation time (1 day or 1/365 years)
neg <- 0.25
pi <- 1
#Sets reproductive number to off.r if off.p = 0.5
off.r <- 1.2
off.p <- 0.5
#From Yincheng's work
w.shape <- 1.05
w.scale <- 1 / (0.0014 * 365)
w.shift <- 0.25
ws.shape <- w.shape
ws.scale <- w.scale
ws.shift <- w.shift
time <- 20
#Do you want multiple outbreaks
multOutbreaks <- TRUE 

#Genotype parameters
rootseq <- NULL
length <- 300
#This gives the equivalent of 0.5 mutations/genome/year
rate <- 0.5 / length
#Thresholds for SNP distance
thresholds <- c(2, 12)

#Parameters used to manually run functions
covariates <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "timeCat")
pTraining <- 0.6
goldStandard <- "transmission"
truth <- "transmission"
pVar <- "pScaled"



####################### Running the Simulations ######################

coeff <- NULL
performance <- NULL

#Running the simuliaton for the outbreaks
for (iteration in 1:nSim){
  
  #Setting the seed for each run so I can re-run specific iterations with errors
  set.seed(iTask * 1000 + iteration)
  
  tryCatch({
    
    #Running the simulation
    res <- simRunEst()
    
    runID <- paste(iTask, iteration, sampleSize, sep = "_")
    cTemp <- res[[1]] %>% mutate(runID = runID)
    pTemp <- res[[2]] %>% mutate(runID = runID)
    
    coeff <- bind_rows(coeff, cTemp)
    performance <- bind_rows(performance, pTemp)
    #Printing message that the run is finished
    print(paste0("Completed run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
    
  }, error = function(e){
    #Printing message that the run is finished
    print(paste0("Error in run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
    cat("ERROR :", conditionMessage(e), "\n")})
}

#Saving dataframes with summary of results for all of the runs
saveRDS(si, file=paste0("../Simulation_Results/coefficients", "_", iTask, ".rds"))
saveRDS(performance, file=paste0("../Simulation_Results/performance", "_", iTask, ".rds"))


