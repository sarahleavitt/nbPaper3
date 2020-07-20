# nbPaper3

This directory contains the code to produce the results for Leavitt et al. ????, 
published in ??? (LINK HERE). It contains code to run the simulations exploring the
relationship between covariate associations with transmission and covariate 
associations with close genetic distance and if the relationship can be improved
by the naive Bayes transmission method. It also finds the contribution of the 
covariates to the transmission probabilities estimated for a TB outbreak in Hamburg,
Germany and TB surveillance data from the Massachusetts department of health.
These programs use the nbTransmission R package available
on CRAN or GitHub (https://github.com/sarahleavitt/nbTransmission)
as well as the simulation programs in the nbSimulation directory
(https://github.com/sarahleavitt/nbSimulation).
 

## Simulation Programs

The following programs were used to run a simulation study to asses how the association
between covariates and close genetic relatedness mirros the true relationship between
covariates and transmission and if the effect estimates from the naive Bayes transmission
method are a better estimate of the truth.

### SimCovariatesNBEffect.R

This program is analogous to SimCovariates.R in the nbSimulation directory but adds
a few extra covariates for this simulation since the emphasis is exploring different
covariate contributions to the probabilities.


### SimRunEst.R

This program contains a function "simRunEst" which performs one iteration of a 
simulation which simulates an outbreak with pathogen genetics and covariates,
estimates the relative transmission probability using naive Bayes, assesses the 
performance, and estimates the contribution of each covariate value to the probabilities 
(in the form of log odds ratios with confidence intervals).


### PerformSimulationEst.R

This program is a wrapper for "simRunEst" which runs the simulation coded in that 
function nSim times. It saves the raw pair-level results for one iteration 
and the effect estimates, reproductive number, and performance for all iterations. 
This is run by "SimQsubEst.qsub" with sample size and observation date as inputs.


### SimQsubEst.qsub

This a qsub file to run PerformSimulationEst.R as a batch job. It runs t parallel jobs
and therefore the total number of simulations is t*nSim times. It needs sample size
and observation date (1 for infection date, 2 for sampling date) as inputs.


### Paper3Results.R

This program reads in all of the results from the the simulation scenarios and then 
calculates all results used in Leavitt et al. ???? including values in the text,
tables and figures.


***


## Application Programs

The following program was used to analyze the Hamburg and TB outbreaks. They
use the data prepared in the HamburgPrep.R program in the nbPaper1 directory and 
the MassPrep.R program in the nbPaper2 directory.


### HamMass.R

This program estimates relative transmission probabilities for a TB outbreak in
Hamburg, Germany and for TB surveillance data from Massachusetts. It then extracts
the contribution of each covariate value from the output and creates figures 
visualizing the effect estimates.


