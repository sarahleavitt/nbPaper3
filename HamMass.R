#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program finds the effect estimates from naive Bayes for Hamburg and MA
################################################################################

setwd("~/Boston University/Dissertation/nbPaper3")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(devtools)
load_all("../nbTransmission")


#Function to find true ORs
findORs <- function(df, covariates, outcome, l = 1){
  df <- as.data.frame(df)
  est <- NULL
  #Looping through all covariates
  for(i in 1:length(covariates)){
    #Extracting covariate name
    var <- covariates[i]
    #Creating a table with proportions in each level of the covariate from training data
    tab <- prop.table(table(df[, var], df[, outcome]) + l, 2)
    odds <- tab[, 2] / tab[, 1]
    orMean <- odds/odds[1]
    num <- table(df[, var], df[, outcome]) + l
    se <- NA
    for(i in 1:(nrow(num) - 1)){
      numTab <- num[c(1, i+1), ]
      se <- c(se, sqrt(sum(1/numTab)))
    }
    orCILB = exp(log(orMean) - 1.96 * se)
    orCIUB = exp(log(orMean) + 1.96 * se)
    level <- paste(var, names(orMean), sep = ":")
    cTemp <- cbind.data.frame(level, orMean, orCILB, orCIUB, label = outcome,
                              stringsAsFactors = FALSE)
    est <- bind_rows(est, cTemp)
  }
  return(est)
}


######################### Hamburg Analysis ###########################

#Reading in datasets from HamburgPrep.R
set.seed(103020)
hamInd <- readRDS("../Datasets/HamburgInd.rds")
hamPair <- readRDS("../Datasets/HamburgPair.rds")

orderedHam <- (hamPair
               %>% filter(!is.na(IsolationDiff) & IsolationDiff >= 0)
               %>% mutate(snpClose = ifelse(snpDist < 2, TRUE,
                                     ifelse(snpDist > 12, FALSE, NA)))
)

covarHam <- c("Sex", "Age", "Nationality", "Study", "SmearPos", "HIV",
                "SubstanceAbuse", "Residence", "Milieu", "TimeCat")

resHamG <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                           pairIDVar = "edgeID", goldStdVar = "snpClose",
                           covariates = covarHam, label = "snpCloseNB",
                           l = 1, n = 10, m = 1, nReps = 20)

resHamC <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                           pairIDVar = "edgeID", goldStdVar = "SameGroup",
                           covariates = covarHam, label = "SameGroupNB",
                           l = 1, n = 10, m = 1, nReps = 20)

resHamCovG <- resHamG$probabilities %>% full_join(orderedHam, by = "edgeID")
resHamCovC <- resHamC$probabilities %>% full_join(orderedHam, by = "edgeID")
estHamG <- resHamG$estimates
estHamC <- resHamC$estimates

#Finding the covariate and genetic ORs
trueHamG <- findORs(orderedHam, covariates = covarHam, outcome = "snpClose", l = 1)
trueHamC <- findORs(orderedHam, covariates = covarHam, outcome = "SameGroup", l = 1)

estHamG$level <- factor(estHamG$level, levels = rev(estHamG$level))
estHamC$level <- factor(estHamC$level, levels = rev(estHamC$level))
trueHamG$level <- factor(trueHamG$level, levels = rev(trueHamG$level))
trueHamC$level <- factor(trueHamC$level, levels = rev(trueHamC$level))

estHam <- bind_rows(estHamC, estHamG, trueHamG, trueHamC)

ggplot(data = estHam, aes(x = level, y = orMean, ymin = orCILB,
                          ymax = orCIUB, color = label)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 360)) +
  scale_y_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5,
                           1, 2, 4, 10, 30, 100, 500)) +
  coord_flip()




######################### Massachusetts Analysis ############################

#Reading in cleaned datasets from MassPrep.R
massInd <- readRDS("../Datasets/MassInd.rds")
massPair <- readRDS("../Datasets/MassPair.rds")


#Creating an ordered dataset that also removes pairs with different lineages
orderedMass <- massPair %>% filter(CombinedDiff >= 0, Lineage == "Same" | is.na(Lineage))


#Estimating the probabilities with time difference
covarMass <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                "SharedResG", "GENType", "TimeCat")

resMass <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairIDVar = "EdgeID",
                           goldStdVar = "ContactTrain", covariates = covarMass,
                           label = "ContactTrainNB", l = 0.5, n = 10, m = 1, nReps = 20)

resMassCov <- orderedMass %>% full_join(resMass$probabilities, by = "EdgeID")
estMassC <- resMass$estimates
estMassC$level <- factor(estMassC$level, levels = rev(estMass$level))

#Finding the covariate and genetic ORs
trueMassC <- findORs(orderedMass, covariates = covarMass, outcome = "ContactTrain", l = 0.5)
trueMassC$level <- factor(trueMassC$level, levels = rev(trueMassC$level))

estMass <- bind_rows(estMassC, trueMassC)

ggplot(data = estMass, aes(x = level, y = orMean, ymin = orCILB,
                          ymax = orCIUB, color = label)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 360)) +
  scale_y_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5,
                           1, 2, 4, 10, 30, 100, 500)) +
  coord_flip()


#### Saving results ####
saveRDS(estHam, "../Datasets/HamburgEst.rds")
saveRDS(estMass, "../Datasets/MassEst.rds")





############################# Covariate Tables ###############################

#### Massachusetts ####
covarM <- CreateTableOne(vars = covarMass, factorVars = covarMass, data = orderedMass)
covarM <- as.data.frame(print(covarM, showAllLevels = TRUE))

#Stratified by contact group
covarCM <- CreateTableOne(vars = covarMass, factorVars = covarMass,
                             data = orderedMass, strata = "ContactTrain", test = FALSE)
covarCM <- as.data.frame(print(covarCM, showAllLevels = TRUE))

covarAllM <- cbind.data.frame(covarM, covarCM)


#### Hamburg ####
covarH <- CreateTableOne(vars = covarHam, factorVars = covarHam, data = orderedHam)
covarH <- as.data.frame(print(covarH, showAllLevels = TRUE))

#Stratified by contact group
covarCH <- CreateTableOne(vars = covarHam, factorVars = covarHam,
                          data = orderedHam, strata = "SameGroup", test = FALSE)
covarCH <- as.data.frame(print(covarCH, showAllLevels = TRUE))

covarGH <- CreateTableOne(vars = covarHam, factorVars = covarHam,
                          data = orderedHam, strata = "snpClose", test = FALSE)
covarGH <- as.data.frame(print(covarGH, showAllLevels = TRUE))

covarAllH <- cbind.data.frame(covarH, covarCH, covarGH)


