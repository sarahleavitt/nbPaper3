#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program finds the effect estimates from naive Bayes for Hamburg and MA
################################################################################

setwd("~/Boston University/Dissertation")
#rm(list = ls())

library(dplyr)
library(tidyr)
library(devtools)
load_all("nbTransmission")


######################### Hamburg Analysis ###########################

#Reading in datasets from HamburgPrep.R
set.seed(103020)
hamInd <- readRDS("Datasets/HamburgInd.rds")
hamPair <- readRDS("Datasets/HamburgPair.rds")

orderedHam <- (hamPair
               %>% filter(!is.na(IsolationDiff) & IsolationDiff >= 0)
               %>% mutate(snpClose = ifelse(snpDist < 2, TRUE,
                                     ifelse(snpDist > 12, FALSE, NA)))
)

covarHam <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV",
                "SubstanceAbuse", "Residence", "Milieu", "TimeCat")

resHamG <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                           pairIDVar = "edgeID", goldStdVar = "snpClose",
                           covariates = covarHam, label = "HamSNPs",
                           l = 1, n = 10, m = 1, nReps = 20)

resHamC <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                           pairIDVar = "edgeID", goldStdVar = "SameGroup",
                           covariates = covarHam, label = "HamContacts",
                           l = 1, n = 10, m = 1, nReps = 20)

resHamCovG <- resHamG$probabilities %>% full_join(orderedHam, by = "edgeID")
resHamCovC <- resHamC$probabilities %>% full_join(orderedHam, by = "edgeID")
estHamG <- resHamG$estimates
estHamC <- resHamC$estimates

estHamG$level <- factor(estHamG$level, levels = rev(estHamG$level))
estHamC$level <- factor(estHamC$level, levels = rev(estHamC$level))

estHam <- bind_rows(estHamC, estHamG)

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
massInd <- readRDS("Datasets/MassInd.rds")
massPair <- readRDS("Datasets/MassPair.rds")

#Creating an ordered dataset that also removes pairs with different lineages
orderedMass <- massPair %>% filter(CombinedDiff >= 0, Lineage == "Same" | is.na(Lineage))


#Estimating the probabilities with time difference
covarMass <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                "SharedResG", "GENType", "TimeCat")

resMass <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairIDVar = "EdgeID",
                           goldStdVar = "ContactTrain", covariates = covarMass,
                           label = "ContactTime", l = 0.5, n = 10, m = 1, nReps = 20)

resMassCov <- orderedMass %>% full_join(resMass$probabilities, by = "EdgeID")
estMass <- resMass$estimates
estMass$level <- factor(estMass$level, levels = rev(estMass$level))

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
saveRDS(estHam, "Datasets/HamburgEst.rds")
saveRDS(estMass, "Datasets/MassEst.rds")


