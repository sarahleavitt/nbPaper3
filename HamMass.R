#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program finds the effect estimates from naive Bayes for Hamburg and MA
################################################################################

setwd("~/Boston University/Dissertation")
options(scipen = 999)
#rm(list = ls())

library(dplyr)
library(tidyr)
library(devtools)
library(ggplot2)
library(tableone)
load_all("nbTransmission")


#Function to find true ORs
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
    logorMean <- log(odds/odds[1])
    num <- table(df[, var], df[, outcome]) + l
    logorSE <- NA
    for(i in 1:(nrow(num) - 1)){
      numTab <- num[c(1, i+1), ]
      logorSE <- c(logorSE, sqrt(sum(1/numTab)))
    }
    logorCILB = logorMean - 1.96 * logorSE
    logorCIUB = logorMean + 1.96 * logorSE
    level <- paste(var, names(logorMean), sep = ":")
    estTemp <- cbind.data.frame(level, logorMean, logorSE, logorCILB, logorCIUB,
                                outcome = outcome, stringsAsFactors = FALSE)
    est <- bind_rows(est, estTemp)
  }
  return(est)
}

######################### Hamburg Analysis ###########################

#Reading in datasets from HamburgAnalysis.R
set.seed(103020)
hamInd <- readRDS("Datasets/HamburgInd.rds")
hamPair <- readRDS("Datasets/HamburgPair.rds")

orderedHam <- (hamPair
               %>% filter(!is.na(IsolationDiff) & IsolationDiff > 0)
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

estHam <- (estHamG
           %>% bind_rows(estHamC)
           %>% mutate(Value = gsub("[A-z0-9]+\\:", "", level),
                      Variable = gsub("\\:[A-z0-9+-<=>]+", "", level))
           %>% arrange(abs(logorMean))
           %>% mutate(labelc = ifelse(label == "SameGroupNB", "Confirmed contact",
                               ifelse(label == "snpCloseNB", "Close genetic relatedness", NA)),
                      Variable = factor(Variable, levels = covarHam, 
                                        labels = c("Sex", "Age Group",
                                                   "Nationality", "City",
                                                   "Smear Result", "HIV Status",
                                                   "Substance Abuse", "Residence",
                                                   "Affiliation with\nlocal drinking scene",
                                                   "Observation time\ndifference")),
                      
                      Value = factor(Value, levels = c("m-f", "f-m", "m-m",
                                                       "Same", "Same-Other", "Diff-Other",
                                                       "Same-Germany", "InfectorSmear+", "InfectorHIV+",
                                                       "Both", "Neither", "BothStable",
                                                       "BothHomeless", "BothNot", "BothAssociated",
                                                       ">4y", "3-4y", "2-3y", "1-2y",
                                                       "f-f", "Different", "Diff-Germany",
                                                       "InfectorSmear-", "InfectorHIV-", "<=1y"),
                                     
                                     labels = c("Male to female", "Female to male", "Male to Male",
                                                "Same", "Same foreign\ncountry",  "Different foreign\ncountries",
                                                "Both German", "Infector smear+", "Infector HIV+",
                                                "Both", "Neither", "Both permanent",
                                                "Both homeless", "Neither affiliated", "Both affiliated",
                                                ">4 years", "3-4 years", "2-3 years", "1-2 years",
                                                "Female to female", "Different",
                                                "One German, one\nforeign country", "Infector smear-",
                                                "Infector HIV-", "<1 year")))
)

#Saving results
saveRDS(estHam, "Datasets/HamburgEst.rds")

ggplot(data = estHam, aes(x = Value, y = exp(logorMean), ymin = exp(logorCILB),
                           ymax = exp(logorCIUB), color = labelc)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  facet_wrap(~Variable, scales = "free_y") +
  ylab("Odds ratio with 95% confidence interval") +
  labs(color = "Training Links") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 360),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "bottom") +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 500)) +
  coord_flip() +
  ggsave("Figures/HamCovar.png", width = 10, height = 6, dpi = 300)





######################### Massachusetts Analysis ############################

#Reading in cleaned datasets from MassPrep.R
massInd <- readRDS("Datasets/MassInd.rds")
massPair <- readRDS("Datasets/MassPair.rds")


#Creating an ordered dataset that also removes pairs with different lineages
orderedMass <- massPair %>% filter(CombinedDiff >= 0, Lineage == "Same" | is.na(Lineage))


#Estimating the probabilities with time difference
covarMass <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                "SharedResG", "GENType", "TimeCat")

set.seed(103020)
resMass <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairIDVar = "EdgeID",
                           goldStdVar = "ContactTrain", covariates = covarMass,
                           label = "ContactTrainNB", l = 0.5, n = 10, m = 1, nReps = 20)

resMassCov <- orderedMass %>% full_join(resMass$probabilities, by = "EdgeID")

estMass <- (estMass
            %>% mutate(Level = factor(level, levels = rev(level)),
                       Value = gsub("[A-z0-9]+\\:", "", level),
                       Variable = gsub("\\:[A-z0-9+-<>=]+", "", level))
            %>% arrange(abs(logorMean))
            %>% mutate(labelc = ifelse(label == "ContactTrainNB", "Confirmed contact", NA),
                       Variable = factor(Variable, levels = covarMass, 
                                         labels = c("Sex", "Age Group",
                                                    "Country of Birth", "MA County of Residence",
                                                    "Smear Result", "Immune-suppressed",
                                                    "Shared resistance to\nhow many drugs",
                                                    "CDC GENType",
                                                    "Observation time\ndifference")),
                       
                       Value = factor(Value, levels = c("m-f", "f-m", "m-m",
                                                        "Same", "Same-Other", "Diff-Other",
                                                        "Same-USA", "Neighbor", "InfectorSmear+", 
                                                        "InfectorImmuneSup+", "3+", "2", "1",
                                                        ">4y", "3-4y", "2-3y", "1-2y",
                                                        "f-f", "Different", "Diff-USA",
                                                        "Other", "InfectorSmear-", "InfectorImmuneSup-",
                                                        "0", "<=1y"),
                                      
                                      labels = c("Male to female", "Female to male", "Male to Male",
                                                 "Same", "Same foreign\ncountry",  "Different foreign\ncountries",
                                                 "Both US born", "Neighboring", "Infector smear+",
                                                 "Infector not\nsuppressed", "3+", "2", "1",
                                                 ">4 years", "3-4 years", "2-3 years", "1-2 years",
                                                 "Female to female", "Different", "One US, one\nforeign country",
                                                 "More distant", "Infector smear-", "Infector\nsuppressed",
                                                 "0", "<1 year")))
)

#Saving results
saveRDS(estMass, "Datasets/MassEst.rds")


ggplot(data = estMass, aes(x = Value, y = exp(logorMean), ymin = exp(logorCILB),
                          ymax = exp(logorCIUB), color = labelc)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  facet_wrap(~Variable, scales = "free_y") +
  ylab("Odds ratio with 95% confidence interval") +
  labs(color = "Training Links") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 360),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "bottom") +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 500)) +
  coord_flip() +
  ggsave("Figures/MassCovar.png", width = 8, height = 5, dpi = 300)


## PRESENTATION VERSION ##
ggplot(data = estMass, aes(x = Value, y = exp(logorMean), ymin = exp(logorCILB),
                           ymax = exp(logorCIUB), color = label)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  facet_wrap(~Variable, scales = "free_y") +
  ylab("Odds ratio with 95% confidence interval") +
  theme_bw(base_size = 16) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 360),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "none") +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 500)) +
  coord_flip() +
  ggsave("Figures/MassCovar_pres.png", width = 11, height = 7.5, dpi = 300)






############################# Covariate Tables ###############################


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

covarAllH <- cbind.data.frame(covarH, covarGH, covarCH)



#### Massachusetts ####

covarM <- CreateTableOne(vars = covarMass, factorVars = covarMass, data = orderedMass)
covarM <- as.data.frame(print(covarM, showAllLevels = TRUE))

#Stratified by contact group
covarCM <- CreateTableOne(vars = covarMass, factorVars = covarMass,
                             data = orderedMass, strata = "ContactTrain", test = FALSE)
covarCM <- as.data.frame(print(covarCM, showAllLevels = TRUE))

covarAllM <- cbind.data.frame(covarM, covarCM)





