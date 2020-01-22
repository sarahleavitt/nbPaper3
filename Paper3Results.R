#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program makes tables and figures for the simulation results of the 
# contribution of covariates paper
################################################################################

#rm(list = ls())
options(scipen = 999)
setwd("~/Boston University/Dissertation/Simulation_ResultsEE_1.14.20")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tableone)
library(gridExtra)

#Initializing dataframes
est <- NULL
perform <- NULL

#Reading in the results
for (file in list.files()){
  
  if(grepl("^est", file)){
    estTemp <- readRDS(file)
    est <- bind_rows(est, estTemp)
  }
  
  if(grepl("^perform", file)){
    pTemp <- readRDS(file)
    perform <- bind_rows(perform, pTemp) 
  }
}



###################### Performance Metrics #####################

longData <- (perform
             %>% select(runID, goldStd, pTraining, aucVal, pCorrect,
                        pTop5, pTop10, pTop25, pTop50)
             %>% gather(metric, value, -goldStd, -pTraining, -runID)
             %>% mutate(metric = factor(metric, levels = c("aucVal", "pCorrect", "pTop5",
                                                           "pTop10", "pTop25", "pTop50"),
                                        labels = c("Area under the ROC",
                                                   "Proportion Correct",
                                                   "Proportion in Top 5%",
                                                   "Proportion in Top 10%",
                                                   "Proportion in Top 25%",
                                                   "Proportion in Top 50%")))
)

ggplot(data = longData, aes(x = factor(pTraining), y = value,
                            fill = goldStd, color = goldStd)) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  scale_x_discrete(name = "Scenario") +
  scale_y_continuous(name = "Value") +
  facet_wrap(~ metric, nrow = 3, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "none")


###################### Contribution of Covariates ########################

#Taking the average OR over all simulations as the 'truth'
trueOR <- (est
           %>% filter(outcome == "transmission")
           %>% group_by(level)
           %>% summarize(orTruth = mean(orMean))
)

estORs <- (est
           %>% filter(outcome != "transmission")
           %>% full_join(trueOR, by = "level")
           %>% mutate(logorMean = log(orMean),
                      logorTruth = log(orTruth),
                      orBias = orMean - orTruth,
                      logorBias = logorMean - logorTruth,
                      coverage = orTruth <= orCIUB & orTruth >= orCILB)
)

#Summary of MAPE and coverage
biasCov <- (estORs
            %>% filter(pTraining == 0.6)
             %>% group_by(outcome, level)
             %>% summarize(n = n(),
                           orTruth = first(orTruth),
                           mse = mean(logorBias, na.rm = TRUE) ^ 2,
                           mape = mean(abs(logorBias/logorTruth)),
                           pCoverage = sum(coverage)/n,
                           width = mean(log(orCIUB) - log(orCILB)))
)

#Plot of bias
ggplot(data = estORs %>% filter(!outcome %in% c("snpClose", "transmissionNB"),
                                !is.na(coverage)),
       aes(x = level, y = logorBias, fill = outcome, color = outcome)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


ggplot(data = estORs %>% filter(!outcome %in% c("snpClose", "transmissionNB"),
                                !is.na(coverage))) +
  geom_boxplot(aes(x = level, y = logorMean, fill = outcome, color = outcome),
               alpha = 0.5) +
  geom_point(data = trueOR, aes(x = level, y = log(orTruth))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


#Plot of MAPE
ggplot(data = biasCov %>% filter(!outcome %in% c("snpClose", "transmissionNB"),
                                 !is.na(pCoverage)),
       aes(x = level, y = mape, fill = outcome, color = outcome, shape = outcome)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

#Plot of MSE
ggplot(data = biasCov %>% filter(!outcome %in% c("snpClose", "transmissionNB"),
                                 !is.na(pCoverage)),
       aes(x = level, y = mse, fill = outcome, color = outcome, shape = outcome)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


#Plot of CI coverage
ggplot(data = biasCov %>% filter(!outcome %in% c("snpClose", "transmissionNB")),
       aes(x = level, y = pCoverage, fill = outcome, color = outcome, shape = outcome)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = 0.95)


#Plot of CI width
ggplot(data = biasCov %>% filter(!outcome %in% c("snpClose", "transmissionNB")),
       aes(x = level, y = width, fill = outcome, color = outcome, shape = outcome)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

