#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program makes tables and figures for the simulation results of the 
# contribution of covariates paper
################################################################################

#rm(list = ls())
options(scipen = 999)
setwd("~/Boston University/Dissertation/Simulation_Results")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tableone)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)

est <- readRDS("estimates_1.rds")
perform <- readRDS("performance_1.rds")




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
snpORs <- (est
            %>% filter(outcome == "snpClose")
            %>% select(runID, level, orMeanS = orMean, orCILBS = orCILB,
                       orCIUBS = orCIUB)
)

estORs <- (est
           %>% filter(outcome != "transmission")
           %>% full_join(trueOR, by = "level")
           %>% mutate(biasTruth = orMean - orTruth,
                      coverage = orTruth <= orCIUB & orTruth >= orCILB)
)

#Plot of bias
ggplot(data = estORs %>% filter(outcome != "snpClose"),
       aes(x = level, y = biasTruth, fill = outcome, color = outcome)) +
  geom_violin() +
  geom_hline(yintercept = 0) +
  facet_wrap(~pTraining, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


#Summary of coverage
coverage <- (estORs
             %>% group_by(outcome, level, pTraining)
             %>% summarize(n = n(),
                           pCoverage = sum(coverage)/n)
)

ggplot(data = coverage %>% filter(outcome != "snpClose"),
       aes(x = level, y = pCoverage, fill = outcome, color = outcome)) +
  geom_jitter() +
  facet_wrap(~pTraining, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

