#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program makes tables and figures for the simulation results of the 
# contribution of covariates paper
################################################################################

#rm(list = ls())
options(scipen = 999)
setwd("~/Boston University/Dissertation/Simulation_ResultsEE_1.30.20")

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
             %>% select(runID, goldStd, threshold = pTraining, aucVal,
                        pCorrect, pTop5, pTop10, pTop25, pTop50)
             %>% gather(metric, value, -goldStd, -threshold, -runID)
             %>% mutate(metric = factor(metric, levels = c("aucVal", "pCorrect", "pTop5",
                                                           "pTop10", "pTop25", "pTop50"),
                                        labels = c("Area under the ROC",
                                                   "Proportion Correct",
                                                   "Proportion in Top 5%",
                                                   "Proportion in Top 10%",
                                                   "Proportion in Top 25%",
                                                   "Proportion in Top 50%")))
)

ggplot(data = longData, aes(x = factor(threshold), y = value,
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
           %>% summarize(logorTruth = mean(logorMean))
)

#Calculating bias and coverage
estORsAll <- (est
              %>% filter(outcome != "transmission")
              %>% full_join(trueOR, by = "level")
              %>% mutate(logorBias = logorMean - logorTruth,
                         coverage = logorTruth <= logorCIUB & logorTruth >= logorCILB)
              %>% filter(!outcome %in% c("snpCloseGS", "transmissionNB"),
                         !is.na(coverage))
              %>% mutate(Estimate = ifelse(outcome == "snpClose", "Conventional SNP Distance",
                                           "Naive Bayes SNP Distance"))
)

#Removing all but threshold 2 for main analysis
estORs <- estORsAll %>% filter(threshold == 2)

#Summary of MAPE and coverage
biasCov <- (estORs
             %>% group_by(outcome, threshold, level)
             %>% summarize(n = n(),
                           logorTruth = first(logorTruth),
                           mse = mean(logorBias, na.rm = TRUE) ^ 2,
                           mape = mean(abs(logorBias/logorTruth)),
                           pCoverage = sum(coverage)/n,
                           width = mean(logorCIUB - logorCILB))
            %>% filter(!outcome %in% c("snpCloseGS", "transmissionNB"))
            %>% mutate(Estimate = ifelse(outcome == "snpClose", "Conventional SNP Distance",
                                         "Naive Bayes SNP Distance"))
)


#Plot of bias
ggplot(data = estORs) +
  geom_boxplot(aes(x = level, y = logorMean, fill = Estimate, color = Estimate),
               alpha = 0.5) +
  geom_point(data = trueOR, aes(x = level, y = logorTruth, shape = "True Log Odds Ratio")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("Covariate and Level") +
  ylab("Mean Estimated Log Odds Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggsave("../Figures/EE_Bias.png", width = 7, height = 5, dpi = 300)


#Plot of MAPE
pe1 <- ggplot(data = biasCov, aes(x = level, y = mape, fill = Estimate,
                           color = Estimate, shape = Estimate)) +
  geom_point() +
  xlab("Covariate and Level") +
  ylab("Mean Absolute Precentage Error") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none")

#Plot of MSE
pe2 <- ggplot(data = biasCov, aes(x = level, y = mse, fill = Estimate,
                           color = Estimate, shape = Estimate)) +
  geom_point() +
  xlab("Covariate and Level") +
  ylab("Mean Standard Error") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none")


#Plot of CI coverage
pe3 <- ggplot(data = biasCov, aes(x = level, y = pCoverage, fill = Estimate,
                           color = Estimate, shape = Estimate)) +
  geom_point() +
  geom_hline(yintercept = 0.95) +
  xlab("Covariate and Level") +
  ylab("Confidence Interval Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none")


#Plot of CI width
pe4 <- ggplot(data = biasCov, aes(x = level, y = width, fill = Estimate,
                           color = Estimate, shape = Estimate)) +
  geom_point() +
  xlab("Covariate and Level") +
  ylab("Confidence Interval Width") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        legend.title = element_blank())


grid.arrange(pe1, pe2, pe3, pe4)
pe <- arrangeGrob(pe1, pe2, pe3, pe4)
ggsave(file = "../Figures/EE_Error.png", plot = pe,
       width = 8, height = 6, units = "in", dpi = 300)




#### Supplementary Figure: Bias by Threshold ####

ggplot(data = estORsAll) +
  facet_wrap(~threshold) +
  geom_boxplot(aes(x = level, y = logorMean, fill = outcome, color = outcome),
               alpha = 0.5) +
  geom_point(data = trueOR, aes(x = level, y = logorTruth)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

