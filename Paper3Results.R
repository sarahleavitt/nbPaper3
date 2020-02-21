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

#### Supplementary Table: True ORs ####

#Taking the average OR over all simulations as the 'truth'
trueOR <- (est
           %>% filter(outcome == "transmission")
           %>% group_by(level)
           %>% summarize(logorTruth = mean(logorMean))
           %>% mutate(orTruth = exp(logorTruth),
                      level = gsub("Y", "Z", level),
                      level = gsub("timeCat", "Time", level))
)

#Calculating bias and coverage
estORsAll <- (est
              %>% filter(outcome != "transmission")
              %>% mutate(level = gsub("Y", "Z", level),
                         level = gsub("timeCat", "Time", level))
              %>% full_join(trueOR, by = "level")
              %>% mutate(logorBias = logorMean - logorTruth,
                         coverage = logorTruth <= logorCIUB & logorTruth >= logorCILB)
              %>% filter(!outcome %in% c("snpCloseGS", "transmissionNB"),
                         !is.na(coverage))
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
)


#### Figure: Boxplots of bias ####
ggplot(data = estORs) +
  geom_boxplot(aes(x = level, y = logorMean, fill = outcome, color = outcome),
               alpha = 0.5) +
  geom_point(data = trueOR, aes(x = level, y = logorTruth, shape = "True log odds ratio")) +
  scale_shape_discrete(labels = c(expression("True log odds ratio (OR"^"T"*")"))) +
  scale_fill_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                 expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  scale_color_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                 expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("Covariate and Level") +
  ylab("Mean Estimated Log Odds Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggsave("../Figures/EE_Bias.png", width = 7, height = 5, dpi = 300)


## PRESENTATION VERSION ##
ggplot(data = estORs) +
  geom_boxplot(aes(x = level, y = logorMean, fill = outcome, color = outcome),
               alpha = 0.5) +
  geom_point(data = trueOR, aes(x = level, y = logorTruth, shape = "True log odds ratio")) +
  scale_shape_discrete(labels = c(expression("True log odds ratio (OR"^"T"*")"))) +
  scale_fill_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                 expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  scale_color_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                  expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("Covariate and Level") +
  ylab("Mean Estimated Log Odds Ratio") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggsave("../Figures/EE_Bias_pres.png", width = 7, height = 5, dpi = 300)


#Plot of MAPE
pe1 <- ggplot(data = biasCov, aes(x = level, y = mape, fill = outcome,
                           color = outcome, shape = outcome)) +
  geom_point() +
  xlab("Covariate and Level") +
  ylab("Mean Absolute Precentage Error") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none")

#Plot of MSE
pe2 <- ggplot(data = biasCov, aes(x = level, y = mse, fill = outcome,
                           color = outcome, shape = outcome)) +
  geom_point() +
  xlab("Covariate and Level") +
  ylab("Mean Standard Error") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none")


#Plot of CI coverage
pe3 <- ggplot(data = biasCov, aes(x = level, y = pCoverage, fill = outcome,
                           color = outcome, shape = outcome)) +
  geom_point() +
  geom_hline(yintercept = 0.95) +
  xlab("Covariate and Level") +
  ylab("Confidence Interval Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none")


#Plot of CI width
pe4 <- ggplot(data = biasCov, aes(x = level, y = width, fill = outcome,
                           color = outcome, shape = outcome)) +
  geom_point() +
  scale_fill_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                 expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  scale_color_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                  expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  scale_shape_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                  expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  xlab("Covariate and Level") +
  ylab("Confidence Interval Width") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(margin = margin(l = 0)),
        legend.spacing.x = unit(0.1, "cm"),
        legend.position = "bottom",
        legend.title = element_blank())



#### Figure: Plot of Error and Coverage ####

grid.arrange(pe1, pe2, pe3, pe4)
pe <- arrangeGrob(pe1, pe2, pe3, pe4)
ggsave(file = "../Figures/EE_Error.png", plot = pe,
       width = 10, height = 7, units = "in", dpi = 300)




#### Supplementary Figure: Bias by Threshold ####

estORsAll <- estORsAll %>% mutate(thresholdc = paste0("<", threshold, " SNPs"))

ggplot(data = estORsAll) +
  facet_wrap(~thresholdc) +
  geom_boxplot(aes(x = level, y = logorMean, fill = outcome, color = outcome),
               alpha = 0.5) +
  geom_point(data = trueOR, aes(x = level, y = logorTruth, shape = "True log odds ratio")) +
  scale_shape_discrete(labels = c(expression("True log odds ratio (OR"^"T"*")"))) +
  scale_fill_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                 expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  scale_color_discrete(labels = c(expression("Close genetic relatedness (OR"^"G"*")"),
                                  expression("Naive Bayes modified\nclose genetic relatedness (OR"^"M"*")"))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("Covariate and Level") +
  ylab("Mean Estimated Log Odds Ratio") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggsave("../Figures/EE_Thresholds.png", width = 8, height = 7, dpi = 300)

