#Sarah V. Leavitt
#Boston University Dissertation
#Paper 3

################################################################################
# This program runs one iteration of a simulation
################################################################################


simRunEst <- function() {

  #Simulate outbreak  
  obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                     w.scale = w.scale, w.shape = w.shape, w.shift = w.shift,
                     ws.scale = ws.scale, ws.shape = ws.shape, ws.shift = ws.shift,
                     sampleSize = sampleSize, time = time, multOutbreaks = multOutbreaks,
                     length = length, rate = rate)
  indData <- obk[[1]]
  pairData <- obk[[2]]
  print(paste0("Simulated outbreak, n = ", nrow(indData)))
  
  #Simulating covariates
  covar <- simCovariates(indData, pairData, observationDate = observationDate)
  covarPair <- covar[[1]]
  covarInd <- covar[[2]]
  print("Simulated covariates")
  
  #Subseting to the pairs with the potential infector observed before the infectee
  covarOrderedPair <- covarPair %>% filter(observationDate.2 > observationDate.1)

  
  #Function to find true ORs
  findORs <- function(outcome, l = 1){
    est <- NULL
    #Looping through all covariates
    for(i in 1:length(covariates)){
      #Extracting covariate name
      var <- covariates[i]
      #Creating a table with proportions in each level of the covariate from training data
      tab <- prop.table(table(covarOrderedPair[, var], covarOrderedPair[, outcome]) + l, 2)
      odds <- tab[, 2] / tab[, 1]
      logorMean <- log(odds/odds[1])
      num <- table(covarOrderedPair[, var], covarOrderedPair[, outcome]) + l
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
  
     
    
  #### Choosing pTraining cases for training set ####
  
  #Make sure that the training dataset has at least 3 true links
  #(relevant for small sample sizes)
  nLinked <- 0
  nTries <- 0
  while(nLinked < 3){
    #Finding all pairs that can be included in the training dataset (no missing variables)
    trainingID <- (covarInd
                   %>% filter(complete == TRUE, !is.na(sampleDate))
                   %>% sample_frac(pTraining)
                   %>% pull(individualID)
    )
    
    covarOrderedPair <- (covarOrderedPair
                         %>% mutate(trainPair = ifelse(individualID.1 %in% trainingID &
                                                         individualID.2 %in% trainingID,
                                                       TRUE, FALSE),
                                    
                                    transmissionGS = ifelse(trainPair == TRUE, transmission, NA))
    )
    nLinked <- sum(covarOrderedPair$trainPair == TRUE & covarOrderedPair$transmission == TRUE)
    nTries <- nTries + 1
  }
  if(nTries > 1){print(nTries)}
  
  #Training: True Transmission
  res1 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                          pairIDVar = "edgeID", goldStdVar = "transmissionGS",
                          covariates = covariates, label = paste0("TruthNB"),
                          n = 10, m = 1, nReps = 20)
  probs1 <- res1$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
  print("Completed transmission gold standard analysis")
  
  #Saving coefficients
  est1 <- res1$estimates
  est1$outcome <- "transmissionNB"
  
  
  #### Looping over lower SNP thresholds ####
  
  thresholds <- c(2, 3, 4, 5)
  estTemp <- NULL
  rTemp <- NULL
  
  for(lowerT in thresholds){
    
    #Creating a label that identifies the set of inputs
    label <- paste0("T", lowerT)
    
    #Creating gold standard variables that take into account if the pair should be in the training dataset
    covarOrderedPair <- (covarOrderedPair
                         %>% mutate(snpClose = ifelse(is.na(snpDist), NA,
                                               ifelse(snpDist < lowerT, TRUE, FALSE)),
                                    snpCloseGS = ifelse(trainPair == TRUE & snpDist < lowerT, TRUE,
                                                 ifelse(trainPair == TRUE & snpDist > thresholds[2],
                                                        FALSE, NA)))
    )
    
    #Training: SNP Distance
    res2 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                            pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                            covariates = covariates, label = paste0("SNPs_", label),
                            n = 10, m = 1, nReps = 20)
    probs2 <- res2$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    
    rTemp <- bind_rows(rTemp, probs1, probs2)
    
    #Saving the covariate estimates
    est2 <- res2$estimates
    est2$threshold <- lowerT
    est2$outcome <- "snpCloseNB"
    
    #Finding the ORs using snpCloseGS
    estGS <- findORs("snpCloseGS")
    estGS$threshold <- lowerT
    
    #Finding the ORs using snpClose
    estS <- findORs("snpClose")
    estS$threshold <- lowerT
    
    #Finding sensitivity and specificity
    sens <- sum(covarOrderedPair$transmission == TRUE & 
                  covarOrderedPair$snpClose == TRUE, na.rm = TRUE) /
      sum(covarOrderedPair$transmission == TRUE, na.rm = TRUE)
    ppv <- sum(covarOrderedPair$transmission == TRUE & 
                 covarOrderedPair$snpClose == TRUE, na.rm = TRUE) /
      sum(covarOrderedPair$snpClose == TRUE, na.rm = TRUE)
    
    estComb <- bind_rows(est2, estGS, estS)
    estComb$sens <- sens
    estComb$ppv <- ppv
    
    estTemp <- bind_rows(estTemp, estComb)
    
    print(paste0("Completed SNP gold standard analysis with threshold ", lowerT))
  }
  
  
  ## Evaluating the performance ##
  pTemp <- (rTemp
            %>% group_by(label)
            %>% do(simEvaluate(., truthVar = "transmission"))
            %>% ungroup()
            %>% mutate(threshold = as.numeric(str_extract(str_extract(label, "T[:digit:]+\\.*[:digit:]*"),
                                                          "[:digit:]+\\.*[:digit:]*")),
                       goldStd = ifelse(grepl("Truth", label), "transmission", 
                                 ifelse(grepl("SNPs", label), "snpClose", NA)))
            
  )
  
  #Finding the true ORs for transmission
  estT <- findORs("transmission")

  allEst <- (estTemp
               %>% bind_rows(est1, estT)
               %>% mutate(goldStd = gsub("[A-Z]{2}$", "", outcome))
  )
    
  return(list(allEst, pTemp))
} 



