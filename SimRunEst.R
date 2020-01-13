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
  covarOrderedPair <- (covarPair
                       %>% filter(observationDate.2 > observationDate.1)
                       %>% mutate(snpClose = ifelse(snpDist < thresholds[1], TRUE,
                                             ifelse(snpDist > thresholds[2], FALSE, NA)))
  )

  
  covariates <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "timeCat")
  #Vector of proportions of cases to use in training dataset 
  trainingSizes <- c(0.4, 0.6, 0.8)
  
  #Initializing dataframes to hold results, and coefficients
  rTemp <- NULL
  cTemp <- NULL
  
      
  ## Looping over proportion in the training dataset ##
  for(pTraining in trainingSizes){
    
    #Creating a label that identifies the set of inputs
    label <- paste0("T", pTraining)
    
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
      
      covarOrderedPair <- covarOrderedPair %>% mutate(trainPair = ifelse(individualID.1 %in% trainingID &
                                                                           individualID.2 %in% trainingID,
                                                                         TRUE, FALSE))
      nLinked <- sum(covarOrderedPair$trainPair == TRUE & covarOrderedPair$transmission == TRUE)
      nTries <- nTries + 1
    }
    if(nTries > 1){print(nTries)}
    
    
    
    #### Estimating probabilities ####
    
    #Creating gold standard variables that take into account if the pair should be in the training dataset
    covarOrderedPair <- (covarOrderedPair
                         %>% mutate(snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA),
                                    transmissionGS = ifelse(trainPair == TRUE, transmission, NA))
    )
    
    #Training: True Transmission
    res1 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                            pairIDVar = "edgeID", goldStdVar = "transmissionGS",
                            covariates = covariates, label = paste0("Truth_", label),
                            n = 10, m = 1, nReps = 10)
    probs1 <- res1$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    print("Completed transmission gold standard analysis")
    
    
    #Training: SNP Distance
    res2 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                            pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                            covariates = covariates, label = paste0("SNPs_", label),
                            n = 10, m = 1, nReps = 10)
    probs2 <- res2$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    print("Completed SNP threshold gold standard analysis")
    
    rTemp <- bind_rows(rTemp, probs1, probs2)
    
    coeff1 <- res1$estimates
    coeff1$pTraining <- pTraining
    coeff1$outcome <- "transmissionNB"
    coeff2 <- res2$estimates
    coeff2$pTraining <- pTraining
    coeff2$outcome <- "snpCloseNB"
    
    cTemp <- bind_rows(cTemp, coeff1, coeff2)
  }
  print(paste0("Completed analysis with training proportion ", pTraining))
  
  
  ## Evaluating the performance ##
  pTemp <- (rTemp
            %>% group_by(label)
            %>% do(simEvaluate(., truthVar = "transmission"))
            %>% ungroup()
            %>% mutate(trainingP = as.numeric(str_extract(str_extract(label, "T[:digit:]+\\.*[:digit:]*"),
                                                          "[:digit:]+\\.*[:digit:]*")))
  )
  
  
  #### Finding the true ORs ####
  
  #Function to find true ORs
  findORs <- function(outcome, l = 1){
    coeff <- NULL
    #Looping through all covariates
    for(i in 1:length(covariates)){
      #Extracting covariate name
      var <- covariates[i]
      #Creating a table with proportions in each level of the covariate from training data
      tab <- prop.table(table(covarOrderedPair[, var], covarOrderedPair[, outcome]) + l, 2)
      odds <- tab[, 2] / tab[, 1]
      orMean <- odds/odds[1]
      num <- table(covarOrderedPair[, var], covarOrderedPair[, outcome]) + l
      se <- NA
      for(i in 1:(nrow(num) - 1)){
        numTab <- num[c(1, i+1), ]
        se <- c(se, sqrt(sum(1/numTab)))
      }
      orCILB = exp(log(orMean) - 1.96 * se)
      orCIUB = exp(log(orMean) + 1.96 * se)
      level <- paste(var, names(orMean), sep = ":")
      cTemp <- cbind.data.frame(level, orMean, orCILB, orCIUB, outcome = outcome,
                                stringsAsFactors = FALSE)
      coeff <- bind_rows(coeff, cTemp)
    }
    return(coeff)
  }
  
  #Finding the true ORs for snpClose and transmission
  coeffT <- findORs("transmission")
  coeffS <- findORs("snpClose")

  allCoeff <- (cTemp
               %>% bind_rows(coeffT, coeffS)
               %>% mutate(goldStd = gsub("[A-Z]{2}$", "", outcome))
  )
    
  return(list(allCoeff, pTemp))
} 



