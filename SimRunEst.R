#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program runs one iteration of a simulation
################################################################################


simRun <- function(multOutbreaks, observationDate) {

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

  
  #Initializing dataframes to hold results, and coefficients
  rTemp <- NULL
  cTemp <- NULL
  
  #List of sets of covariates to test
  covariates <- c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "timeCat")
  
  #Vector of proportions of cases to use in training dataset 
  trainingSizes <- 0.6
  
      
  ## Looping over proportion in the training dataset ##
  for(pTraining in trainingSizes){
    
    #Creating a label that identifies the set of inputs
    label <- paste0("S", pSampling, "T", pTraining, "C", covarI)
    
    
    #### Choosing pTraining cases for training set ####
    
    #Make sure that the training dataset has at least 3 true links (relevant for small sample sizes)
    nLinked <- 0
    nTries <- 0
    while(nLinked < 3){
      #Finding all pairs that can be included in the training dataset (no missing variables)
      trainingID <- (finalInd
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
                         %>% mutate(snpClose = ifelse(snpDist < thresholds[1], TRUE,
                                                      ifelse(snpDist > thresholds[2], FALSE, NA)),
                                    snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA),
                                    transmissionGS = ifelse(trainPair == TRUE, transmission, NA))
    )
    
    #Training: True Transmission
    res1 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                            pairIDVar = "edgeID", goldStdVar = "transmissionGS",
                            covariates = covarList[[covarI]], label = paste0("Truth_", label),
                            n = 10, m = 1, nReps = 10)
    probs1 <- res1$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    print("Completed transmission gold standard analysis")
    
    
    #Training: SNP Distance
    res2 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                            pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                            covariates = covarList[[covarI]], label = paste0("SNPs_", label),
                            n = 10, m = 1, nReps = 10)
    probs2 <- res2$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    print("Completed SNP threshold gold standard analysis")
    
    rTemp <- bind_rows(rTemp, probs1, probs2)
    cTemp <- bind_rows(cTemp, res1[[2]], res2[[2]])
  }
  print(paste0("Completed analysis with training proportion ", pTraining))
  
  
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
      level <- paste(Var, names(orMean), sep = ":")
      cTemp <- cbind.data.frame(level, orMean, orCILB, orCIUB, stringsAsFactors = FALSE)
      coeff <- bind_rows(coeff, cTemp)
    }
    return(coeff)
  }
  
  #Finding the true ORs for snpClose and transmission
  coeffT <- findORs("transmission")
  coeffTGS <- findORs("transmissionGS")
  coeffS <- findORs("snpClose")
  coeffSGS <- findORs("snpCloseGS")
  
  
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  rTemp <- (rTemp
            %>% group_by(label, individualID.2)
            %>% arrange(desc(pScaled))
            %>% mutate(pRank = rank(desc(pScaled), ties.method = "min"))
            %>% ungroup()
  )
  
  ## Evaluating the performance ##
  pTemp <- (rTemp
            %>% group_by(label)
            %>% do(simEvaluate(., truthVar = "transmission"))
            %>% ungroup()
  )
    
  return(list(rTemp, cTemp))
} 



