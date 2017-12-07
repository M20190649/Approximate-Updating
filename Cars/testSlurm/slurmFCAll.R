rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(RcppEigen, lib.loc = 'packages')
library(rstan, lib.loc = 'packages')
source('slurmRFuns.R')
sourceCpp('slurmCppFuns.cpp')

id <- readRDS('carsID.RDS')
datafc <- readRDS('ForecastData.RDS')
prior <- readRDS('prior.RDS')
H <- 30
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
results <- data.frame()
methods <- c('None', 'Single Hierarchy', 'Finite Mixture')

starting <- list(matrix(c(-5, -5, rep(0, 4), c(chol(diag(0.5, 6)))), ncol = 1),
                 prior[[3]])


# MCMC Hyper parameters
hyper <- list()
for(k in 1:2){
  hyper[[k]] <- list()
  hyper[[k]]$mean <- prior[[k]][1:6]
  uinv <- solve(matrix(prior[[k]][7:42], 6))
  hyper[[k]]$varInv <- t(uinv) %*% uinv
}
# Mixture Prior
mean <- prior[[3]][1:36]
varinv <- NULL
for(k in 1:6){
  uinv <- matrix(prior[[3]][k*36 + 1:36], 6)
  varinv <- rbind(varinv, t(uinv) %*% uinv)
}
weights <- prior[[3]][6*7*6 + 1:6]
weights <- exp(weights) / sum(exp(weights))
hyper[[3]] <- list(mean = mean, varInv = varinv, weights = weights)

# Extract Data
for(iter in 1:3){
i <- sample(873, 1)
  
data <- datafc[[i]]

# Set forcast supports
aLower <- min(data[,1])
if(aLower < 0){
  aLower <- 1.5 * aLower
} else {
  aLower <- 0.5 * aLower
}
dLower <- min(data[,2])
if(dLower < 0){
  dLower <- 1.5 * dLower
} else {
  dLower <- 0.5 * dLower
}

asup <- seq(aLower, 1.5*max(data[,1]), length.out=1000)
dsup <- seq(dLower, 1.5*max(data[,2]), length.out=1000)

# Incrementally add data to VB fits
for(s in 1:15){
  if(sSeq[s] > nrow(data)){
    break
  }
  # Get data for stream
  if(s == 1){
    dataSub <- data[1:sSeq[s],]
  } else {
    dataSub <- data[(sSeq[s-1]+1):sSeq[s],]
  }
  # Update posterior approximations - Or re-estimate new ones from scratch for streaming data
  if(s == 1){
    meanPrior <- matrix(prior[[3]][1:36], 6)
    linvPrior <- array(0, dim = c(6, 6, 6))
    detsPrior <- NULL
    for(k in 1:6){
      uinv <- matrix(prior[[3]][k*36 + 1:36], 6)
      linvPrior[,,k] <- t(uinv)
      detsPrior <- c(detsPrior, prod(diag(uinv)))
    }
    weightsPrior <- prior[[3]][253:258]
    weightsPrior <- exp(weightsPrior) / sum(exp(weightsPrior))
    fitOnline <- carsVBMixScore(dataSub, starting[[2]],
                                   priorMix = list(mean = meanPrior, linv = linvPrior, dets = detsPrior, weights = weightsPrior),
                                   S = 100, maxIter = 10000, threshold = 0.02)$lambda

  } else {
    fitOnline <- carsVBMixScore(dataSub, fitOnline,
                                priorMix = list(mean = meanOn, linv = linvOn, dets = detsOn, weights = weightsOn),
                                S = 100, maxIter = 10000, threshold = 0.05)$lambda
  }

  # Get offline VB posteriors
  if(s == 1){
    fitOffline <- fitOnline
  } else {
    fitOffline <- carsVBMixScore(data[1:sSeq[s],], fitOffline,
                                 priorMix = list(mean = meanPrior, linv = linvPrior, dets = detsPrior, weights = weightsPrior),
                                 S = 100, maxIter = 10000, threshold = 0.02)$lambda
  }
  
  # Get MCMC posteriors
  MCMC <- singleMCMCallMH(data[1:sSeq[s],], 5000, c(-5, -5, 0, 0, 0, 0), hyper[[3]],
                                 stepsize = 0.05, mix = TRUE)$draws
  
  # Extract Lower Triangular Matrices from VB
  LOn <-  array(0, dim = c(6, 6, 6))
  linvOn <- array(0, dim = c(6, 6, 6))
  meanOn <- matrix(fitOnline[1:36], 6)
  detsOn <- NULL
  for(k in 1:6){
    uinv <- matrix(fitOnline[36*k + 1:36], 6)
    l <- invertTri(uinv, lower = FALSE)
    linvOn[,,k] <- t(uinv)
    LOn[,,k] <- l
    detsOn <- c(detsOn, prod(diag(uinv)))
  }
  weightsOn <- exp(fitOnline[253:258]) / sum(exp(fitOnline[253:258]))
  densitiesOnline <- evalVBMix(data[(sSeq[s]-1):sSeq[s], ], meanOn, LOn, weightsOn, 200, H, asup, dsup)
  
  # Repeat for offline
  if(s == 1){
    densitiesOffline <- densitiesOnline
  } else {
    LOff <- array(0, dim = c(6, 6, 6))
    meanOff <- matrix(fitOffline[1:36], 6)
    for(k in 1:6){
      uinv <- matrix(fitOffline[36*k + 1:36], 6)
      l <- invertTri(uinv, lower = FALSE)
      LOff[,,k] <- l
    }
    weightsOff <- exp(fitOffline[253:258]) / sum(exp(fitOffline[253:258]))
    densitiesOffline <- evalVBMix(data[(sSeq[s]-1):sSeq[s], ], meanOff, LOff, weightsOff, 200, H, asup, dsup)
  }
  # Finally MCMC
  densitiesMCMC <- evalMCMCDens(data[(sSeq[s]-1):sSeq[s], ], 1000, H, asup, dsup, MCMC)
  
 
  # Grab logscores for each method, h, and variable.
  #for(k in 1:2){
    for(h in 1:H){
      aindex <- min(which(asup > data[sSeq[s]+h,1]))
      dindex <- min(which(dsup > data[sSeq[s]+h,2]))
      
      alogscoreVBOff <- log(densitiesOffline[(h-1)*1000 + aindex, 1])#, k])
      dlogscoreVBOff <- log(densitiesOffline[(h-1)*1000 + dindex, 2])#, k])
      alogscoreVBOn <- log(densitiesOnline[(h-1)*1000 + aindex, 1])#, k])
      dlogscoreVBOn <- log(densitiesOnline[(h-1)*1000 + dindex, 2])#, k])
      alogscoreMCMC <- log(densitiesMCMC[(h-1)*1000 + aindex, 1])#, k])
      dlogscoreMCMC <- log(densitiesMCMC[(h-1)*1000 + dindex, 2])#, k])
     
    # Attach results
      results <- rbind(results, 
                       data.frame(logscore = c(alogscoreVBOff, dlogscoreVBOff, alogscoreVBOn, dlogscoreVBOn, alogscoreMCMC, dlogscoreMCMC),
                                  variable = rep(c('a', 'd'), 3),
                                  method = rep(c('VB-Offline', 'VB-Stream', 'MCMC'), rep(2, 3)),
                                  prior = 'Finite Mixture',#methods[k],
                                  S = sSeq[s],
                                  h = h,
                                  id = id$idfc[i]))
    }
  #}
  print(paste(i, s))
}
}

results %>%
  filter(logscore > -8) %>% 
  ggplot() + geom_boxplot(aes(x = factor(ceiling(S/30)), y = logscore, colour = method)) + facet_wrap(~variable, scales = 'free')

write.csv(results, paste0('eval/car', id$idfc[[i]], '.csv'), row.names=FALSE)
