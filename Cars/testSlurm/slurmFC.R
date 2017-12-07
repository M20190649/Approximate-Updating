rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(RcppEigen, lib.loc = 'packages')
library(rstan, lib.loc = 'packages')
source('slurmRFuns.R')
sourceCpp('slurmCppFuns.cpp')

id <- readRDS('carsID.RDS')
datafc <- readRDS('ForecastData.RDS')
#prior <- readRDS('prior.RDS')
priorMix <- readRDS('priorK6.RDS')
increment <- FALSE
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
results <- data.frame()
methods <- c('None', 'Hierarchy')

#starting <- matrix(c(-5, -5, rep(0, 4), c(chol(diag(0.5, 6)))), ncol = 1)
starting <- matrix(c(rep(c(-5, -5, 0, 0, 0, 0), 6), 
			rep(c(diag(0.5, 6)), 6), rep(1, 6)), ncol = 1)



#fit <- prior

#prior <- list(prior[[1]], prior[[3]])

#hyper <- list()
#for(k in 1:2){
#  hyper[[k]] <- list()
#  hyper[[k]]$mean <- prior[[k]][1:6]
#  uinv <- solve(matrix(prior[[k]][7:42], 6))
#  hyper[[k]]$varInv <- t(uinv) %*% uinv
#}

# Extract Data
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
for(s in seq_along(sSeq)){
  if(sSeq[s] > nrow(data)){
    break
  }
  if(s == 1 | !increment){
    dat <- data[1:sSeq[s],]
  } else {
    dat <- data[(sSeq[s-1]+1):sSeq[s],]
  }
  # Update posterior approximations - Or re-estimate new ones from scratch
 # if(increment){
 #   fit <- fitCarMods(dat, fit, increment, NULL)
 # } else {
 #   fit <- fitCarMods(dat, prior, increment, starting)
 # }
  fit <- carsVBMixScore(dat, starting, priorMix, K = 6, S = 5, maxIter = 10000, threshold = 0.05)$lambda
  # Run MCMC for each method
  #MCMC <- list()
  #for(k in 1:2){
  #  MCMC[[k]] <- singleMCMCallMH(dat, 5000, c(-5, -5, 0, 0, 0, 0), hyper[[k]],
  #                               stepsize = ifelse(k == 1, 0.1, 0.01))$draws
  #}
  MCMC <- singleMCMCallMH(dat, 5000, c(-5, -5, 0, 0, 0, 0), priorMix, mix = TRUE)$draws

  
  
  # Extract Lower Triangular Matrices from VB
  # L <- NULL
  #for(k in 1:2){
  #  L <- rbind(L, t(matrix(fit[[k]][7:42], 6)))
  #}
  #means <- cbind(fit[[1]][1:6],
  #               fit[[2]][1:6])
  #fit[[3]][1:6]) 
  means <- fit[1:36]
  L <- matrix(fit[37:252], ncol = 6, byrow = TRUE)
  w <- exp(fit[253:258]) / sum(exp(fit[253:258]))

  
  #densities <- evalFcDens(data[(sSeq[s]-1):sSeq[s], ], means, L, 1000, S, asup, dsup, MCMC)
  densities <- evalMixDens(data[(sSeq[s]-1):sSeq[s], ], means, L, cumsum(w), 1000, S, asup, dsup, MCMC)
  # Grab logscores for each method, h, and variable.
  #for(k in 1:2){
    for(h in 1:S){
      aindex <- min(which(asup > data[sSeq[s]+h,1]))
      alogscoreVB <- log(densities[(h-1)*1000 + aindex, 1])#, k])
      alogscoreMCMC <- log(densities[(h-1)*1000 + aindex, 3])#, k])
      dindex <- min(which(dsup > data[sSeq[s]+h,2]))
      dlogscoreVB <- log(densities[(h-1)*1000 + dindex, 2])#, k])
      dlogscoreMCMC <- log(densities[(h-1)*1000 + dindex, 4])#, k])
      # Attach results
      results <- rbind(results, 
                       data.frame(logscore = c(alogscoreVB, alogscoreMCMC, dlogscoreVB, dlogscoreMCMC),
                                  variable = c('a', 'a', 'd', 'd'),
                                  method = c('VB', 'MCMC', 'VB', 'MCMC'),
                                  prior = 'Finite Mixture',
                                  S = sSeq[s],
                                  h = h,
                                  id = id$idfc[i]))
    }
  #}
}

write.csv(results, paste0('evalMix/car', id$idfc[[i]], '.csv'), row.names=FALSE)
