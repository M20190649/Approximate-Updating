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
prior <- readRDS('priorStream.RDS')
increment <- TRUE
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
results <- data.frame()
methods <- c('None', 'Hierarchy', 'Finite Mixture')

starting <- list(matrix(c(-5, -5, rep(0, 4), c(chol(diag(0.5, 6)))), ncol = 1),
                 matrix(c(rep(c(-5, -5, 0, 0, 0, 0), 6), rep(c(diag(0.5, 6)), 6), rep(1, 6)), ncol = 1))

# Offline VB Prior components
meanOff <- prior[[3]][1:36]
linvOff <- NULL
detsOff <- NULL
for(k in 1:6){
  uinvOff <- matrix(prior[[3]][k*36 + 1:36], 6)
  linvOff <- rbind(linvOff, t(uinvOff))
  detsOff <- c(detsOff, det(uinvOff))
}
weightsOff <- prior[[3]][6*7*6 + 1:6]
weightsOff <- exp(weightsOff) / sum(exp(weightsOff))


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
  if(s == 1 | !increment){
    fit <- fitCarMods(dat, prior, starting)
  } else {
    fit <- fitCarMods(dat, fit, starting)
  }

  # Get offline VB Mixture model estimate that broke earlier
  if(s == 1){
    fitOffline <- fit[[3]]
  } else {
    fitOffline <- carsVBMixScore(data[1:sSeq[s],], starting[[2]],  priorMix = list(mean = meanOff, linv = linvOff, dets = detsOff, weights = weightsOff),
                                  S = 40, maxIter = 10000, threshold = 0.05)$lambda
  }
  # Extract Lower Triangular Matrices from VB
  L <- list()
  means <- list()
  for(k in 1:2){
    means[[k]] <- fit[[k]][1:6]
    L[[k]] <- matrix(fit[[k]][7:42], 6)
  }
  means[[3]] <- fit[[3]][1:36]
  # Convert Inverse from mixture lambda to L
  L[[3]] <- matrix(0, ncol = 6, nrow = 36)
  for(k in 1:6){
    L[[3]][(k-1)*6 + 1:6,] <- t(solve(matrix(fit[[3]][36*k + 1:36], 6)))
  }
  w <- exp(fit[[3]][253:258]) / sum(exp(fit[[3]][253:258]))

  densities <- evalVBDens(data[(sSeq[s]-1):sSeq[s], ], means, L, cumsum(w), 1000, S, asup, dsup)
 
  # Grab logscores for each method, h, and variable.
  for(k in 1:3){
    for(h in 1:S){
      aindex <- min(which(asup > data[sSeq[s]+h,1]))
      alogscoreVB <- log(densities[(h-1)*1000 + aindex, 1, k])
      dindex <- min(which(dsup > data[sSeq[s]+h,2]))
      dlogscoreVB <- log(densities[(h-1)*1000 + dindex, 2, k])
     
    # Attach results
      results <- rbind(results, 
                       data.frame(logscore = c(alogscoreVB, dlogscoreVB),
                                  variable = c('a', 'd'),
                                  method = c('VB-Stream'),
                                  prior = methods[k],
                                  S = sSeq[s],
                                  h = h,
                                  id = id$idfc[i]))
    }
  }

  # Get densities and logscores for offline VB model
  means <- matrix(0, 6, 6)
  L <- array(0, dim = c(6, 6, 6))

  # Convert Inverse from mixture lambda to L
  for(k in 1:6){
    means[,k] <- fitOffline[(k-1)*6 + 1:6]
    L[,,k] <- t(solve(matrix(fitOffline[36*k + 1:36], 6)))
  }
  w <- exp(fitOffline[253:258]) / sum(exp(fitOffline[253:258]))
  densities <- evalVBMix(data[(sSeq[s]-1):sSeq[s], ], means, L, cumsum(w), 1000, S, asup, dsup)

  for(h in 1:S){
      aindex <- min(which(asup > data[sSeq[s]+h,1]))
      alogscoreVB <- log(densities[(h-1)*1000 + aindex, 1])
      dindex <- min(which(dsup > data[sSeq[s]+h,2]))
      dlogscoreVB <- log(densities[(h-1)*1000 + dindex, 2])

    # Attach results
      results <- rbind(results,
                       data.frame(logscore = c(alogscoreVB, dlogscoreVB),
                                  variable = c('a', 'd'),
                                  method = 'VB',
                                  prior = 'Finite Mixture',
                                  S = sSeq[s],
                                  h = h,
                                  id = id$idfc[i]))
    }



}

write.csv(results, paste0('evalStream/car', id$idfc[[i]], '.csv'), row.names=FALSE)
