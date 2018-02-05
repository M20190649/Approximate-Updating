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
H <- 10
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
results <- data.frame()
methods <- c('None', 'Single Hierarchy', 'Finite Mixture')

# MCMC Hyper parameters
hyper <- list()
for(k in 1:2){
  hyper[[k]] <- list()
  hyper[[k]]$mean <- prior[[k]][1:6]
  uinv <- solve(matrix(prior[[k]][7:42], 6))
  hyper[[k]]$varInv <- t(uinv) %*% uinv
}
# Mixture Prior - The prior.RDS object has mean / inverse u however this is converted to mean / log(sd) for the diagonal mixture model in VB
mean <- prior[[3]][1:36]
priorDiag <- mean
varinv <- NULL
for(k in 1:6){
  uinv <- matrix(prior[[3]][k*36 + 1:36], 6)
  vari <- t(uinv) %*% uinv
  logsd <- log(diag(solve(vari)))
  varinv <- rbind(varinv, vari)
  priorDiag <- c(priorDiag, logsd)
}
weights <- prior[[3]][6*7*6 + 1:6]
priorDiag <- c(priorDiag, weights)
weights <- exp(weights) / sum(exp(weights))
hyper[[3]] <- list(mean = mean, varInv = varinv, weights = weights)

prior[[3]] <- matrix(priorDiag, ncol = 1)

starting <- list(matrix(c(-5, -5, rep(0, 4), c(chol(diag(0.5, 6)))), ncol = 1),
                 prior[[3]])



# Extract Data
data <- datafc[[i]]

# Set forcast supports
aLower <- min(data[,1])
if(aLower < 0){
  aLower <- 2 * aLower
} else {
  aLower <- aLower - 0.5
}
dLower <- min(data[,2])
if(dLower < 0){
  dLower <- 2 * dLower
} else {
  dLower <- 0.5 * dLower
}

asup <- seq(aLower, 2*max(data[,1]), length.out = 200)
dsup <- seq(dLower, 2*max(data[,2]), length.out=200)
grid <- cbind(asup, dsup)

# Incrementally add data to VB fits
for(s in seq_along(sSeq)){
  if((sSeq[s] + H) > nrow(data)){
    break
  }
  # Get data for stream
  if(s > 1){
    dataSub <- data[(sSeq[s-1]-1):sSeq[s], 1:2]
  }
  # Update posterior approximations - Or re-estimate new ones from scratch for streaming data
  if(s == 1){
    fitOnline <- fitCarMods(data[1:sSeq[s], 1:2], prior, starting)
    fitOffline <- fitOnline
  } else {
    fitOnline <- fitCarMods(dataSub, fitOnline, list(fitOnline[[1]], fitOnline[[3]]))
    fitOffline <- fitCarMods(data[1:sSeq[s], 1:2], prior, list(fitOffline[[1]], fitOffline[[3]]))
  }
  
  # Get MCMC posteriors
  MCMC <- list()
  for(k in 1:3){
    MCMC[[k]] <- singleMCMCallMH(data[1:sSeq[s],1:2], 5000, c(-5, -5, 0, 0, 0, 0), hyper[[k]],
                                 stepsize = 0.05, mix = (k == 3))$draws
  }
  
  # Evaluate predictive densities
  resultsOnline <- lapply(1:3, function(x) VBDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                                  fit = fitOnline[[x]],
                                                  grid = grid, 
                                                  H = H,
                                                  mix = (x == 3)))
  if(s == 1){
    resultsOffline <- resultsOnline
  } else {
    resultsOffline <- lapply(1:3, function(x) VBDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                                     fit = fitOffline[[x]],
                                                     grid = grid, 
                                                     H = H,
                                                     mix = (x == 3)))
  }
  resultsMCMC <- lapply(MCMC, function(x) MCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                                   N = 200,
                                                   H = H,
                                                   grid = grid,
                                                   MCMCdraws = x))
  # Predict constant angle / velocity distances
  dConst <- data[sSeq[s], 2]
  vConst <- data[sSeq[s] + 1, 3] # Velocity is lagged in the data
  xConst <- vConst * cos(3.141593/2 + dConst)
  yConst <- vConst * sin(3.141593/2 + dConst)
  dist <- numeric(0)
  for(h in 1:H){
    dist <- c(dist, sqrt((h*xConst - sum(data[sSeq[s] + 1:h, 4]))^2 + (h*yConst - sum(data[sSeq[s] + 1:h, 5]))^2))
  }
  
  # Grab logscores etc. 
  for(k in 1:3){
    results <- rbind(results,
                     data.frame(logscore = c(resultsOnline[[k]]$logscore, resultsOffline[[k]]$logscore, resultsMCMC[[k]]$logscore),
                                xCDF = c(resultsOnline[[k]]$xCDF, resultsOffline[[k]]$xCDF, resultsMCMC[[k]]$xCDF),
                                mapDist = c(resultsOnline[[k]]$mapDist, resultsOffline[[k]]$mapDist, resultsMCMC[[k]]$mapDist),
                                consDist = dist,
                                h = rep(1:H, 3),
                                prior = methods[k],
                                method = rep(c('VB-Stream', 'VB-Standard', 'MCMC'), rep(H, 3)),
                                S = sSeq[s],
                                id = id$idfc[i]))
  }
  print(paste(i, s))
}

write.csv(results, paste0('eval/car', id$idfc[[i]], '.csv'), row.names=FALSE)