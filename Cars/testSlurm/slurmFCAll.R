rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp)#, lib.loc = 'packages')
library(RcppArmadillo)#, lib.loc = 'packages')
library(RcppEigen)#, lib.loc = 'packages')
library(rstan)#, lib.loc = 'packages')
source('slurmRFuns.R')
sourceCpp('slurmCppFuns.cpp')

ids <- readRDS('carsID.RDS')

if(i <= 873){
  datafc <- readRDS('ForecastData.RDS')
  data <- datafc[[i]]
  id <- ids$idfc[i]
} else {
  datafc <- readRDS('fcDataChanged.RDS')
  data <- datafc[[i-873]]
  id <- datafc[[501]][i-873]
}
homogDraws <- readRDS('noHierN2000.RDS')$draws
prior <- readRDS('prior.RDS')
H <- 30
S <- 10
maxT <- 450

sSeq <- seq(100, maxT, S)
results <- data.frame()
methods <- c('Non-Informative', 'Single Hierarchy', 'Finite Mixture')

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

asup <- seq(aLower, 2*max(data[,1]), length.out = 300)
dsup <- seq(dLower, 2*max(data[,2]), length.out= 300)
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
                                                   N = 300,
                                                   H = H,
                                                   grid = grid,
                                                   MCMCdraws = x))
  
  # Grab logscores etc. for heterogenous models
  for(k in 1:3){
    results <- rbind(results,
                     data.frame(logscore = c(resultsOnline[[k]]$logscore, resultsOffline[[k]]$logscore, resultsMCMC[[k]]$logscore),
                                xCDF = c(resultsOnline[[k]]$xCDF, resultsOffline[[k]]$xCDF, resultsMCMC[[k]]$xCDF),
                                dist = c(resultsOnline[[k]]$mapDist, resultsOffline[[k]]$mapDist, resultsMCMC[[k]]$mapDist),
                                h = rep(1:H, 3),
                                model = rep(paste(methods[k], c('VB-Updating', 'VB-Standard', 'MCMC')), rep(H, 3)),
                                S = sSeq[s],
                                id = id))
  }
  homogResults <-  MCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                      N = 200,
                      H = H,
                      grid = grid,
                      MCMCdraws = homogDraws)
  
  results <- rbind(results,
                   data.frame(logscore = homogResults$logscore,
                              xCDF = homogResults$xCDF,
                              dist = homogResults$mapDist,
                              h = 1:H,
                              model = 'Homogenous MCMC',
                              S = sSeq[s],
                              id = id))
  
  # Naive Model Forecasts
  dConst <- data[sSeq[s], 2]
  aConst <- data[sSeq[s] +1, 3] - data[sSeq[s] , 3]
  aAvg <- 0.1 * (data[sSeq[s] + 1, 3] - data[sSeq[s] -9, 3])
  dAvg <- mean(data[-9:0 + sSeq[s], 2])
  # Model 1: Past A, Const D
  # Model 2: Past A, Zero D
  # Model 3: Past A, Mean D
  # Model 4-6: Mean A, same order of D
  # Model 7-9: Const V (Zero A), same order of D
  x0 <- rep(data[sSeq[s], 4], 9)
  y0 <- rep(data[sSeq[s], 5], 9)
  v0 <- rep(data[sSeq[s] + 1, 3], 9) 
  for(h in 1:H){
    v0 <- v0 + c(rep(aConst, 3), rep(aAvg, 3), rep(0, 3))
    x0 <- x0 + v0 * rep(c(cos(dConst + pi/2), cos(pi/2), cos(pi/2 + dAvg)), 3)
    y0 <- y0 + v0 * rep(c(sin(dConst + pi/2), sin(pi/2), sin(pi/2 + dAvg)), 3)
    dist <- sqrt((x0 - sum(data[sSeq[s] + 1:h, 4]))^2 + (y0 - sum(data[sSeq[s] + 1:h, 5]))^2)
    results <- rbind(results,
                     data.frame(logscore = NA,
                                xCDF = NA,
                                dist = dist,
                                h = h,
                                model = paste('Naive', 1:9),
                                S = sSeq[s],
                                id = id))
  }
  
  print(paste(i, s))
}

write.csv(results, paste0('evalJoint/car', id, '.csv'), row.names=FALSE)
