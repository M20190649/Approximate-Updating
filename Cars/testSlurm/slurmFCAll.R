rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- 5#as.numeric(repenv)
set.seed(1000 + i)

library(Rcpp)#, lib.loc = 'packages')
library(RcppArmadillo)#, lib.loc = 'packages')
library(RcppEigen)#, lib.loc = 'packages')
library(rstan)#, lib.loc = 'packages')
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
data <- datafc[[i]]

# Set forcast supports
xLower <- min(data[,4])
if(xLower < 0){
  xLower <- 2 * xLower
} else {
  xLower <- xLower - 0.5
}
yLower <- min(data[,5])
if(yLower < 0){
  yLower <- 1.5 * yLower
} else {
  yLower <- 0.5 * yLower
}

xsup <- seq(min(xLower, -2), max(2, 2*max(data[,4])), length.out=300)
ysup <- seq(yLower, 1.5*max(data[,5]), length.out=300)
grid <- as.matrix(expand.grid(xsup, ysup))

# Incrementally add data to VB fits
for(s in seq_along(sSeq)){
  if(sSeq[s] > nrow(data)){
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
  start <- Sys.time()
  densitiesOnline <- lapply(1:3, function(x) VBDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                                    fit = fitOnline[[x]],
                                                    grid = grid, 
                                                    H = H,
                                                    mix = (x == 3)))
  if(s == 1){
    densitiesOffline <- densitiesOnline
  } else {
    densitiesOffline <- lapply(1:3, function(x) VBDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                                       fit = fitOffline[[x]],
                                                       grid = grid, 
                                                       H = H,
                                                       mix = (x == 3)))
  }
  densitiesMCMC <- lapply(MCMC, function(x) evalMCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
                                                     N = sqrt(nrow(grid)),
                                                     H = H,
                                                     grid = grid,
                                                     MCMCdraws = x))
  # Grab logscores for each method, h, and variable.
  for(k in 1:3){
    for(h in 1:H){
      xindex <- min(which(xsup > data[sSeq[s]+h,4]))
      yindex <- min(which(ysup > data[sSeq[s]+h,5]))
      
      scoreVBOff <- densitiesOffline[[k]][yindex, xindex, h]
      scoreVBOn <- densitiesOnline[[k]][yindex, xindex, h]
      scoreMCMC <- densitiesMCMC[[k]][yindex, xindex, h]
      
      offCDF <- sum(densitiesOffline[[k]][,1:xindex,h]) * (xsup[2] - xsup[1]) * (ysup[2] - ysup[1])
      onCDF <- sum(densitiesOnline[[k]][,1:xindex,h]) * (xsup[2] - xsup[1]) * (ysup[2] - ysup[1])
      MCMCCDF <- sum(densitiesMCMC[[k]][,1:xindex,h]) * (xsup[2] - xsup[1]) * (ysup[2] - ysup[1])
   
      
    # Attach results
      results <- rbind(results, 
                       data.frame(logscore = c(log(scoreVBOff), log(scoreVBOn), log(scoreMCMC)),
                                  xCDF = c(offCDF, onCDF, MCMCCDF),
                                  method = c('VB-Offline', 'VB-Stream', 'MCMC'),
                                  prior = methods[k],
                                  S = sSeq[s],
                                  h = h,
                                  id = id$idfc[i]))
    }
  }
  print(paste(i, s))
}

write.csv(results, paste0('eval/car', id$idfc[[i]], '.csv'), row.names=FALSE)



