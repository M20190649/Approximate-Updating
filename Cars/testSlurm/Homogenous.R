rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
iter <- as.numeric(repenv)

library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
library(RcppEigen, lib.loc = 'packages')
library(rstan, lib.loc = 'packages')
source('slurmRFuns.R')
sourceCpp('slurmCppFuns.cpp')


ids <- readRDS('carsID.RDS')

for(j in 1:10){
 results <- data.frame()
 if(iter == 134 & j > 7){
   break
 }
 i <- (iter-1)*10 + j

if(i <= 873){
  datafc <- readRDS('ForecastData.RDS')
  data <- datafc[[i]]
  id <- ids$idfc[i]
} else {
  datafc <- readRDS('fcDataChanged.RDS')
  data <- datafc[[i-873]]
  id <- datafc[[501]][i-873]
}

draws <- readRDS('noHierN2000.RDS')$draws

H <- 30
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
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


# Incrementally add data
for(s in seq_along(sSeq)){
  if((sSeq[s] + H) > nrow(data)){
    break
  }
  output <-  MCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
           N = 200,
           H = H,
           grid = grid,
           MCMCdraws = draws)

 results <- rbind(results,
                 data.frame(logscore = output$logscore,
                            xCDF = output$xCDF,
                            mapDist = output$mapDist,
                            consDist = NA,
                            h = 1:H,
                            prior = 'Homogenous',
                            method = 'MCMC',
                            S = sSeq[s],
                            id = id))
}


write.csv(results, paste0('homog/car', id, '.csv'), row.names=FALSE)
}
