library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(ggplot2)
library(microbenchmark)
library(mvtnorm)
library(tidyr)
sourceCpp("DLM_MCMC.cpp")
sourceCpp("DLM_SGA_FR.cpp")

mu = 2
phi = 0.95
sigmaSqY = 1
sigmaSqX = 1
alphaY = 1
betaY = 1
alphaX = 1
betaX = 1
muBar = 0
muVar = 10

VBmodelfit = function(Tvec, reps, minReps, meanfield=FALSE){
  output=data.frame()
  
  for(i in 1:length(Tvec)){
    for(j in 1:reps){
      T=Tvec[i]
      MCMCreps = min(minReps, T*10)
      x0 = rnorm(1, 0, sqrt(sigmaSqX/(1-phi^2)))
      x = rep(0, T)
      y = rep(0, T)
      for(t in 1:T){
        if(t == 1){
          x[t] = phi*x0 + rnorm(1, 0, sqrt(sigmaSqX)) 
        } else {
          x[t] = phi*x[t-1] + rnorm(1, 0, sqrt(sigmaSqX))
        }
        y[t] = mu + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
      }
      initM = c(0, 0, cor(y[2:T],y[2:T-1]), mean(y), rep(0, T+1))
      MCMC = DLM_MCMC(y, MCMCreps)
      output = rbind(output, c(colMeans(cbind(log(MCMC$theta[(MCMCreps/2+1):MCMCreps, 1]),
                                              log(MCMC$theta[(MCMCreps/2+1):MCMCreps, 2]),
                                              MCMC$theta[(MCMCreps/2+1):MCMCreps, 3],
                                              MCMC$theta[(MCMCreps/2+1):MCMCreps, 4])), T, 1))
      if(!meanfield){
        FRVB = DLM_SGA(y=y[1:T], S=T, M=25, maxIter=5000, initialM=initM, initialL=diag(0.1, T+5))
        output = rbind(output, c(FRVB$Mu[1:4], T, 3))
      }
      
      MFVB = DLM_SGA(y=y[1:T], S=T, M=5, maxIter=5000, initialM=initM, initialL = diag(0.1, T+5), meanfield=TRUE)
      output = rbind(output, c(MFVB$Mu[1:4], T, 2))
      
      if(j %% 1 == 0){
        print(paste0('T = ', T, '; rep = ', j))
      }
    }
  }
  colnames(output) = c('MeanSigSgY', 'MeanSigSqX', 'MeanPhi', 'MeanGamma', 'T', 'Method')
  output = dplyr::mutate(output, Method=ifelse(Method==1, 'MCMC', ifelse(Method==2, 'Diagonal', 'Non-Diagonal')))
  return(output)
}

#ADVI = VBmodelfit(c(50, 100), 100, 5000)
#saveRDS(ADVI, 'ADVIp1.rds')

#testFit %>% gather(Variable, Mean, -Method, -T) -> testFit

#ggplot(testFit) + geom_boxplot(aes(x=Method, y=Mean)) + facet_grid(Variable~T, scales='free')

#ADVI2 = VBmodelfit(c(250, 500), 100, 5000)
#saveRDS(ADVI2, 'ADVIp2.rds')

#ADVI3 = VBmodelfit(c(1000, 2000), 100, 10000, meanfield=TRUE)
#saveRDS(ADVI3, 'ADVIp3.rds')

ADVI1 = readRDS('ADVIp1.rds')
ADVI2 = readRDS('ADVIp2.rds')
ADVI3 = readRDS('ADVIp3.rds')

ADVI = rbind(ADVI1, ADVI2, ADVI3) %>% gather(Variable, Mean, -Method, -T)
ggplot(ADVI) + geom_boxplot(aes(x=Method, y=Mean)) + facet_grid(Variable~T, scales='free')

