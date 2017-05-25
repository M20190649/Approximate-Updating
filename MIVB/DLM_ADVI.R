library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
sourceCpp("DLM_MCMC.cpp")
#sourceCpp("DLMSplit.cpp")

# Parameters
gamma = 0.5
phi = 0.8
sigmaSqY = 1
sigmaSqX = 1

T = 250
MCMCreps = 10000
reps = 4
i=1

VB = data.frame()
MCMCdf = data.frame()
for(i in 1:reps){
  x0 = rnorm(1, gamma / (1-phi), sqrt(sigmaSqX/(1-phi^2)))
  x = rep(0, T)
  y = rep(0, T)
  for(t in 1:T){
    if(t == 1){
      x[t] = gamma + phi*(x0 - gamma) + rnorm(1, 0, sqrt(sigmaSqX)) 
    } else {
      x[t] = gamma + phi*(x[t-1] - gamma) + rnorm(1, 0, sqrt(sigmaSqX))
    }
    y[t] = x[t] + rnorm(1, 0, sqrt(sigmaSqY))
  }
  MCMCfit = DLM_MCMC(y, MCMCreps)
  colnames(MCMCfit$theta) = c('sigma^2[y]', 'sigma^2[x]', 'phi', 'gamma')
  MCMCl = gather(as.data.frame(MCMCfit$theta), variable, draw)
  MCMCl$dataset = i
  MCMCdf = rbind(MCMCdf, MCMCl)
  
  initMuTheta = colMeans(cbind(log(MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 1:2]), MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 3:4]))
  initLTheta = apply(cbind(log(MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 1:2]), MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 3:4]), 2, sd)
  initMuX = colMeans(MCMCfit$x[(MCMCreps/2 + 1):MCMCreps, 1:(T+1)])
  initLX = apply(MCMCfit$x[(MCMCreps/2 + 1):MCMCreps, 1:(T+1)], 2, sd)
  
  initMu = c(initMuTheta, initMuX)
  initL = diag(c(initLTheta, initLX))

  
}

     
      # Initial value of Mean vector is mostly 0, but the first autocorrelation should be a decent starting value
      # for phi, and the mean should be a decent starting value for gamma
      initM = c(0, 0, cor(y[2:T],y[2:T-1]), mean(y), rep(0, T+1))
      
      # Run MCMC, take mean of second half of draws
      MCMC = DLM_MCMC(y, MCMCreps)
      #output = rbind(output, c(colMeans(cbind(log(MCMC$theta[(MCMCreps/2+1):MCMCreps, 1]),
       #                                       log(MCMC$theta[(MCMCreps/2+1):MCMCreps, 2]),
        #                                      MCMC$theta[(MCMCreps/2+1):MCMCreps, 3],
         #                                     MCMC$theta[(MCMCreps/2+1):MCMCreps, 4])), T, 1))
      
      # Non-Diagonal VB, estimating values for all latent states, with 25 simulations per iteration
      # Stopping if it fails to converge after 5000 replications
      # Initial value for the lower triangular matrix of Sigma is 0.1 * I
      #if(!meanfield){
        NDVB = DLM_SGA(y=y[1:T], S=T, M=50, maxIter=5000, initialM=initM, initialL=diag(0.1, T+5))
     #   output = rbind(output, c(NDVB$Mu[1:4], T, 3))
      #}
      
      # Diagonal VB, same inputs except less simulations per iteration
      MFVB = DLM_SGA(y=y[1:T], S=T, M=5, maxIter=5000, initialM=initM, initialL = diag(0.1, T+5), meanfield=TRUE)
     # output = rbind(output, c(MFVB$Mu[1:4], T, 2))
      
      MCMCtheta = MCMC$theta[5001:10000,]
      colnames(MCMCtheta) = c('SigSqY', 'SigSqX', 'Phi', 'Gamma')
      MCMCtheta = gather(as.data.frame(MCMCtheta), variable, draw)
      NDSig = sqrt(NDVB$L %*% t(NDVB$L))
      
      sigxsup = seq(0, 3, length.out=1000)
      sigysup = seq(0, 5, length.out=1000)
      phisup = seq(0, 1, length.out=1000)
      gammasup = seq(1, 3.5, length.out=1000)
      SigSqYd = dlnorm(sigysup, MFVB$Mu[1], MFVB$Sd[1])
      SigSqXd = dlnorm(sigxsup, MFVB$Mu[2], MFVB$Sd[2])
      phid = dnorm(phisup, MFVB$Mu[3], MFVB$Sd[3])
      gammad = dnorm(gammasup, MFVB$Mu[4], MFVB$Sd[4])
      
      MF = data.frame(support = c(sigysup, sigxsup, gammasup, phisup), density = c(SigSqYd, SigSqXd, gammad, phid))
      MF$type = 'Diag'
      MF$variable = rep(c('SigSqY', 'SigSqX', 'Gamma', 'Phi'), rep(1000, 4))
      
      SigSqYd = dlnorm(sigysup, NDVB$Mu[1], NDSig[1, 1])
      SigSqXd = dlnorm(sigxsup, NDVB$Mu[2], NDSig[2, 2])
      phid = dnorm(phisup, NDVB$Mu[3], NDSig[3, 3])
      gammad = dnorm(gammasup, NDVB$Mu[4], NDSig[4, 4])
      
      ND = data.frame(support = c(sigysup, sigxsup, gammasup, phisup), density = c(SigSqYd, SigSqXd, gammad, phid))
      ND$type = 'Non-Diag'
      ND$variable = rep(c('SigSqY', 'SigSqX', 'Gamma', 'Phi'), rep(1000, 4))
      
      
      VB = rbind(MF, ND)
      
      ggplot(VB) + geom_line(aes(support, density, colour=type)) + geom_density(data=MCMCtheta, aes(draw)) + facet_wrap(~variable, scales='free')
      
      
      
      
      
      
      # Print progress
  #    if(j %% 5 == 0){
  #      print(paste0('T = ', T, '; rep = ', j))
  #    }
  #  }
  #}
  #colnames(output) = c('sigma[y]^2','sigma[x]^2', 'phi', 'gamma', 'T', 'Method')
  #output = dplyr::mutate(output, Method=ifelse(Method==1, 'MCMC', ifelse(Method==2, 'ADVI: Diag', 'ADVIL Non0Diag')))
  #return(output)
#}

#ADVI = VBmodelfit(c(50, 100), 100, 5000)
#saveRDS(ADVI, 'ADVIp1.rds')

#testFit %>% gather(Variable, Mean, -Method, -T) -> testFit

#ggplot(testFit) + geom_boxplot(aes(x=Method, y=Mean)) + facet_grid(Variable~T, scales='free')

#ADVI2 = VBmodelfit(c(250, 500), 100, 5000)
#saveRDS(ADVI2, 'ADVIp2.rds')

#ADVI3 = VBmodelfit(c(1000, 2000), 100, 10000, meanfield=TRUE)
#saveRDS(ADVI3, 'ADVIp3.rds')

#ADVI1 = readRDS('ADVIp1.rds')
#ADVI2 = readRDS('ADVIp2.rds')
#ADVI3 = readRDS('ADVIp3.rds')

#ADVI = rbind(ADVI1, ADVI2, ADVI3) %>% gather(Variable, Mean, -Method, -T)

#ggplot(ADVI) + geom_boxplot(aes(x=Method, y=Mean)) + facet_grid(Variable~T, scales='free', labeller = label_parsed) +
 # theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5), strip.text.y = element_text(angle = 0))

