library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(msm)
sourceCpp('MCMC.cpp')
sourceCpp('DLM.cpp')
sourceCpp('SVM.cpp')

MIVB = function(T, S, reps, truev, model = 'DLM'){
  VB = data.frame()
  MCMC = data.frame()
  Xt = data.frame()
  for(i in 1:reps){
    x0 = rnorm(1, truev[4], sqrt(truev[2]/(1-truev[3]^2)))
    x = rep(0, T+S)
    y = rep(0, T+S)
    for(t in 1:(T+S)){
      if(t == 1){
        x[t] = truev[4] + truev[3]*(x0 - truev[4]) + rnorm(1, 0, sqrt(truev[2])) 
      } else {
        x[t] = truev[4] + truev[3]*(x[t-1] - truev[4]) + rnorm(1, 0, sqrt(truev[2])) 
      }
      if(model == 'DLM'){
        y[t] = x[t] + rnorm(1, 0, sqrt(truev[1]))
      } else {
        y[t] = x[t] + log(rnorm(1)^2)
      }
    }
    if(model == 'DLM'){
      MCMCfit = DLM_MCMC(y[1:T], 100000)
      colnames(MCMCfit$theta) = c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma')
    } else {
      MCMCfit = SVM_MCMC(y[1:T], 100000)
      colnames(MCMCfit$theta) = c('sigma^2', 'phi', 'gamma')
    }
    
    MCMCl = gather(as.data.frame(MCMCfit$theta[50001:100000,]), Variable, Draws)
    MCMCl$Dataset = i
    MCMCl <- MCMCl %>% group_by(Variable) %>%
      filter(Draws > quantile(Draws, 0.025) & Draws < quantile(Draws, 0.975)) %>% ungroup()
    grid = MCMCl %>% group_by(Variable) %>% summarise(min=min(Draws), max=max(Draws)) %>% select(min, max) %>% as.matrix()
    
    MCMC = rbind(MCMC, MCMCl)
    if(model == 'DLM'){
      supportPos = seq(min(grid[3:4, 1])-1.5, max(grid[3:4, 2])+1.5, length.out = 500)
    } else {
      supportPos = seq(grid[3, 1]-1.5, grid[3, 2]+1.5, length.out = 500)
    }
    supportReal = seq(min(grid[1:2, 1])-1.5, max(grid[1:2, 2])+1.5, length.out = 500)
    
    if(model == 'DLM'){
      states = apply(MCMCfit$x[50001:100000,1:(T+1)], 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))) %>%
        t() %>% as.data.frame()
      colnames(states) = c('Mean', 'L95', 'U95')
      states$Method = 'MCMC'
      
      initMuTheta = cbind(log(MCMCfit$theta[50001:100000, 1:2]), MCMCfit$theta[50001:100000, 3:4]) %>% colMeans
      initLTheta = cbind(log(MCMCfit$theta[50001:100000, 1:2]), MCMCfit$theta[50001:100000, 3:4]) %>% var %>% chol %>% t
      initMuX = states$Mean
      initLX = apply(MCMCfit$x[50001:100000, 1:(T+1)], 2, sd)
      
      initMu = c(initMuTheta, initMuX, rep(0, S))
      initL = matrix(0, T+S+5, T+S+5)
      initL[1:4, 1:4] = initLTheta
      diag(initL)[5:(T+S+5)] = c(initLX, rep(0.1, S))
      
      VBfit = SGA_DLM(y=y, M=50, maxIter=5000, Mu=initMu, L=initL, S=S, alpha=0.1, meanfield=TRUE, xderiv=TRUE, variance=TRUE)
      
      VBsd = sqrt(diag(VBfit$L %*% t(VBfit$L)))
      dSigmaSqY = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
      dSigmaSqX = dlnorm(supportPos, VBfit$Mu[2], VBsd[2])
      dPhi = dnorm(supportReal, VBfit$Mu[3], VBsd[3])
      dGamma = dnorm(supportReal, VBfit$Mu[4], VBsd[4])
      df = data.frame(c(rep(supportPos, 2), rep(supportReal, 2)),
                    c(dSigmaSqY, dSigmaSqX, dPhi, dGamma), 
                    rep(c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma'), rep(500, 4)), i)
      VB = rbind(VB, df)
      
      vbstates = data.frame(mu=VBfit$Mu[6:(S+5) + T], sd=VBsd[6:(S+5) + T]) %>%
        apply(1, function(x) c(x[1], x[1]-1.96*x[2], x[1]+1.96*x[2])) %>% t() %>% as.data.frame()
    } else {
      weights = thetaWeights(y, MCMCfit$theta[50001:100000, ], MCMCfit$x[50001:100000, ], MCMCfit$s[50001:100000, ])
      
      states = apply(weights * MCMCfit$x[50001:100000,1:(T+1)], 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))) %>%
        t() %>% as.data.frame()
      colnames(states) = c('Mean', 'L95', 'U95')
      states$Method = 'MCMC'
     
      initMuTheta = apply(weights * cbind(log(MCMCfit$theta[50001:100000, 1]), MCMCfit$theta[50001:100000, 2:3]), 2, mean)
      initLTheta = apply(weights * cbind(log(MCMCfit$theta[50001:100000, 1]), MCMCfit$theta[50001:100000, 2:3]), 2, sd)
      initMuX = states$Mean
      initLX = apply(MCMCfit$x[50001:100000, 1:(T+1)], 2, sd)
      initMu = c(initMuTheta, initMuX, rep(0, S))
      initL = diag(c(initLTheta, initLX, rep(0.1, S)))
      
      VBfit = SGA_SVM(y=y, M=50, maxIter=5000, Mu=initMu, L=initL, S=S, alpha=0.1, meanfield=TRUE)
      
      VBsd = diag(VBfit$L)
      dSigmaSq = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
      dPhi = dnorm(supportReal, VBfit$Mu[2], VBsd[2])
      dGamma = dnorm(supportReal, VBfit$Mu[3], VBsd[3])
      df = data.frame(c(supportPos, rep(supportReal, 2)),
                      c(dSigmaSq, dPhi, dGamma),
                      rep(c('sigma^2', 'phi', 'gamma'), rep(500, 3)), 
                      i)
      VB = rbind(VB, df)
      vbstates = data.frame(mu=VBfit$Mu[5:(S+4)+T], sd=VBsd[5:(S+4)+T]) %>%
        apply(1, function(x) c(x[1], x[1]-1.96*x[2], x[1]+1.96*x[2])) %>% t() %>% as.data.frame()
    }
    colnames(vbstates) = c('Mean', 'L95', 'U95')
    vbstates$Method = 'ADVI'
    st = rbind(states, vbstates)
    st$t = 0:(T+S)
    st$Dataset = i
    st$Actual = c(x0, x)
    Xt = rbind(Xt, st)
  }
  colnames(VB) = c('Support', 'Density', 'Variable', 'Dataset')
  return(list(MCMC=MCMC, VB=VB, Xt=Xt))
}



