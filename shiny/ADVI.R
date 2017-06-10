library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(msm)
library(RcppEigen)
library(rstan)
sourceCpp('Variance.cpp')
#sourceCpp('AR1.cpp')
#sourceCpp('AR1Mean.cpp')
sourceCpp('MCMC.cpp')
#sourceCpp('DLM.cpp')
sourceCpp('SVM.cpp')
#sourceCpp('AutoDiffSGA.cpp')

ldgamma = function(x, alpha, beta){
  alpha * log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
}

ADVI = function(T, reps, truev, model, meanfield=TRUE, xderiv=FALSE, var=FALSE){
  VB = data.frame()
  True = data.frame()
  supportReal = seq(min(truev)-2, max(truev)+2, length.out=500)
  supportPos = seq(0.01, max(truev)+2, length.out=500)
  
  if(model=='Normal'){
    palpha=1
    pbeta=1
    for(i in 1:reps){
      y = rnorm(T, 0, truev)
      # Fit from arbitary starting values
      lnmu = 0
      lndelta = -0.75
      iga = 0.1
      igb = 0.5
      # Fit the VB approx for both distributions
      VBln = SGA_Var(y=y, M=50, lognormal=TRUE, initPar1=lnmu, initPar2=lndelta, alpha=0.15)
      VBig = SGA_Var(y=y, M=50, lognormal=FALSE, initPar1=iga, initPar2=igb, alpha=0.75)
      # Evaluate density over a grid for fitted densities, starting densities and true density
      dln = dlnorm(supportPos, VBln$Params[1], exp(VBln$Params[2]))
      dig = exp(ldgamma(supportPos, exp(VBig$Params[1]), exp(VBig$Params[2])))
      dtrue = exp(ldgamma(supportPos, 1+palpha+T/2, pbeta + 0.5*sum(y^2)))
      dlnStart = dlnorm(supportPos, lnmu, exp(lndelta))
      digStart = exp(ldgamma(supportPos, exp(iga), exp(igb)))
      # Collect results
      df = data.frame(supportPos, c(dln, dlnStart, dig, digStart),
                    rep(c('Lognormal', 'Inverse-Gamma'), rep(1000, 2)),
                    rep(c('Converged', 'Starting'), rep(500, 2)), i)
      VB = rbind(VB, df)
      df = data.frame(supportPos, dtrue, rep(c('Lognormal', 'Inverse-Gamma'), rep(500, 2)), i)
      True = rbind(True, df)
    }
    colnames(VB) = c('Support', 'Density', 'Variable', 'Version', 'Dataset')
    colnames(True) = c('Support', 'Density', 'Variable', 'Dataset')
    return(list(ADVI=VB, True=True))
    
  } else if(model=='AR1-Zero Mean') {
    for(i in 1:reps){
      T = 200
      truev = c(1, 0.8, 0.5)
      x = rep(0, T+1)
      x[1] = rnorm(1, 0, sqrt(truev[1] / (1-truev[2]^2)))
      for(t in 2:(T+1)){
        x[t] = rnorm(1, truev[2]*x[t-1], sqrt(truev[1]))
      }
      MCMCfit = as.data.frame(MCMC_AR1(x, 100000))
      colnames(MCMCfit) = c('sigma^2', 'phi')
      MCMCl = tidyr::gather(MCMCfit[50001:100000,], Variable, Draws)
      MCMCl$Dataset = i
      True = rbind(True, MCMCl)
      smu = -1#mean(log(MCMCfit[1:1000, 1]))
      ssd = 1#sd(log(MCMCfit[1:1000, 1]))
      pmu = 0#mean(MCMCfit[1001:2000, 1])
      psd = 0.5#sd(MCMCfit[1001:2000, 1])
      
      initMu = c(smu, pmu)
      initL = diag(c(ssd, psd))
      lM = cbind(initMu, initL)
      xM = matrix(x, T+1)
      VBfit = AD_SGA(xM, lM, lrows=2, M=100, maxIter=0, alpha=0.1)
      VBsd = sqrt(diag(t(VBfit$U) %*% VBfit$U))
      dSigmaSq = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
      dPhi = dnorm(supportReal, VBfit$Mu[2], VBsd[2])
      dsStart = dlnorm(supportPos, smu, ssd)
      dpStart = dnorm(supportReal, pmu, psd)
      
      df = data.frame(c(supportPos, supportReal), c(dSigmaSq, dPhi, dsStart, dpStart), 
                      rep(c('sigma^2', 'phi'), rep(500, 2)), rep(c('Converged', 'Starting'), rep(1000, 2)), i)
      VB = rbind(VB, df)
    }
    colnames(VB) = c('Support', 'Density', 'Variable', 'Version', 'Dataset')
    True <- True %>% group_by(Variable, Dataset) %>%
      filter(Draws > quantile(Draws, 0.025) & Draws < quantile(Draws, 0.975)) %>% ungroup()
    return(list(ADVI=VB, True=True))
    
  } else if (model=='AR1 with mean') {
    for(i in 1:reps){
      x = rep(0, T+1)
      x[1] = rnorm(1, truev[3], sqrt(truev[1] / (1-truev[2]^2)))
      for(t in 2:(T+1)){
        x[t] = rnorm(1, truev[3] + truev[2]*(x[t-1]-truev[3]), sqrt(truev[1]))
      }
      MCMCfit = as.data.frame(MCMC_AR1M(x, 100000))
      colnames(MCMCfit) = c('sigma^2', 'phi', 'gamma')
      MCMCl = tidyr::gather(MCMCfit[50001:100000,], Variable, Draws)
      MCMCl$Dataset = i
      True = rbind(True, MCMCl)
      smu = -1#mean(log(MCMCfit[1:1000, 1]))
      ssd = 0.5#sd(log(MCMCfit[1:1000, 1]))
      pmu = 0#mean(MCMCfit[1001:2000, 1])
      psd = 0.3#sd(MCMCfit[1001:2000, 1])
      gmu = 0#mean(MCMCfit[2001:3000, 1])
      gsd = 1#sd(MCMCfit[2001:3000, 1])
      
      initMu = c(smu, pmu, gmu)
      initL = diag(c(ssd, psd, gsd))
      lM = cbind(initMu, initL)
      xM = matrix(x, T+1)
      VBfit = AD_SGA(xM, lM, lrows=3, M=100, maxIter=5000)
      
      VBsd = sqrt(diag(t(VBfit$U) %*% VBfit$U))
      dSigmaSq = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
      dPhi = dnorm(supportReal, VBfit$Mu[2], VBsd[2])
      dGamma = dnorm(supportReal, VBfit$Mu[3], VBsd[3])
      dsStart = dlnorm(supportPos, smu, ssd)
      dpStart = dnorm(supportReal, pmu, psd)
      dgStart = dnorm(supportReal, gmu, gsd)
      
      df = data.frame(c(supportPos, supportReal, supportReal), c(dSigmaSq, dPhi, dGamma, dsStart, dpStart, dgStart), 
                      rep(c('sigma^2', 'phi', 'gamma'), rep(500, 3)), rep(c('Converged', 'Starting'), rep(1500, 2)), i)
      VB = rbind(VB, df)
    }
    colnames(VB) = c('Support', 'Density', 'Variable', 'Version', 'Dataset')
    True <- True %>% group_by(Variable, Dataset) %>%
      filter(Draws > quantile(Draws, 0.025) & Draws < quantile(Draws, 0.975)) %>% ungroup()
    return(list(ADVI=VB, True=True))
    
  } else if(model == 'DLM'){
    VB = data.frame()
    True = data.frame()
    Xt = data.frame()
    for(i in 1:reps){
      x0 = rnorm(1, truev[4] / (1-truev[3]), sqrt(truev[2]/(1-truev[3]^2)))
      x = rep(0, T)
      y = rep(0, T)
      for(t in 1:T){
        if(t == 1){
          x[t] = truev[4] + truev[3]*(x0 - truev[4]) + rnorm(1, 0, sqrt(truev[2])) 
        } else {
          x[t] = truev[4] + truev[3]*(x[t-1] - truev[4]) + rnorm(1, 0, sqrt(truev[2])) 
        }
        y[t] = x[t] + rnorm(1, 0, sqrt(truev[1]))
      }
      MCMCfit = DLM_MCMC(y, 100000)
      colnames(MCMCfit$theta) = c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma')
      MCMCl = gather(as.data.frame(MCMCfit$theta[50001:100000,]), Variable, Draws)
      MCMCl$Dataset = i
      True = rbind(True, MCMCl)
      
      states = apply(MCMCfit$x[50001:100000,1:(T+1)], 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))) %>%
        t() %>% as.data.frame()
      colnames(states) = c('Mean', 'L95', 'U95')
      states$Method = 'MCMC'
  
      initMuTheta = colMeans(cbind(log(MCMCfit$theta[50001:100000, 1:2]), MCMCfit$theta[50001:100000, 3:4]))
      initLTheta = apply(cbind(log(MCMCfit$theta[50001:100000, 1:2]), MCMCfit$theta[50001:100000, 3:4]), 2, sd)
      initMuX = states$Mean
      initLX = apply(MCMCfit$x[50001:100000, 1:(T+1)], 2, sd)
      initMuTheta = c(-1.5, -1.5, 0, -1)
      if(xderiv) initMuX = rep(0, T+1)
      if(var) initLTheta = rep(0.1, 4)
      if(xderiv & var) initLX = rep(0.1, T+1)
      initMu = c(initMuTheta, initMuX)
      initL = diag(c(initLTheta, initLX))
      lM = cbind(initMu, initL)
      yM = matrix(y, T)
      
      VBfit = AD_SGA(yM, lM, T+5, M=100, maxIter=5000, meanfield=TRUE)
      
      VBsd = sqrt(diag(VBfit$L %*% t(VBfit$L)))
      dSigmaSqY = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
      dSigmaSqX = dlnorm(supportPos, VBfit$Mu[2], VBsd[2])
      dPhi = dnorm(supportReal, VBfit$Mu[3], VBsd[3])
      dGamma = dnorm(supportReal, VBfit$Mu[4], VBsd[4])
      dsyStart = dlnorm(supportPos, initMuTheta[1], initLTheta[1])
      dsxStart = dlnorm(supportPos, initMuTheta[2], initLTheta[2])
      dpStart = dnorm(supportReal, initMuTheta[3], initLTheta[3])
      dgStart = dnorm(supportReal, initMuTheta[4], initLTheta[4])
      
      vbstates = data.frame(mu=VBfit$Mu[5:(T+5)], sd=VBsd[5:(T+5)]) %>%
        apply(1, function(x) c(x[1], x[1]-1.96*x[2], x[1]+1.96*x[2])) %>% t() %>% as.data.frame()
      colnames(vbstates) = c('Mean', 'L95', 'U95')
      vbstates$Method = 'ADVI'
      st = rbind(states, vbstates)
      st$t = 0:T
      st$Dataset = i
      st$Actual = c(x0, x)
      Xt = rbind(Xt, st)

      df = data.frame(c(rep(supportPos, 2), rep(supportReal, 2)),
                        c(dSigmaSqY, dSigmaSqX, dPhi, dGamma, dsyStart, dsxStart, dpStart, dgStart), 
                        rep(c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma'), rep(500, 4)), 
                        rep(c('converged', 'starting'), rep(2000, 2)), i)
      VB = rbind(VB, df)
    }
    colnames(VB) = c('Support', 'Density', 'Variable', 'Version', 'Dataset')
    True <- True %>% group_by(Variable, Dataset) %>%
      filter(Draws > quantile(Draws, 0.025) & Draws < quantile(Draws, 0.975)) %>% ungroup()
    return(list(ADVI=VB, True=True, Xt=Xt))
    
  } else if(model == 'SVM'){
    VB = data.frame()
    True = data.frame()
    Xt = data.frame()
    for(i in 1:reps){
      x0 = rnorm(1, truev[3] / (1-truev[2]), sqrt(truev[1]/(1-truev[2]^2)))
      x = rep(0, T)
      y = rep(0, T)
      for(t in 1:T){
        if(t == 1){
          x[t] = truev[3] + truev[2]*(x0 - truev[3]) + rnorm(1, 0, sqrt(truev[1])) 
        } else {
          x[t] = truev[3] + truev[2]*(x[t-1] - truev[3]) + rnorm(1, 0, sqrt(truev[1])) 
        }
        y[t] = x[t] + log(rnorm(1)^2)
      }
      MCMCfit = SVM_MCMC(y, 100000)
      colnames(MCMCfit$theta) = c('sigma^2', 'phi', 'gamma')
      MCMCl = gather(as.data.frame(MCMCfit$theta[50001:100000,]), Variable, Draws)
      MCMCl$Dataset = i
      True = rbind(True, MCMCl)
      
      states = apply(MCMCfit$x[50001:100000,1:(T+1)], 2, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))) %>%
        t() %>% as.data.frame()
      colnames(states) = c('Mean', 'L95', 'U95')
      states$Method = 'MCMC'
      
      weights = c(thetaWeights(y, MCMCfit$theta[50001:100000, ], MCMCfit$x[50001:100000, ], MCMCfit$s[50001:100000, ]))
      initMuTheta = apply(cbind(log(MCMCfit$theta[50001:100000, 1]), MCMCfit$theta[50001:100000, 2:3]), 2, 
                          function(x) sum(x * weights))
      initLTheta  = apply(cbind(log(MCMCfit$theta[50001:100000, 1]), MCMCfit$theta[50001:100000, 2:3]), 2, 
                          function(x) sqrt( sum((x - mean(x))^2 * weights) ))
      initMuX = states$Mean
      initLX = apply(MCMCfit$x[50001:100000, 1:(T+1)], 2, sd)
      #nitMuTheta = c(-1.5, 0, -1)
      initMu = c(initMuTheta, initMuX)
      initL = diag(c(initLTheta, initLX))
      
      VBfit = SGA_SVM(y=y, M=50, maxIter=5000, Mu=initMu, L=initL, meanfield=meanfield, alpha=0.1)
      
      VBsd = sqrt(diag(VBfit$L %*% t(VBfit$L)))
      dSigmaSq = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
      dPhi = dnorm(supportReal, VBfit$Mu[2], VBsd[2])
      dGamma = dnorm(supportReal, VBfit$Mu[3], VBsd[3])
      dsStart = dlnorm(supportPos, initMuTheta[1], initLTheta[1])
      dpStart = dnorm(supportReal, initMuTheta[2], initLTheta[2])
      dgStart = dnorm(supportReal, initMuTheta[3], initLTheta[3])
      
      vbstates = data.frame(mu=VBfit$Mu[4:(T+4)], sd=VBsd[4:(T+4)]) %>%
        apply(1, function(x) c(x[1], x[1]-1.96*x[2], x[1]+1.96*x[2])) %>% t() %>% as.data.frame()
      colnames(vbstates) = c('Mean', 'L95', 'U95')
      vbstates$Method = 'ADVI'
      st = rbind(states, vbstates)
      st$t = 0:T
      st$Dataset = i
      st$Actual = c(x0, x)
      Xt = rbind(Xt, st)
      
      df = data.frame(c(rep(supportPos, 1), rep(supportReal, 2)),
                      c(dSigmaSq, dPhi, dGamma, dsStart, dpStart, dgStart), 
                      rep(c('sigma^2', 'phi', 'gamma'), rep(500, 3)), 
                      rep(c('converged', 'starting'), rep(1500, 2)), i)
      VB = rbind(VB, df)
    }
    colnames(VB) = c('Support', 'Density', 'Variable', 'Version', 'Dataset')
    True <- True %>% group_by(Variable, Dataset) %>%
      filter(Draws > quantile(Draws, 0.025) & Draws < quantile(Draws, 0.975)) %>% ungroup()
    return(list(ADVI=VB, True=True, Xt=Xt))  
  }
}
