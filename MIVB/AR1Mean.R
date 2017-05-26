library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(msm)
sourceCpp('AR1Mean.cpp')

logPhiDens = function(x, phi, sigmaSq, gamma){
  term1 = 1/2 * log(1-phi^2)
  term2 = -1/(2*sigmaSq)
  term3 = sum((x[1:T]-gamma - phi*(x[2:(T+1)]-gamma))^2)
  term4 = (1-phi^2)*(x[1]-gamma)^2
  return(term1 + term2 * (term3 + term4))
}
logGammaDens = function(x, phi, sigmaSq, gamma){
  term1 = -gamma^2 / 200
  term2 = -1/(2*sigmaSq) * sum((x[2:(T+1)]-gamma-phi*(x[1:T]-gamma))^2)
  term3 = -(1-phi^2)/(2*sigmaSq) * (x[1] - gamma)^2
  return(term1 + term2 + term3)
}
MCMC = function(x, reps=10000, dataset){
  phiD = rep(0, reps)
  sigD = rep(0, reps)
  gammaD = rep(0, reps)
  acceptP = 0
  acceptG = 0
  for(i in 2:reps){
    sigD[i] = 1/rgamma(1, (T+3)/2, 1 + 1/2 * sum((x[2:(T+1)] - gammaD[i-1] - phiD[i-1]*(x[1:T]-gammaD[i-1]))^2)
                       + 1/2 * (1 - phiD[i-1]^2) * (x[1]-gamma)^2)
    
    candidate = rtnorm(1, phiD[i-1], 0.08, -1, 1)
    canQDens = dtnorm(candidate, phiD[i-1], 0.08, -1, 1, log=TRUE)
    oldQDens = dtnorm(phiD[i-1], candidate, 0.08, -1, 1, log=TRUE)
    canPDens = logPhiDens(x, candidate, sigD[i], gammaD[i-1])
    oldPDens = logPhiDens(x, phiD[i-1], sigD[i], gammaD[i-1])
    ratio = min(1, exp(canPDens - oldPDens - canQDens + oldQDens))
    if(runif(1) < ratio){
      phiD[i] = candidate
      acceptP = acceptP + 1
    } else {
      phiD[i] = phiD[i-1]
    }
    
    candidate = rnorm(1, gammaD[i-1], 0.5)
    canPDens = logGammaDens(x, phiD[i], sigD[i], candidate)
    oldPDens = logGammaDens(x, phiD[i], sigD[i], gammaD[i-1])
    ratio = min(1, exp(canPDens - oldPDens))
    if(runif(1) < ratio){
      gammaD[i] = candidate
      acceptG = acceptG + 1
    } else {
      gammaD[i] = gammaD[i-1]
    }
  }
  print(acceptP/reps)
  print(acceptG/reps)
  return(data.frame(c(sigD[(reps/2+1):reps], phiD[(reps/2+1):reps], gammaD[(reps/2+1):reps]), 
                    rep(c('sigma^2', 'phi', 'gamma'), rep(reps/2, 3)),
                    dataset))
}

T = 500
sigmaSq = 1
phi = 0.8
gamma = 0.5
reps = 6

VB = data.frame()
MCMCdf = data.frame()
for(i in 1:reps){
  x = rep(0, T+1)
  x[1] = rnorm(1, gamma, sqrt(sigmaSq / (1-phi^2)))
  for(t in 2:(T+1)){
    x[t] = rnorm(1, gamma + phi*(x[t-1]-gamma), sqrt(sigmaSq))
  }
  MCMCfit = MCMC(x, 2000, i)
  MCMCdf = rbind(MCMCdf, MCMCfit)
  smu = -1#mean(log(MCMCfit[1:1000, 1]))
  ssd = 0.5#sd(log(MCMCfit[1:1000, 1]))
  pmu = 0#mean(MCMCfit[1001:2000, 1])
  psd = 0.3#sd(MCMCfit[1001:2000, 1])
  gmu = 0#mean(MCMCfit[2001:3000, 1])
  gsd = 1#sd(MCMCfit[2001:3000, 1])
  
  initMu = c(smu, pmu, gmu)
  initL = diag(c(ssd, psd, gsd))
  VBfit = SGA_AR1M(x, 100, 5000, initMu, initL)
  
  supportPhi = seq(-0.8, 1.05, length.out=500)
  supportSigmaSq = seq(0.01, 2, length.out=500)
  supportGamma = seq(-1.5, 2, length.out=500)
  VBsd = sqrt(diag(VBfit$L %*% t(VBfit$L)))
  dSigmaSq = dlnorm(supportSigmaSq, VBfit$Mu[1], VBsd[1])
  dPhi = dnorm(supportPhi, VBfit$Mu[2], VBsd[2])
  dGamma = dnorm(supportGamma, VBfit$Mu[3], VBsd[3])
  dsStart = dlnorm(supportSigmaSq, smu, ssd)
  dpStart = dnorm(supportPhi, pmu, psd)
  dgStart = dnorm(supportGamma, gmu, gsd)
  
  df = data.frame(c(supportSigmaSq, supportPhi, supportGamma), c(dSigmaSq, dPhi, dGamma, dsStart, dpStart, dgStart), 
                  rep(c('sigma^2', 'phi', 'gamma'), rep(500, 3)), rep(c('converged', 'starting'), rep(1500, 2)), i)
  VB = rbind(VB, df)
}

colnames(MCMCdf) = c('draws', 'variable', 'dataset')
colnames(VB) = c('support', 'density', 'variable', 'lambda', 'dataset')

ggplot() + geom_density(data=MCMCdf, aes(x=draws)) + 
  geom_line(data=VB, aes(x=support, y=density, colour=lambda)) +
  facet_grid(variable ~ dataset, scales='free', labeller=label_parsed) + 
  theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank())
