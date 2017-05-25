library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(truncnorm)
sourceCpp('AR1Mean.cpp')
logPhiDens = function(x, phi, sigmaSq, gamma){
  term1 = 1/2 * log(1-phi^2)
  term2 = -1/(2*sigmaSq)
  term3 = sum(phi^2*(x[1:T]-gamma)^2 - 2*phi*(x[2:(T+1)]-gamma)*(x[1:T]-gamma))
  term4 = - (1-phi^2)*(x[1]-gamma/(1-phi))^2
  return(term1 + term2 * (term3 + term4))
}
logGammaDens = function(x, phi, sigmaSq, gamma){
  term1 = -gamma^2 / 200
  term2 = -1/(2*sigmaSq) * sum((x[2:(T+1)]-gamma-phi*(x[1:T]-gamma))^2)
  term3 = -(1-phi^2)/(2*sigmaSq) * (x[1] - gamma/(1-phi))^2
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
                       + 1/2 * (1 - phiD[i-1]^2) * (x[1]-gamma/(1-phi))^2)
    
    candidate = rtruncnorm(1, -1, 1, phiD[i-1], 0.1)
    canQDens = log(dtruncnorm(candidate, -1, 1, phiD[i-1], 0.1))
    oldQDens = log(dtruncnorm(phiD[i-1], -1, 1, candidate, 0.1))
    canPDens = logPhiDens(x, candidate, sigD[i], gammaD[i-1])
    oldPDens = logPhiDens(x, phiD[i-1], sigD[i], gammaD[i-1])
    ratio = min(1, exp(canPDens - oldPDens - canQDens + oldQDens))
    if(runif(1) < ratio){
      phiD[i] = candidate
      acceptP = acceptP + 1
    } else {
      phiD[i] = phiD[i-1]
    }
    
    candidate = rnorm(1, gammaD[i-1], 0.1)
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
                    rep(c('SigmaSq', 'Phi', 'Gamma'), rep(reps/2, 3)),
                    dataset))
}

T = 500
sigmaSq = 1
phi = 0.9
gamma = 2
reps = 4

VB = data.frame()
MCMCdf = data.frame()

for(i in 1:reps){
  x = rep(0, T+1)
  x[1] = rnorm(1, gamma/(1-phi), sqrt(sigmaSq / (1-phi^2)))
  for(t in 2:(T+1)){
    x[t] = rnorm(1, gamma + phi*(x[t-1]-gamma), sqrt(sigmaSq))
  }
  MCMCdf = rbind(MCMCdf, MCMC(x, 10000, i))
  
  initMu = rep(0, 0, 0)
  initL = diag(0.1, 3)
  VBfit = SGA_AR1(x, 25, 5000, initMu, initL)
  
  supportPhi = seq(-0.99, 0.99, length.out=500)
  supportSigmaSq = seq(0.01, 3, length.out=500)
  supportGamma = seq(1, 5, length.out=500)
  VBsd = diag(sqrt(VBFit$L %*% t(VBFit$L)))
  dSigmaSq = dlnorm(supportSigmaSq, VBfit$Mu[1], VBsd[1])
  dPhi = dnorm(supportPhi, VBfit$Mu[2], VBsd[2])
  dGamma = dnorm(supportGamma, VbfiT$Mu[3], VBsd[3])
  
  df = data.frame(c(supportSigmaSq, supportPhi, supportGamma), c(dSigmaSq, dPhi, dGamma), 
                  rep(c('SigmaSq', 'Phi', 'Gamma'), rep(500, 3)),
                  i)
  VB = rbind(VB, df)
}

colnames(MCMCdf) = c('draws', 'variable', 'dataset')
colnames(VB) = c('support', 'density', 'variable', 'dataset')

ggplot() + geom_density(data=MCMCdf, aes(x=draws)) + 
  geom_line(data=VB, aes(x=support, y=density), colour='red') +
  facet_grid(dataset ~ variable)