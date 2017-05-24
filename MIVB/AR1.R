library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(truncnorm)
sourceCpp('AR1.cpp')
logPhiDens = function(x, phi, sigmaSq){
  term1 = 1/2 * log(1-phi^2)
  term2 = -1/(2*sigmaSq)
  term3 = phi^2 * (sum(x[1:T]^2) - x[1]^2)
  term4 = -2 * phi * sum(x[1:T]*x[2:(T+1)])
  return(term1 + term2 * (term3 + term4))
}
MCMC = function(x, reps=10000, dataset){
  phiD = rep(0, reps)
  sigD = rep(0, reps)
  accept = 0
  for(i in 2:reps){
    sigD[i] = 1/rgamma(1, (T+3)/2, 1 + 1/2 * sum((x[2:(T+1)] - phiD[i-1]*x[1:T])^2) + 1/2 * (1 - phiD[i-1]^2) * x[1]^2)
    
    candidate = rtruncnorm(1, -1, 1, phiD[i-1], 0.1)
    canQDens = log(dtruncnorm(candidate, -1, 1, phiD[i-1], 0.1))
    oldQDens = log(dtruncnorm(phiD[i-1], -1, 1, candidate, 0.1))
    canPDens = logPhiDens(x, candidate, sigD[i])
    oldPDens = logPhiDens(x, phiD[i-1], sigD[i])
    ratio = min(1, exp(canPDens - oldPDens - canQDens + oldQDens))
    if(runif(1) < ratio){
      phiD[i] = candidate
      accept = accept + 1
    } else {
      phiD[i] = phiD[i-1]
    }
  }
  print(accept/reps)
  return(data.frame(c(sigD[(reps/2+1):reps], phiD[(reps/2+1):reps]), 
                    rep(c('SigmaSq', 'Phi'), rep(reps/2, 2)),
                    dataset))
}

T = 500
sigmaSq = 1
phi = 0.9
reps = 4

VB = data.frame()
MCMCdf = data.frame()

for(i in 1:reps){
  x = rep(0, T+1)
  x[1] = rnorm(1, 0, sqrt(sigmaSq / (1-phi^2)))
  for(t in 2:(T+1)){
    x[t] = rnorm(1, phi*x[t-1], sqrt(sigmaSq))
  }
  MCMCdf = rbind(MCMCdf, MCMC(x, 10000, i))
  
  initMu = rep(0, 0)
  initL = diag(0.1, 2)
  VBfit = SGA_AR1(x, 25, 5000, initMu, initL)
  
  supportPhi = seq(-0.99, 0.99, length.out=500)
  supportSigmaSq = seq(0.01, 3, length.out=500)
  dPhi = dnorm(supportPhi, VBfit$Mu[2], sqrt(sum(VBfit$L[2, ]^2)))
  dSigmaSq = dlnorm(supportSigmaSq, VBfit$Mu[1], Vbfit$L[1, 1])
  
  df = data.frame(c(supportSigmaSq, supportPhi), c(dSigmaSq, dPhi), 
                  rep(c('SigmaSq', 'Phi'), rep(500, 2)),
                  i)
  VB = rbind(VB, df)
}

colnames(MCMCdf) = c('draws', 'variable', 'dataset')
colnames(VB) = c('support', 'density', 'variable', 'dataset')

ggplot() + geom_density(data=MCMCdf, aes(x=draws)) + 
  geom_line(data=VB, aes(x=support, y=density), colour='red') +
   facet_grid(dataset ~ variable)
