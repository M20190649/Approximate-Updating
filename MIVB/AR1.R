library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(msm)
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
    
    candidate = rtnorm(1, phiD[i-1], 0.1, -1, 1)
    canQDens = dtnorm(candidate, phiD[i-1], 0.1, -1, 1, log=TRUE)
    oldQDens = dtnorm(phiD[i-1], candidate, 0.1, -1, 1, log=TRUE)
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
                    rep(c('sigma^2', 'phi'), rep(reps/2, 2)),
                    dataset))
}

T = 500
sigmaSq = 1
phi = 0.8
reps = 6
VB = data.frame()
MCMCdf = data.frame()

for(i in 1:reps){
  x = rep(0, T+1)
  x[1] = rnorm(1, 0, sqrt(sigmaSq / (1-phi^2)))
  for(t in 2:(T+1)){
    x[t] = rnorm(1, phi*x[t-1], sqrt(sigmaSq))
  }
  MCMCfit = MCMC(x, 2000, i)
  MCMCdf = rbind(MCMCdf, MCMCfit)
  smu = -1#mean(log(MCMCfit[1:1000, 1]))
  ssd = 0.5#sd(log(MCMCfit[1:1000, 1]))
  pmu = 0#mean(MCMCfit[1001:2000, 1])
  psd = 0.5#sd(MCMCfit[1001:2000, 1])
  
  initMu = c(smu, pmu)
  initL = diag(c(ssd, psd))
  VBfit = SGA_AR1(x, 100, 5000, initMu, initL, alpha=0.1)
  VBsd = sqrt(diag(VBfit$L %*% t(VBfit$L)))
  supportPhi = seq(-1.2, 1.2, length.out=500)
  supportSigmaSq = seq(0.01, 1.5, length.out=500)
  dSigmaSq = dlnorm(supportSigmaSq, VBfit$Mu[1], VBsd[1])
  dPhi = dnorm(supportPhi, VBfit$Mu[2], VBsd[2])
  dsStart = dlnorm(supportSigmaSq, smu, ssd)
  dpStart = dnorm(supportPhi, pmu, psd)
  
  df = data.frame(c(supportSigmaSq, supportPhi), c(dSigmaSq, dPhi, dsStart, dpStart), 
                  rep(c('sigma^2', 'phi'), rep(500, 2)), rep(c('converged', 'starting'), rep(1000, 2)), i)
  VB = rbind(VB, df)
}

colnames(MCMCdf) = c('draws', 'variable', 'dataset')
colnames(VB) = c('support', 'density', 'variable', 'lambda', 'dataset')

ggplot() + geom_density(data=MCMCdf, aes(x=draws)) + 
  geom_line(data=VB, aes(x=support, y=density, colour=lambda)) +
  facet_grid(variable ~ dataset, scales='free', labeller=label_parsed) + 
  theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank())
