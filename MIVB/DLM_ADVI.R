library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
sourceCpp("DLM_MCMC.cpp")
sourceCpp("DLMSplit.cpp")

# Parameters
gamma = 0.5
phi = 0.8
sigmaSqY = 1
sigmaSqX = 1

T = 250
MCMCreps = 10000
reps = 1
xderiv = FALSE
supportReal = seq(-3, 3, length.out=500)
supportPos = seq(0.01, 3, length.out=500)

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
  colnames(MCMCfit$theta) = c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma')
  if(xderiv){
    colnames(MCMCfit$x) = c(paste0('X[',0:(T-1), ']'), 'X[T]', 'aTT', 'pTT')
    MCMCfit$x = as.data.frame(MCMCfit$x)
    MCMCl = gather(as.data.frame(cbind(MCMCfit$theta, select(MCMCfit$x, 1, T+1))), variable, draw)
  } else {
    MCMCl = gather(as.data.frame(MCMCfit$theta), variable, draw)
  }
  MCMCl$dataset = i
  MCMCdf = rbind(MCMCdf, MCMCl)
  
  initMuTheta = colMeans(cbind(log(MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 1:2]), MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 3:4]))
  initLTheta = apply(cbind(log(MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 1:2]), MCMCfit$theta[(MCMCreps/2 + 1):MCMCreps, 3:4]), 2, sd)
  initMuX = colMeans(MCMCfit$x[(MCMCreps/2 + 1):MCMCreps, 1:(T+1)])
  initLX = apply(MCMCfit$x[(MCMCreps/2 + 1):MCMCreps, 1:(T+1)], 2, sd)
  
  initMuTheta[1:3] = c(-1.5, -1.5, 0)
  #initLTheta[3:4] = c(0.3, 1)
  
  initMu = c(initMuTheta, initMuX)
  initL = diag(c(initLTheta, initLX))
  
  VBfit = SGA_DLM(y, 50, 5000, initMu, initL, meanfield=TRUE, xderiv=xderiv, alpha=0.1)
  
  VBsd = sqrt(diag(VBfit$L %*% t(VBfit$L)))
  dSigmaSqY = dlnorm(supportPos, VBfit$Mu[1], VBsd[1])
  dSigmaSqX = dlnorm(supportPos, VBfit$Mu[2], VBsd[2])
  dPhi = dnorm(supportReal, VBfit$Mu[3], VBsd[3])
  dGamma = dnorm(supportReal, VBfit$Mu[4], VBsd[4])

  dsyStart = dlnorm(supportPos, initMuTheta[1], initLTheta[1])
  dsxStart = dlnorm(supportPos, initMuTheta[2], initLTheta[2])
  dpStart = dnorm(supportReal, initMuTheta[3], initLTheta[3])
  dgStart = dnorm(supportReal, initMuTheta[4], initLTheta[4])
  
  if(xderiv){
    dX0 = dnorm(supportReal, VBfit$Mu[5], VBsd[5])
    dXT = dnorm(supportReal, VBfit$Mu[T+5], VBsd[T+5])
    dX0Start = dnorm(supportReal, initMuX[1], initLX[1])
    dXTStart = dnorm(supportReal, initMuX[T+1], initLX[T+1])
  
    df = data.frame(c(rep(supportPos, 2), rep(supportReal, 4)),
                  c(dSigmaSqY, dSigmaSqX, dPhi, dGamma, dX0, dXT, dsyStart, dsxStart, dpStart, dgStart, dX0Start, dXTStart), 
                  rep(c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma', 'X[0]', 'X[T]'), rep(500, 6)), 
                  rep(c('converged', 'starting'), rep(3000, 2)), i)
  } else {
    df = data.frame(c(rep(supportPos, 2), rep(supportReal, 2)),
                    c(dSigmaSqY, dSigmaSqX, dPhi, dGamma, dsyStart, dsxStart, dpStart, dgStart), 
                    rep(c('sigma[y]^2', 'sigma[x]^2', 'phi', 'gamma'), rep(500, 4)), 
                    rep(c('converged', 'starting'), rep(2000, 2)), i)
    
  }
  VB = rbind(VB, df)
}

colnames(MCMCdf) = c('variable', 'draws', 'dataset')
colnames(VB) = c('support', 'density', 'variable', 'lambda', 'dataset')

MCMCfil <- MCMCdf %>% group_by(variable, dataset) %>%
  filter(draws > quantile(draws, 0.005) & draws < quantile(draws, 0.995))

ggplot() + geom_density(data=MCMCfil, aes(x=draws)) + 
  geom_line(data=VB, aes(x=support, y=density, colour=lambda)) +
  facet_grid(variable~dataset, scales='free', labeller=label_parsed) + 
  theme(strip.text.x = element_blank(), strip.text.y = element_text(angle=0),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.title=element_blank())

  
     
  