library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(MASS)
library(pscl)
sourceCpp("DLM_MCMC.cpp")
source('qSelection.R')
source('DLM_CopVB.R')


MIVB = function(T, S, h, sigmaSqY=1, sigmaSqX=1, phi=0.95, gamma=2, MCMCreps=T*10){

  # Simulate Data
  y = SimDLM(T, S, h, sigmaSqY, sigmaSqX, phi, gamma)
  
  # Run MCMC
  MCMC = DLM_MCMC(y[1:T], MCMCreps)
  
  # S + h step ahead forecast
  ysupport = seq(min(y)-sd(y), max(y)+sd(y), length.out=1000)
  obs = min(which((y[T+S+h] < ysupport) == TRUE))
  ydens = rep(0, 1000)
  for(i in 1:500){
    u = sample((MCMCreps/2+1):MCMCreps, 1)
    gammaDraw = MCMC$theta[u, 4]
    phiDraw = MCMC$theta[u, 3]
    sigxDraw = MCMC$theta[u, 2]
    sigyDraw = MCMC$theta[u, 1]
    xTmean = MCMC$x[u, T+2]
    xTvar = MCMC$x[u, T+3]
    phiprod = rep(1, S+h)
    for(i in 2:(S+h)){
      phiprod[i] = phiDraw^2 * phiprod[i-1]
    }
    ydens = ydens + dnorm(ysupport, gammaDraw + phiDraw^(S+h)*xTmean, 
                          sqrt(sigyDraw + sum(sigxDraw*phiprod) + phiDraw^(2*(S+h))*xTvar))/500
    
  }
  logscoreMCMC = log(ydens[obs])
  # Filter MCMC forward to T+S, perform h step ahead forecast
  
  ydensf = rep(0, 1000)
  for(i in 1:500){
    u = sample((MCMCreps/2+1):MCMCreps, 1)
    gammaDraw = MCMC$theta[u, 4]
    phiDraw = MCMC$theta[u, 3]
    sigxDraw = MCMC$theta[u, 2]
    sigyDraw = MCMC$theta[u, 1]
    xTmean = MCMC$x[u, T+2]
    xTvar = MCMC$x[u, T+3]
    XTS = FFUpdatercpp(y[(T+1):(T+S)], phiDraw, gammaDraw, sigyDraw, sigxDraw, xTmean, xTvar)
    phiprod = rep(1, h)
    for(i in 2:(h)){
      phiprod[i] = phiDraw^2 * phiprod[i-1]
    }
    ydensf = ydensf + dnorm(ysupport, gammaDraw + phiDraw^(h)*xTmean, 
                          sqrt(sigyDraw + sum(sigxDraw*phiprod) + phiDraw^(2*(h))*xTvar))/500
  }
  logscoreMCMCf = log(ydensf[obs])
 
  # Fit Marginals
  # Theta
  SigSqYFit = FitPositiveMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 1])
  SigSqXFit = FitPositiveMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 2])
  PhiFit = FitRealMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 3])
  GammaFit = FitRealMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 4])
  thetaDist = c(SigSqYFit$Dist, SigSqXFit$Dist, PhiFit$Dist, GammaFit$Dist)
  thetaParams = matrix(0, 3, 4)
  thetaParams[1:2, 1] = SigSqYFit$params
  thetaParams[1:2, 2] = SigSqXFit$params
  thetaParams[1:length(PhiFit$params), 3] = PhiFit$params
  thetaParams[1:length(GammaFit$params), 4] = GammaFit$params

  # Choose distribution for the last few X marginals
  subset = min(T/10, 15)
  Xfit = list()
  XfitTable = matrix(0, 2, subset)
  for(i in 1:subset){
    Xfit[[i]] = FitRealMarginals(MCMC$x[(MCMCreps/2+1):MCMCreps,T+2-i])
    XfitTable[,i] = Xfit[[i]]$AIC
  }
  xDist = which.min(rowSums(XfitTable))
  # Fit parameters to whole X set
  if(xDist == 1){
    xParams = matrix(0, 2, T+1)
    xParams[1,] = apply(MCMC$x[(MCMCreps/2+1):MCMCreps,1:(T+1)], 2, mean)
    xParams[2,] = apply(MCMC$x[(MCMCreps/2+1):MCMCreps,1:(T+1)], 2, var)
  } else {
    xParams = matrix(0, T, 3)
    for(i in 1:T){
      fitT = fitdistr(MCMC$x[(MCMCreps/2+1):MCMCreps,i], 't')
      Xparams[i, ] = fitT$estimate
    }
  }
  
  # Select and Estimate the Vine
  Vine = RestrictedVineSelect(T, subset, 4, MCMCreps, MCMC)
  
  
  #VineObj = RVineMatrix(Vine$Matrix, Vine$family, Vine$par, Vine$par2)
  # Run the VB updater to time T+S
  VBfit = CopulaVB(y[1:(T+S)], S, thetaDist, xDist, thetaParams, xParams, Vine, 15, 100)

  # Perform h step ahead forecast
  unifs = RVineSim(VBfit$Vine, 500)
  simul = MarginalTransform(unifs, thetaDist, Dist, VBfit$thetaParams, VBfit$xParams)
  ydensVB = rep(0, 1000)
  for(i in 1:500){
    phiprod = rep(1, h)
    for(i in 2:(h)){
      phiprod[i] = phiDraw^2 * phiprod[i-1]
    }
    ydensVB = ydensVB + dnorm(ysupport, simul[i, 4] +  simul[i, 3]^(h)*simul[i,T+S+5],
                             sqrt(simul[i,1] + sum(simul[i,2]*phiprod)))/500
  }
  logscoreVB = log(ydensVB[obs])
  return(data.frame(futureY=y[1:h+T+S], logscoreMCMC, logscoreMCMCf, logscoreVB))
}
