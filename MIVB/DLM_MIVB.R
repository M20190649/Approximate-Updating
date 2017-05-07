library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(MASS)
library(pscl)
library(distr)
sourceCpp("DLM_MCMC.cpp")
source('qSelection.R')
#source('DLM_CopVB.R')

# CHANGE NORMAL MIXTURE TO T(loc, scale) IN FITREALMARGINALS AND MIVB

SimDLM = function(T, S, h, sigmaSqY=1, sigmaSqX=2, phi=0.95, gamma=2){
  x0 = rnorm(1, 0, sqrt(sigmaSqX/(1-phi^2)))
  x = rep(0, T+S+h)
  y = rep(0, T+S+h)
  for(t in 1:(T+S+h)){
    if(t == 1){
      x[t] = phi*x0 + rnorm(1, 0, sqrt(sigmaSqX)) 
    } else {
      x[t] = phi*x[t-1] + rnorm(1, 0, sqrt(sigmaSqX))
    }
    y[t] = mu + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
  }
  return(y)
}

MIVB = function(T, S, h, sigmaSqY=1, sigmaSqX=2, phi=0.95, gamma=2, MCMCreps=T*10){

  # Simulate Data
  y = SimDLM(T, S, h, sigmaSqY, sigmaSqX, phi, gamma)
  
  # Run MCMC
  MCMC = DLM_MCMC(y[1:T], MCMCreps)
  
  # Filter MCMC forward and calculate logscores
  ydens = matrix(0, h, 1000)
  for(i in 1:500){
    u = sample((MCMCreps/2+1):MCMCreps, 1)
    gammaDraw = MCMC$theta[u, 4]
    phiDraw = MCMC$theta[u, 3]
    sigxDraw = MCMC$theta[u, 2]
    sigyDraw = MCMC$theta[u, 1]
    xTmean = MCMC$x[u, T+2]
    xTvar = MCMC$x[u, T+3]
    XTS = FFUpdatercpp(y[(T+1):(T+S)], phiDraw, gammaDraw, sigyDraw, sigxDraw, xTmean, xTvar)
    ydens[,1] = ydens[,1] + dnorm(ysupport, gammaDraw + phiDraw*XTS[1], sqrt(sigyDraw + phiDraw^2 * XTS[2] + sigxDraw))/500
    if(h > 1){
      for(t in 2:h){
        XTS = FFUpdatercpp(y[T+S+t-1], phiDraw, gammaDraw, sigyDraw, sigxDraw, XTS[1], XTS[2])
        ydens[,t] = ydens[,t] + dnorm(ysupport, gammaDraw + phiDraw*XTS[1], sqrt(sigyDraw + phiDraw^2 * XTS[2] + sigxDraw))/500
      }
    }
  }
  obs = rep(0, h)
  logscoreMCMC = rep(0, h)
  for(t in 1:h){
    obs[t] = min(which((y[T+S+t] < ysupport) == TRUE))
    logscoreMCMC[t] = log(ydens[obs, t])
  }
  
  # Fit Marginals
  # Theta
  SigSqYFit = FitPositiveMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 1])
  SigSqXFit = FitPositiveMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 2])
  PhiFit = FitRealMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 3])
  GammaFit = FitRealMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 4])
  thetaDist = c(SigSqYFit$AIC, SigSqXFit$AIC, PhiFit$AIC, GammaFit$AIC)
  thetaParams = c(SigSqYFit$params, SigSqXFit$params, PhiFit$params, GammaFit$params)
  
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
    xParams = matrix(0, T+1, 2)
    xParams[,1] = apply(MCMC$x[(MCMCreps/2+1):MCMCreps,1:(T+1)], 2, mean)
    xParams[,2] = apply(MCMC$x[(MCMCreps/2+1):MCMCreps,1:(T+1)], 2, var)
  } else {
    xParams = matrix(0, T, 3)
    for(i in 1:T){
      fitT = fitdistr(MCMC$x[(MCMCreps/2+1):MCMCreps,i],
                      't', list(m=mean(x), s=var(x), df=5))
      Xparams[i, ] = fitT$estimate
    }
  }
  
  # Select the Vine Structure
  subset = min(15, T/10)
  Vine = RestrictedVineSelect(T, subset, 4, MCMCreps, MCMC)
  
  # Estimate the Full Vine Parameters
  v = sample((MCMCreps/2+1):MCMCreps, min(1000, MCMCreps/2))
  pobtheta = pobs(MCMC$theta[v,])
  pobx = pobs(MCMC$x[v, 0:(T+1)])
  VineObject = RVineMatrix(Vine$Matrix, Vine$family)
  VineFit = RVineSeqEst(cbind(pobtheta, pobx), VineObject)
  
  # Run the VB updater to time T+S
  VBfit = CopulaVB(y, S, thetaDist, xDist, thetaParams, xParams, Vine, M, maxIter)
  unifs = RVineSim(VBfit$Vine, 500)
  simul = MarginalTransform(unifs, thetaDist, Dist, VBfit$thetaParams, VBfit$xParams)
  ydensVB = matrix(0, h, 1000)
  
  # Calculate Density
  for(i in 1:500){
    ydensVB[,1] = ydensVB[,1] + dnorm(ysupport, simul[i, 4] + simul[i, 3]*simul[i,T+S+5], 
                                    sqrt(simul[i,1] + simul[i,2]))/500
  }
   
  logscoreVB = rep(0, h)
  logscoreVB[1] = log(ydensVB[obs[1], 1])
  
  # Update for forcasting y_{T+S+2:T+S+h}
  if(h > 1){
    for(t in 2:h){
      VBfit = CopulaVB(y[T+S+t], 1, thetaDist, xDist, VBfit$thetaParams, VBfit$xParams, VBfit$Vine, M, maxIter)
      unifs = RVineSim(VBfit$Vine, 500)
      simul = MarginalTransform(unifs, thetaDist, Dist, VBfit$thetaParams, VBfit$xParams)
      for(i in 1:500){
        ydensVB[,t] = ydensVB[,t] + dnorm(ysupport, simul[i, 4] + simul[i, 3]*simul[i,T+S+5], 
                                        sqrt(simul[i,1] + simul[i,2]))/500
        
      }
      logscoreVB[2] = log(ydensVB[obs[t], t])
    }
  }
  
  return(data.frame(futureY=ysupport[obs], logscoreMCMC, logscoreVB))
}
