library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(MASS)
library(pscl)
library(distr)
library(mixtools)
sourceCpp("DLM_MCMC.cpp")

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

FitRealMarginals = function(x){
  n = length(x)
  AIC = rep(0, 2)
  # Normal
  NormMu = mean(x)
  NormVar = (n-1) / n * var(x)
  AIC[1] = 4  - 2*sum(dnorm(x, NormMu, sqrt(NormVar), log=TRUE))
  # Student T
  #fitT = fitdistr(x, 't', list(m=mean(x), s=var(x), df=5))
  #AIC[2] = 6 - 2*fitT$loglik
  # NormMix
  fitMix = normalmixEM(x, k = 2)
  AIC[2] = 10 - 2*fitMix$loglik
  if(AIC[1] < AIC[2]){
    return(list(Dist=1, params=c(NormMu, NormVar)))
  } else {
    return(list(Dist=2), params=c(fitMix$mu, fitMix$sigma^2, fitMix$lambda[1]))
  }
}

FitPositiveMarginals = function(x){
  AIC = rep(0, 4)
  # Weibull
  fitW = fitdistr(x, 'weibull')
  AIC[1] = 4 - 2*fitW$loglik
  # Gamma
  fitG = fitdistr(x, 'gamma')
  AIC[2] = 4 - 2*fitG$loglik
  # Inverse Gamma
  mu = mean(x)
  v = var(x)
  a = mu^2/v + 2
  b = mu*(mu^2/v + 1)
  AIC[3] = 4 - 2*sum(log(densigamma(x, a, b)))
  # Log Normal
  fitLN = fitdistr(x, 'lognormal')
  AIC[4] = 4 - 2*fitLN$loglik
  distChoice = which.min(AIC)
  if(distChoice==1){
    return(list(Dist=1, params=fitW$estimate))
  } else if(distChoice==2){
    return(list(Dist=2, params=fitG$estimate))
  } else if(distChoice==3){
    return(list(Dist=3, params=c(a, b)))
  } else {
    return(list(Dist=4, params=fitLN$estimate))
  }
}

Mode = function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

RestrictedVineSelect = function(T, subset, MCMCreps, MCMC){
  # Create Vine Matrix for X subset
  Vine = list(Matrix = matrix(0, subset+4, subset+4), 
              Family = matrix(0, subset+4, subset+4),
              Par1 = matrix(0, subset+4, subset+4),
              Par2 = matrix(0, subset+4, subset+4))
  diag(Vine$Matrix)[1:subset] = 1:subset+4
  for(i in 2:(subset)){
    for(j in 1:(i-1)){
      Vine$Matrix[i, j] = subset + j - i + 5
    }
  }
  
  # Fill in Theta Subtree, fit to a subsample of MCMC draws
  v = sample((MCMCreps/2+1):MCMCreps, min(1000, MCMCreps/2))
  pobtheta = pobs(MCMC$theta[v,])
  pobx = pobs(MCMC$x[v,-subset:0+T+1])
  thetaVine = RVineStructureSelect(pobtheta, cores=4)
  Vine$Matrix[1:4+subset,1:4+subset] = thetaVine$Matrix
  Vine$Family[1:4+subset,1:4+subset] = thetaVine$family
  Vine$Par1[1:4+subset,1:4+subset] = thetaVine$par
  Vine$Par2[1:4+subset,1:4+subset] = thetaVine$par2
  
  # Code for Tree 1 X_t using restricted Dissmann's Algorithm
  tree1tau = TauMatrix(cbind(pobtheta, pobx))
  Tree1edge = which.max(rowSums(abs(tree1tau[1:4, 1:subset+4])))
  Vine$Matrix[subset+4, 1:subset] = Tree1edge
  Tree1fit = rep(0, subset)
  for(i in 1:subset){ 
    Tree1fit[i] = BiCopSelect(pobx[,i], pobtheta[,Tree1edge], indeptest=TRUE, level=0.01)$family
  }
  Vine$Family[subset+4,1:subset] = Mode(Tree1fit)
  pobxT2 = matrix(0, length(v), subset) 
  for(i in 1:subset){
    CopEstim = BiCopEst(pobx[,i], pobtheta[,Tree1edge], family=Mode(Tree1fit))
    pobxT2[,i] = BiCopHfunc(pobx[,i], pobtheta[,Tree1edge], family=Mode(Tree1fit), par=CopEstim$par, par2=CopEstim$par2)$hfunc2
  }
  positions = vector(length=0)
  if(any(diag(thetaVine$Matrix)[1:3]==Tree1edge)){
    positions[1,] = which(diag(thetaVine$Matrix)[1:3]==Tree1edge)
  }
  if(any(thetaVine$Matrix[4,1:3]==Tree1edge)){
    positions = c(positions, which(thetaVine$Matrix[4,1:3]==Tree1edge))
  }
  edges = c(thetaVine$Matrix[4,which(diag(thetaVine$Matrix)[1:3]==Tree1edge)],
            diag(thetaVine$Matrix)[which(thetaVine$Matrix[4,1:3]==Tree1edge)])
  pobthetaT2 = matrix(0, length(v), length(edges))
  for(i in 1:length(edges)){
    pobthetaT2[,i] = BiCopHfunc(pobtheta[,edges[i]], pobtheta[,Tree1edge], family=thetaVine$family[4, positions[i]],
                                par=thetaVine$par[4, positions[i]], par2=thetaVine$par2[4], positions[i])$hfunc2
  }
  
  # Tree 2
  if(length(edges)==1){
    Tree2edge = edges
  } else {
    tree2tau = TauMatrix(cbind(pobthetaT2, pobxT2))
    Tree2edge =  which.max(rowSums(abs(tree2tau[1:length(edges), 1:subset+length(edges)])))
  }
  Vine$Matrix[subset+3, 1:subset] = Tree2edge
  Tree2fit = rep(0, subset)
  for(i in 1:subset){ 
    Tree2fit[i] = BiCopSelect(pobxT2[,i], pobthetaT2[,Tree2edge], indeptest=TRUE, level=0.01)$family
  }
  Vine$Family[subset+3,1:subset] = Mode(Tree2fit)
  pobxT3 = matrix(0, length(v), subset)
  for(i in 1:subset){
    CopEstim = BiCopEst(pobxT2[,i], pobthetaT2[,Tree2edge], family=Mode(Tree2fit))
    pobxT3[,i] = BiCopHfunc(pobxT2[,i], pobthetaT2[,Tree2edge], family=Mode(Tree2fit), par=CopEstim$par, par2=CopEstim$par2)$hfunc2
  }
  edges = vector(length=0)
  positions = vector(length=0)
  if(thetaVine$Matrix[1,1]==Tree2edge){
    edges = thetaVine$Matrix[3,1]
    positions = 1
  } else if(thetaVine$Matrix[3,1]==Tree2edge){
    edges = thetaVine$Matrix[1,1]
    positions = 1
  }
  if(thetaVine$Matrix[2,2]==Tree2edge){
    edges = c(edges, thetaVine$Matrix[3,2])
    positions = c(positions, 2)
  } else if(thetaVine$Matrix[3,2]==Tree2edge){
    edges = c(edges, thetaVine$Matrix[2,2])
    positions = c(positions, 2)
  }
  pobthetaT3 = matrix(0, length(v), length(edges))
  for(i in 1:length(edges)){
    pobthetaT3[,i] =  BiCopHfunc(pobthetaT2[,edges[i]], pobthetaT2[,Tree2edge], family=thetaVine$family[3, positions[i]],
                                 par=thetaVine$par[3, positions[i]], par2=thetaVine$par2[3, positions[i]])$hfunc2
  }
  
  # Tree 3
  if(length(edges)==1) {
    Tree3edge = edges
    Tree4edge = 10 - sum(Tree1edge, Tree2edge, Tree3edge)
  } else {
    tree3tau = TauMatrix(cbind(pobthetaT3, pobxT3))
    Tree3edge = which.max(rowSums(abs(tree3tau[1:length(edges), 1:subset+length(edges)])))
    Tree4edge = 10 - sum(Tree1edge, Tree2edge, Tree3edge)
  }
  Vine$Matrix[1:2 + subset, 1:subset] = c(Tree4edge, Tree3edge)
  Tree3fit = rep(0, subset)
  for(i in 1:subset){
    Tree3fit[i] = BiCopSelect(pobxT3[,i], pobthetaT3[,Tree3edge], indeptest=TRUE, level=0.01)$family
  }
  Vine$Family[subset+2, 1:subset] = Mode(Tree3fit)
  pobxT4 = matrix(0, length(v), subset)
  for(i in 1:subset){
    CopEstim = BiCopEst(pobxT3[,i], pobthetaT3[,Tree3edge], family=Mode(Tree3fit))
    pobxT4[,i] = BiCopHfunc(pobxT3[,i], pobthetaT3[,Tree3edge], family=Mode(Tree3fit), par=CopEstim$par, par2=CopEstim$par2)$hfunc2
  }
  pobthetaT4 = BiCopHfunc(pobthetaT2[,Tree4edge], pobthetaT2[,Tree3edge], family=thetaVine$family[2, 1],
                          par=thetaVine$par[2, 1], par2=thetaVine$par2[2, 1])$hfunc2
  
  # Tree 4
  Tree4fit = rep(0, subset)
  for(i in 1:subset){
    Tree4fit[i] = BiCopSelect(pobxT4[,i], pobthetaT4, indeptest=TRUE, level=0.01)$family
  }
  Vine$Family[subset+1, 1:subset] = Mode(Tree4fit)
  pobxT5 = matrix(0, length(v), subset)
  for(i in 1:subset){
    CopEstim = BiCopEst(pobxT4[,i], pobthetaT4, family=Mode(Tree4fit))
    pobxT5[,i] = BiCopHfunc(pobxT4[,i], pobthetaT4, family=Mode(Tree4fit), par=CopEstim$par, par2=CopEstim$par2)$hfunc2
  }
  
  # Tree 5
  Tree5fit = rep(0, subset-1)
  for(i in 1:(subset-1)){
    Tree5fit[i] = BiCopSelect(pobxT5[,i], pobxT5[,i+1], indeptest=TRUE, level=0.01)$family
  }
  Vine$Family[subset, 1:(subset-1)] = Mode(Tree5fit)
  
  # Extrapolate structure to vine over whole x parameter vector
  FullVine = list(Matrix = matrix(0, T+5, T+5),
                  family = matrix(0, T+5, T+5))
  
  diag(FullVine$Matrix)[1:(T+1)] = 5:(T+5)
  for(i in 2:(T+1)){
    for(j in 1:i-1){
      FullVine$Matrix[i, j] = T + j - i + 6
    }
  }
  
  FullVine$Matrix[(T+2):(T+5),(T+2):(T+5)] = thetavine$Matrix
  FullVine$Matrix[(T+2):(T+5),1:(T+1)] = Vine$Matrix[1:4+subset, 1]
  FullVine$family[(T+2):(T+5),(T+2):(T+5)] = thetavine$family
  FullVine$family[(T+2):(T+5),1:(T+1)] = Vine$family[1:4+subset, 1]
  FullVine$family[T+1, 1:T] = Vine$family[subset, 1]  
  
  return(FullVine)
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
        ydens[,2] = ydens[,2] + dnorm(ysupport, gammaDraw + phiDraw*XTS[1], sqrt(sigyDraw + phiDraw^2 * XTS[2] + sigxDraw))/500
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
  
  # Choose distribution for the last few X marginals
  subset = min(T/20, 25)
  Xfit = list()
  XfitTable = matrix(0, 2, subset)
  for(i in 1:subset){
    Xfit[[i]] = FitRealMarginals(MCMC$x[(MCMCreps/2+1):MCMCreps,T+2-i])
    XfitTable[,i] = Xfit[[i]]$AIC
  }
  Xdist = which.min(rowSums(XfitTable))
  # Fit parameters to whole X set
  if(Xdist == 1){
    Xparams = matrix(0, T+1, 2)
    Xparams[,1] = apply(MCMC$x[(MCMCreps/2+1):MCMCreps,1:(T+1)], 2, mean)
    Xparams[,2] = apply(MCMC$x[(MCMCreps/2+1):MCMCreps,1:(T+1)], 2, var)
  } else { # Hoping this one doesn't happen. Will try to get T estimation stable and use it instead of mixtures
    Xparams = matrix(0, T, 5)
    for(i in 1:T){
      Mixfit = normalmixEM(MCMC$x[(MCMCreps/2+1):MCMCreps,i], k = 2)
      Xparams[i, 1:2] = Mixfit$mu
      Xparams[i, 3:4] = Mixfit$sigma^2
      Xparams[i, 5] = Mixfit$lambda[1]
    }
  }
  
  # Select the Vine Structure
  Vine = RestrictedVineSelect(T, subset, MCMCreps, MCMC)
  
  # Estimate the Full Vine Parameters
  v = sample((MCMCreps/2+1):MCMCreps, min(1000, MCMCreps/2))
  pobtheta = pobs(MCMC$theta[v,])
  pobx = pobs(MCMC$x[v, 0:(T+1)])
  VineObject = RVineMatrix(Vine$Matrix, Vine$family)
  VineFit = RVineSeqEst(cbind(pobtheta, pobx), VineObject)
  
  # Update the VB algorithm
  
}