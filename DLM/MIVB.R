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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MIVB = function(T, S, h, sigmaSqY=1, sigmaSqX=2, phi=0.95, gamma=2, MCMCreps=T*10){

  # Simulate Data
  y = SimDLM(T, S, h, sigmaSqY, sigmaSqX, phi, gamma)
  
  # Run MCMC
  MCMC = DLM_MCMC(y[1:T], MCMCreps)
  
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
  } else { # Hoping this one doesn't happen. Would rather choose between Normal and T
    Xparams = matrix(0, T, 5)
    for(i in 1:T){
      Mixfit = normalmixEM(MCMC$x[(MCMCreps/2+1):MCMCreps,i], k = 2)
      Xparams[i, 1:2] = Mixfit$mu
      Xparams[i, 3:4] = Mixfit$sigma^2
      Xparams[i, 5] = Mixfit$lambda[1]
    }
  }
  
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
  
  # Find where to attach X using restricted Dissmann's Algorithm
  tree1tau = TauMatrix(cbind(pobtheta, pobx))
  Tree1edge = which.max(rowSums(abs(tree1tau[1:4, 1:subset+4])))
  Vine$Matrix[subset+4, 1:subset] = Tree1edge
  if(thetavine$type == 'D-vine' & Tree1edge %in% thetaVine$Matrix[1:2,1]){
    if(Tree1edge == thetaVine$Matrix[1,1]){
      Vine$Matrix[1:3+subset,1:subset] = thetaVine$Matrix[1:3,1]
    } else {
      Vine$Matrix[1:3+subset,1:subset] = thetaVine$Matrix[c(3,4,1),1]
    }
    XCop = RVineCopSelect(cbind(pobtheta, pobx), Matrix=Vine$Matrix, indeptest=TRUE, level=0.01, trunclevel=5, cores=4)
    Vine$Family[subset+4,1:subset] = Mode(XCop$family[subset+4,1:subset]) 
    Vine$Family[subset+3,1:subset] = Mode(XCop$family[subset+3,1:subset])
    Vine$Family[subset+2,1:subset] = Mode(XCop$family[subset+2,1:subset])
    Vine$Family[subset+1,1:subset] = Mode(XCop$family[subset+1,1:subset])
    Vine$Family[subset,1:(subset-1)] = Mode(XCop$family[subset,1:(subset-1)])
  } else {
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

    edges = c(thetaVine$Matrix[4,which(diag(thetaVine$Matrix)[1:3]==Tree1edge)],
                   diag(thetaVine$Matrix)[which(thetaVine$Matrix[4,1:3]==Tree1edge)])
    pobthetaT2 = matrix(0, length(v), length(edges))
    
    positions = matrix(0, length(edges), 2)
    counter=1
    if(any(diag(thetaVine$Matrix)[1:3]==Tree1edge)){
      positions[1,] = c(4,which(diag(thetaVine$Matrix)[1:3]==Tree1edge))
      counter=2
    }
    if(any(thetaVine$Matrix[4,1:3]==Tree1edge)){
      mat = matrix(which(thetaVine$Matrix[4,1:3]==Tree1edge),ncol=1)
      mat = cbind(4, mat)
      positions[counter:length(edges),] = mat
    }
    
    for(i in 1:length(edges)){
      pobthetaT2[,i] = BiCopHfunc(pobtheta[,edges[i]], pobtheta[,Tree1edge], family=thetaVine$family[positions[i,1], positions[i,2]],
                                  par=thetaVine$par[positions[i,1], positions[i,2]], par2=thetaVine$par2[positions[i,1], positions[i,2]])$hfunc2
    }
    
    #Next step: Calculate Tau and pick edge for second tree
    #Should add early check to see if length(edges) = 1
  } 
  
  
  
  
}