library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(MASS)
library(pscl)
library(distr)
library(mixtools)
sourceCpp("DLM_MCMC.cpp")


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

FitRealMarginals = function(x){
  n = length(x)
  AIC = rep(0, 2)
  # Normal
  NormMu = mean(x)
  NormVar = (n-1) / n * var(x)
  AIC[1] = 4  - 2*sum(dnorm(x, NormMu, sqrt(NormVar), log=TRUE))
  # Student T - MAKE THIS WORK
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

RestrictedVineSelect = function(T, subset, k, MCMCreps, MCMC){
  # Create Vine Matrix for X subset
  Vine = list(Matrix = matrix(0, subset+k, subset+k), 
              Family = matrix(0, subset+k, subset+k))
  diag(Vine$Matrix)[1:subset] = 1:subset+k
  for(i in 2:(subset)){
    for(j in 1:(i-1)){
      Vine$Matrix[i, j] = subset + j - i + k + 1
    }
  }
  
  # Fill in Theta Subtree, fit to a subsample of MCMC draws
  v = sample((MCMCreps/2+1):MCMCreps, min(500, MCMCreps/2))
  pobtheta = pobs(MCMC$theta[v,])
  pobx = pobs(MCMC$x[v,-subset:0+T+1])
  thetaVine = RVineStructureSelect(pobtheta)
  Vine$Matrix[1:k+subset,1:k+subset] = thetaVine$Matrix
  Vine$Family[1:k+subset,1:k+subset] = thetaVine$family

  # Misc Variables to initialise
  TreeEdge = rep(0, k)
  TreeFit = rep(0, subset)
  TreeCop - rep(0, k)
  xCond = vector(length=0)
  candidates = 1:k
  # Sequential selection of edge, selection of copula family and calculation of data for next tree
  for(tree in k:1){
    # Select next edge
    if(tree == k){
      treeTau = TauMatrix(cbind(pobtheta, pobx))
      TreeEdge[tree] = candidates[which.max(rowSums(abs(treeTau[1:k, 1:subset+k])))]
    } else {
      sameCond = which(pobthetaLabels[2,] == as.numeric(paste0(sort(xCond), collapse='')))
      candidates = pobthetaLabels[1, sameCond]
      if(length(candidates)==1){
        TreeEdge[tree] = candidates
      } else {
        treeTau = TauMatrix(cbind(pobtheta[,sameCond], pobx))
        TreeEdge[tree] = candidates[which.max(rowSums(abs(treeTau[1:length(candidates), 1:subset+length(candidates)])))]
      }
    }
    
    # Select copula family for that edge
    for(i in 1:subset){ 
      TreeFit[i] = BiCopSelect(pobx[,i], pobtheta[,TreeEdge[tree]], indeptest=TRUE, level=0.01)$family
    }
    TreeCop[tree] = Mode(TreeFit)
    
    # Calculate next level of xt | theta
    pobxNext = matrix(0, length(v), subset) 
    for(i in 1:subset){
      if(tree == k){
        whichTheta = TreeEdge[tree]
      } else {
        whichTheta = which(apply(pobthetaLabels, 2, function(x) x[1] == TreeEdge[tree] &
                                   x[2] == as.numeric(paste0(sort(xCond), collapse=''))))
      }
      CopEstim = BiCopEst(pobx[,i], pobtheta[,whichTheta], family=TreeCop[tree])
      pobxNext[,i] = BiCopHfunc(pobx[,i], pobtheta[,whichTheta], family=TreeCop[tree], par=CopEstim$par, par2=CopEstim$par2)$hfunc2
    }
    pobx = pobxNext
    
    # Calculate next level of theta_i | theta_(!i)
    if(tree > 1){
      xCond = c(xCond, TreeEdge[tree])
      pobthetaNext = matrix(0, length(v), (tree-1)*2)
      pobthetaLabels = matrix(0, 2, (tree-1)*2)
      for(i in 1:(tree-1)){
        Var1 = diag(thetaVine$Matrix)[i]
        Var2 = thetaVine$Matrix[tree, i]
        Hfunc = BiCopHfunc(pobtheta[,Var1], pobtheta[,Var2], family=thetaVine$family[tree, i],
                   par=thetaVine$par[tree, i], par2=thetaVine$par2[tree, i])
        pobthetaNext[,2*i-1] = Hfunc$hfunc1
        pobthetaNext[,2*i] = Hfunc$hfunc2
        if(tree == k){
          pobthetaLabels[,2*i-1] = c(Var2, Var1)
          pobthetaLabels[,2*i] = c(Var1, Var2)
        } else {
          pobthetaLabels[,2*i-1] = c(Var2, as.numeric(paste0(sort(c(Var1, thetaVine$Matrix[(tree+1):k, i])), collapse='')))
          pobthetaLabels[,2*i] = c(Var1, as.numeric(paste0(sort(c(Var2, thetaVine$Matrix[(tree+1):k, i])), collapse='')))
        }
      }
      pobtheta = pobthetaNext
      
    }
  }
   
  # X_t - X_t+1 Tree
  XFit = rep(0, subset-1)
  for(i in 1:(subset-1)){
    XFit[i] = BiCopSelect(pobx[,i], pobx[,i+1], indeptest=TRUE, level=0.01)$family
  }
  XCop = Mode(XFit)
  
  # Extrapolate structure to vine over whole x parameter vector
  FullVine = list(Matrix = matrix(0, T+1+k, T+1+k),
                  family = matrix(0, T+1+k, T+1+k))
  
  diag(FullVine$Matrix)[1:(T+1)] = 0:T+1+k
  for(i in 2:(T+1)){
    for(j in 1:i-1){
      FullVine$Matrix[i, j] = T + j - i + k + 2
    }
  }
  
  FullVine$Matrix[2:(1+k)+T, 2:(1+k)+T] = thetavine$Matrix
  FullVine$Matrix[2:(1+k)+T, 1:(T+1)] = TreeEdge
  FullVine$family[2:(1+k)+T, 2:(1+k)+T] = thetavine$family
  FullVine$family[2:(1+k)+T, 1:(T+1)] = TreeCop
  FullVine$family[T+1, 1:T] = XCop
  
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
  PhiFit = FitRealMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 3], truncate=TRUE)
  GammaFit = FitRealMarginals(MCMC$theta[(MCMCreps/2 + 1):MCMCreps, 4])
  thetaDist = c(SigSqYFit$AIC, SigSqXFit$AIC, PhiFit$AIC, GammaFit$AIC)
  thetaParams = c(SigSqYFit$params, SigSqXFit$params, PhiFit$params, GammaFit$params)
  
  # Choose distribution for the last few X marginals
  subset = min(T/20, 25)
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
  } else { # CHANGE THIS
    xParams = matrix(0, T, 3)
    for(i in 1:T){
      Mixfit = normalmixEM(MCMC$x[(MCMCreps/2+1):MCMCreps,i], k = 2)
      Xparams[i, 1:2] = Mixfit$mu
      Xparams[i, 3:4] = Mixfit$sigma^2
      Xparams[i, 5] = Mixfit$lambda[1]
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
  
  # Run the VB updater to time S
  
}
