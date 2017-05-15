FitRealMarginals = function(x){
  n = length(x)
  AIC = rep(0, 2)
  # Normal
  NormMu = mean(x)
  NormVar = (n-1) / n * var(x)
  AIC[1] = 4  - 2*sum(dnorm(x, NormMu, sqrt(NormVar), log=TRUE))
  # Student T
  fitT = fitdistr(x, 't')
  AIC[2] = 6 - 2*fitT$loglik
  if(AIC[1] < AIC[2]){
    return(list(Dist=1, params=c(NormMu, NormVar), AIC=AIC))
  } else {
    return(list(Dist=2, params=fitT$estimate, AIC=AIC))
  }
}

FitPositiveMarginals = function(x){
  require(MASS)
  require(pscl)
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
    return(list(Dist=1, params=fitW$estimate, AIC=AIC))
  } else if(distChoice==2){
    return(list(Dist=2, params=fitG$estimate, AIC=AIC))
  } else if(distChoice==3){
    return(list(Dist=3, params=c(a, b), AIC=AIC))
  } else {
    return(list(Dist=4, params=fitLN$estimate, AIC=AIC))
  }
}

Mode = function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

BiCopAIC = function(x, y, family){
  output = rep(0, length(family))
  tau = cor(x, y, method='kendall')
  for(i in seq_along(family)){
    if(family[i] %in% c(23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40,
                        124, 134, 224, 234) & tau > 0){
      output[i] = 99999
    } else if(family[i] %in% c(3, 4, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20,
                               104, 114, 204, 214) & tau < 0){
      output[i] = 99999
    } else
    output[i] = BiCopEst(x, y, family[i])$AIC
  }
  return(output)
}


RestrictedVineSelect = function(T, subset, k, MCMCreps, MCMC){
  require(VineCopula)
  # Create Vine Matrix for X subset
  Vine = list(Matrix = matrix(0, subset+k, subset+k), 
              Family = matrix(0, subset+k, subset+k),
              par = matrix(0, subset+k, subset+k),
              par2 = matrix(0, subset+k, subset+k))
  diag(Vine$Matrix)[1:subset] = subset:1+k
  for(i in 2:(subset)){
    for(j in 1:(i-1)){
      Vine$Matrix[i, j] = i - j + k 
    }
  }
  
  # Fill in Theta Subtree, fit to a subsample of MCMC draws
  v = sample((MCMCreps/2+1):MCMCreps, min(500, MCMCreps/2))
  pobtheta = array(0, dim = c(length(v), 2*(k-1), k))
  pobtheta[,1:k,1] = pobs(MCMC$theta[v,])
  pobthetaLabels = array(0, dim=c(2, 2*(k-1), k))
  pobx = pobs(MCMC$x[v,-subset:0+T+1])
  thetaVine = RVineStructureSelect(pobtheta[,1:k,1])
  Vine$Matrix[1:k+subset,1:k+subset] = thetaVine$Matrix
  Vine$Family[1:k+subset,1:k+subset] = thetaVine$family
  
  
  # Misc Variables to initialise
  TreeEdge = rep(0, k)
  family = c(0, 1, 2, 2, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 
             27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234)
  TreeAIC = matrix(0, length(family), subset)
  TreeCop = rep(0, k)
  candidates = 1:k
  whichTheta = rep(0, k)
  # Sequential selection of edge, selection of copula family and calculation of data for next tree
  for(tree in 1:k){
    # Select next edge
    if(tree == 1){
      treeTau = TauMatrix(cbind(pobtheta[,1:k, tree], pobx))
      TreeEdge[tree] = candidates[which.max(rowSums(abs(treeTau[1:k, 1:subset+k])))]
    } else {
      sameCond = which(pobthetaLabels[2,,tree] == as.numeric(paste0(sort(TreeEdge[1:(tree-1)]), collapse='')))
      candidates = pobthetaLabels[1, sameCond, tree]
      if(length(candidates)==1){
        TreeEdge[tree] = candidates
      } else {
        treeTau = TauMatrix(cbind(pobtheta[,sameCond, tree], pobx))
        TreeEdge[tree] = candidates[which.max(rowSums(abs(treeTau[1:length(candidates), 1:subset+length(candidates)])))]
      }
    }
    
    # Select copula family for that edge
    if(tree == 1){
      whichTheta[tree] = TreeEdge[tree]
    } else {
      whichTheta[tree] = which(apply(pobthetaLabels[,,tree], 2, function(x) x[1] == TreeEdge[tree] &
                                 x[2] == as.numeric(paste0(sort(TreeEdge[1:(tree-1)]), collapse=''))))
    }
    for(i in 1:subset){ 
      TreeAIC[,i] = BiCopAIC(pobx[,i], pobtheta[,whichTheta[tree], tree], family)
    }
    TreeCop[tree] = family[which.min(rowSums(TreeAIC))]
    
    # Calculate next level of xt | theta
    pobxNext = matrix(0, length(v), subset) 
    for(i in 1:subset){
      CopEstim = BiCopEst(pobx[,i], pobtheta[,whichTheta[tree], tree], family=TreeCop[tree])
      pobxNext[,i] = BiCopHfunc2(pobx[,i], pobtheta[,whichTheta[tree], tree], family=TreeCop[tree], par=CopEstim$par, par2=CopEstim$par2)
      Vine$par[subset+k+1-tree, subset+1-i] = CopEstim$par
      Vine$par2[subset+k+1-tree, subset+1-i] = CopEstim$par2
    }
    pobx = pobxNext
    
    # Calculate next level of theta_i | theta_(!i)
    if(tree < k){
      for(i in 1:(k-tree)){
        Var1 = diag(thetaVine$Matrix)[i]
        Var2 = thetaVine$Matrix[k-tree+1, i]
        Hfunc = BiCopHfunc(pobtheta[,Var1, tree], pobtheta[,Var2, tree], family=thetaVine$family[k-tree+1, i],
                           par=thetaVine$par[k-tree+1, i], par2=thetaVine$par2[k-tree+1, i])
        pobtheta[,2*i-1, tree+1] = Hfunc$hfunc1
        pobtheta[,2*i, tree+1] = Hfunc$hfunc2
        if(tree == 1){
          pobthetaLabels[,2*i-1, tree+1] = c(Var2, Var1)
          pobthetaLabels[,2*i, tree+1] = c(Var1, Var2)
        } else {
          pobthetaLabels[,2*i-1, tree+1] = c(Var2, as.numeric(paste0(sort(c(Var1, thetaVine$Matrix[(k-tree+2):k, i])), collapse='')))
          pobthetaLabels[,2*i, tree+1] = c(Var1, as.numeric(paste0(sort(c(Var2, thetaVine$Matrix[(k-tree+2):k, i])), collapse='')))
        }
      }
    }
  }
  
  # X_t - X_t+1 Tree
   XAIC = matrix(0, length(family), subset-1)
  for(i in 1:(subset-1)){
    XAIC[,i] = BiCopAIC(pobx[,i], pobx[,i+1])
  }
  XCop = family[which.min(rowSums(XAIC))]
  for(i in 1:(subset-1)){
    CopEstim = BiCopEst(pobx[,i], pobx[,i+1], family=XCop)
    Vine$par[subset, subset-i] = CopEstim$par
    Vine$par2[subset, subset-i] = CopEstim$par2
  }
  
  # Extrapolate structure to vine over whole x parameter vector
  FullVine = list(Matrix = matrix(0, T+1+k, T+1+k),
                  family = matrix(0, T+1+k, T+1+k),
                  par = matrix(0, T+1+k, T+1+k), 
                  par2 = matrix(0,T+1+k, T+1+k))
  
  diag(FullVine$Matrix)[1:(T+1)] = T:0+1+k
  for(i in 2:(T+1)){
    for(j in 1:i-1){
      FullVine$Matrix[i, j] = i - j + k
    }
  }
  
  FullVine$Matrix[2:(1+k)+T, 2:(1+k)+T] = thetaVine$Matrix
  FullVine$Matrix[2:(1+k)+T, 1:(T+1)] = TreeEdge
  FullVine$family[2:(1+k)+T, 2:(1+k)+T] = thetaVine$family
  FullVine$family[2:(1+k)+T, 1:(T+1)] = TreeCop
  FullVine$family[T+1, 1:T] = XCop
  FullVine$par[(2:(1+k)+T), 2:(1+k)+T] = thetaVine$par
  FullVine$par2[(2:(1+k)+T), 2:(1+k)+T] = thetaVine$par2
  FullVine$par[0:k+T+1, 1:subset] = Vine$par[0:k + subset, 1:subset]
  
  
  # Estimate remaining par values
  pobxNew = pobs(MCMC$x[v, 1:(T-subset+1)])
  for(tree in 1:k){
    pobxNext = matrix(0, length(v), T-subset+1)
    for(i in 1:(T-subset+1)){
      CopEstim = BiCopEst(pobxNew[,i], pobtheta[,whichTheta[tree], tree], family=TreeCop[tree])
      pobxNext[,i] = BiCopHfunc2(pobxNew[,i], pobtheta[,whichTheta[tree], tree], family=TreeCop[tree], par=CopEstim$par, par2=CopEstim$par2)
      FullVine$par[T+2+k-tree, T+2-i] = CopEstim$par
      FullVine$par2[T+2+k-tree, T+2-i] = CopEstim$par2
    }
    pobxNew = pobxNext
  }
  for(i in 1:(T-subset)){
    CopEstim = BiCopEst(pobxNew[,i], pobxNew[,i+1], family=XCop)
    FullVine$par[T+1, T+1-i] = CopEstim$par
    FullVine$par2[T+1, T+1-i] = CopEstim$par2
  }
  # Final Copula between first subset X and last non-subset X
  CopEstim = BiCopEst(pobxNew[,T-subset], pobx[,1], family=XCop)
  FullVine$par[T+1, subset] = CopEstim$par
  FullVine$par2[T+1, subset] = CopEstim$par2
  
  return(FullVine)
}

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
    y[t] = gamma + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
  }
  return(y)
}
