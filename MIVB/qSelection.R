FitRealMarginals = function(x){
  n = length(x)
  AIC = rep(0, 2)
  # Normal
  NormMu = mean(x)
  NormVar = (n-1) / n * var(x)
  AIC[1] = 4  - 2*sum(dnorm(x, NormMu, sqrt(NormVar), log=TRUE))
  # Student T
  fitT = fitdistr(x, 't', list(m=mean(x), s=var(x), df=5))
  AIC[2] = 6 - 2*fitT$loglik
  if(AIC[1] < AIC[2]){
    return(list(Dist=1, params=c(NormMu, NormVar)))
  } else {
    return(list(Dist=2), params=fitT$estimate)
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

MarginalTransform = function(unifs, thetaDist, xDist, thetaParams, xParams){
  output = matrix(0, nrow(unifs), ncol(unifs))
  for(i in 1:2){
    if(thetaDist[i]==1){
      output[,i] = qweibull(unifs[,i], thetaParams[1, i], thetaParams[2, i])
    } else if(thetaDist[i]==2){
      output[,i] = qgamma(unifs[,i], thetaParams[1, i], thetaParams[2, i])
    } else if(thetaDist[i]==3){
      output[,i] = qigamma(unifs[,i], thetaParams[1, i], thetaParams[2, i])
    } else {
      output[,i] = qlnorm(unifs[,i], thetaParams[1, i], sqrt(thetaParams[2, i]))
    }
  }
  for(i in 3:4){
    if(thetaDist[i]==1){
      output[,i] = qnorm(unifs[,i], tthetaParams[1, i], sqrt(thetaParams[2, i]))
    } else {
      output[,i] = qt(unifs[,i], thetaParams[1, i], thetaParams[2, i], thetaParams[3, i])
    }
  }
  for(i in 5:length(unifs)){
    if(xDist == 2){
      output[,i] = qnorm(unifs[,i], xParams[1, i-4], sqrt(xParams[2, i-4]))
    } else {
      output[,i] = qt(unifs[,i], xParams[1, i-4], xParams[2, i-4], xParams[3, i-4])
    }
  }
  return(output)
}
