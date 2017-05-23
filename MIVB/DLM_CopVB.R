dt_ls <- function(x, df, mu, a) {1/a * dt((x - mu)/a, df)}
pt_ls <- function(x, df, mu, a) {pt((x - mu)/a, df)}
qt_ls <- function(prob, df, mu, a) {qt(prob, df)*a + mu}
rt_ls <- function(n, df, mu, a) {rt(n,df)*a + mu}

MarginalTransform = function(unifs, thetaDist, xDist, thetaParams, xParams, T, j=0){
  if(j > 0){
    if(j <= 2){
      if(thetaDist[j]==1){
        output = qweibull(unifs, thetaParams[1, j], thetaParams[2, j])
      } else if(thetaDist[j]==2){
        output = qgamma(unifs, thetaParams[1, j], thetaParams[2, j])
      } else if(thetaDist[j]==3){
        output = qigamma(unifs, thetaParams[1, j], thetaParams[2, j])
      } else {
        output = qlnorm(unifs, thetaParams[1, j], sqrt(thetaParams[2, j]))
      }
    } else if(j <= 4){
      if(thetaDist[j]==1){
        output = qnorm(unifs, thetaParams[1, j], sqrt(thetaParams[2, j]))
      } else {
        output = qt_ls(unifs, thetaParams[3, j], thetaParams[1, j], thetaParams[2, j])
      }
    } else {
      if(xDist == 2){
        output = qnorm(unifs, xParams[1, T+6-j], sqrt(xParams[2, T+6-j]))
      } else {
        output = qt_ls(unifs, xParams[3,T+6-j], xParams[1, T+6-j], xParams[2, T+6-j])
      }
    }
  } else {
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
        output[,i] = qnorm(unifs[,i], thetaParams[1, i], sqrt(thetaParams[2, i]))
      } else {
        output[,i] = qt_ls(unifs[,i], thetaParams[3, i], thetaParams[1, i], thetaParams[2, i])
      }
    }
    for(i in 5:ncol(unifs)){
      if(xDist == 2){
        output[,i] = qnorm(unifs[,i], xParams[1, T+6-i], sqrt(xParams[2, T+6-i]))
      } else {
        output[,i] = qt_ls(unifs[,i], xParams[3, T+6-i], xParams[1, T+6-i], xParams[2, T+6-i])
      }
    }
  }
  return(output)
}

PLogDens = function(y, sims){
  priordens = log(densigamma(sims[1], 1, 1)) + log(densigamma(sims[2], 1, 1)) + 
    log(1/2) + dnorm(sims[4], 0, 10, log=TRUE)
  xdens = dnorm(sims[5], 0, sims[1]/(1-0.95^2), log=TRUE)
  for(t in 6:length(sims)){
    xdens = xdens + dnorm(sims[t], sims[3]*sims[t-1], sqrt(sims[2]),log=TRUE)
  }
  ydens = 0
  for(t in 1:length(y)){
    ydens = ydens + dnorm(y[t], sims[4] + sims[5+t], sqrt(sims[1]), log=TRUE)
  }
  return(priordens + xdens + ydens)
}

QLogDens = function(unifsdep, sims, thetaDist, xDist, thetaParams, xParams, Vine){
  margins = 0
  for(i in 1:2){
    if(thetaDist[i]==1){
      margins = margins + dweibull(sims[i], thetaParams[1, i], thetaParams[2, i], log=TRUE)
    } else if(thetaDist[i]==2){
      margins = margins + dgamma(sims[i], thetaParams[1, i], thetaParams[2, i], log=TRUE)
    } else if(thetaDist[i]==3){
      margins = margins + log(densigamma(sims[i], thetaParams[1, i], thetaParams[2, i]))
    } else {
      margins = margins + dlnorm(sims[i], thetaParams[1, i], sqrt(thetaParams[2, i]), log=TRUE)
    }
  }
  for(i in 3:4){
    if(thetaDist[i]==1){
      margins = margins + dnorm(sims[i], thetaParams[1, i], sqrt(thetaParams[2, i]), log=TRUE)
    } else {
      margins = margins + log(dt_ls(sims[i], thetaParams[3, i], thetaParams[1, i], thetaParams[2, i]))
    }
  }
  for(i in 5:length(sims)){
    if(xDist == 2){
      margins = margins + dnorm(sims[i], xParams[1, i-4], sqrt(xParams[2, i-4]), log=TRUE)
    } else {
      margins = margins + log(dt_ls(sims[i], xParams[3, i-4], xParams[1, i-4], xParams[2, i-4]))
    }
  }
  copulas = 0
  for(i in 1:5+T){
    for(j in 1:(T+5)){
      if(Vine$family[i, j] != 0){
        u1 = diag(Vine$Matrix)[j]
        u2 = Vine$Matrix[i, j]
        copulas = copulas + log(BiCopPDF(unifsdep[u1], unifsdep[u2], Vine$family[i, j], Vine$par[i, j], Vine$par2[i, j]))
      }
    }
  }
  return(margins + copulas)
}

VineSim = function(unifs, Vine, T){
  n = ncol(unifs)
  nr = nrow(unifs)
  m = Vine$Matrix
  mTheta = diag(Vine$Matrix)[2:5+T]
  for(i in 1:4){
    for(j in 1:4){
      m[T+1+i, T+1+j] = ifelse(m[T+1+i, T+1+j] == 0, 0, 
                               switch(which(mTheta == m[T+1+i, T+1+j]), 4, 3, 2, 1))
    }
  }
  mmax = matrix(0, n, n)
  for(i in 1:n){
    for(j in 1:n){
      mmax[i, j] = max(m[i:n, j])
    }
  }
  vd = array(0, dim = c(n, n, nr))
  vind = array(0, dim = c(n, n, nr))
  vd[n,,] = unifs
  z1 = array(0, dim = c(n, n, nr))
  z2 = array(0, dim = c(n, n, nr))
  x = unifs[,1]
  for(k in (n-1):1){
    for(i in max(k+1, n-4):n){
      if(mmax[i, k] == m[i, k]){
        z2[i, k,] = vd[i, n-mmax[i,k]+1,]
      } else {
        z2[i, k,] = vind[i, n-mmax[i,k]+1,]
      }
      vd[n, k,] = BiCopHinv1(vd[n,k,], z2[i,k,], Vine$family[i, k], Vine$par[i, k], Vine$par2[i, k])
    }
    x = cbind(x, vd[n, k,])
    for(i in n:max(k+1, n-5)){
      z1[i, k,] = vd[i, k,]
      vd[i-1, k,] = BiCopHfunc1(z1[i, k,], z2[i, k,], Vine$family[i, k], Vine$par[i, k], Vine$par2[i, k])
      vind[i-1, k,] = BiCopHfunc1(z1[i, k,], z2[i, k,], Vine$family[i, k], Vine$par[i, k], Vine$par2[i, k])
    }
  }
  return(cbind(x[,rev(mTheta)], x[,5:n]))
}

ELBO = function(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, Vine, N){
  elbo = 0
  for(i in 1:N){
    elbo = elbo + PLogDens(y, sims[i,]) -
      QLogDens(unifsdep[i,], sims[i,], thetaDist, xDist, thetaParams, xParams, Vine)
  }
  return(elbo/N)
}

Partial = function(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams, xParams, 
                   Vine, M, what, i, j, par=NULL){
  h = 0.0000001
  if(what == 'LamT'){
    thetaParams2 = thetaParams
    thetaParams2[i, j] = thetaParams2[i, j] + h
    newsim = sims
    newsim[,j] = MarginalTransform(unifsdep[,j], thetaDist, xDist, thetaParams2, xParams, T, j)
    deriv = 0
    for(i in 1:M){
      deriv = deriv + (PLogDens(y, newsim[i,]) - PLogDens(y, sims[i,]) -
        QLogDens(unifsdep[i,], newsim[i,], thetaDist, xDist, thetaParams2, xParams, Vine) +
        QLogDens(unifsdep[i,], sims[i,], thetaDist, xDist, thetaParams, xParams, Vine))/h
    }
  } else if(what == 'LamX'){
    xParams2 = xParams
    xParams2[i, j] = xParams2[i, j] + h
    newsim = sims
    newsim[,j+5+T-S] = MarginalTransform(unifsdep[,j+5+T-S], thetaDist, xDist, thetaParams, xParams2, T, j+5+T-S)
    deriv = 0
    for(i in 1:M){
      deriv = deriv + (PLogDens(y, newsim[i,]) - PLogDens(y, sims[i,]) -
                           QLogDens(unifsdep[i,], newsim[i,], thetaDist, xDist, thetaParams, xParams2, Vine) +
                           QLogDens(unifsdep[i,], sims[i,], thetaDist, xDist, thetaParams, xParams, Vine))/h
    }
  } else {
    Vine2 = Vine
    if(what == 'EtaT') {
      if(par==1){
        Vine2$par[T+2+i,T+1+j] = Vine2$par[T+2+i,T+1+j] + h
      } else {
        Vine2$par2[T+2+i,T+1+j] = Vine2$par2[T+2+i,T+1+j] + h
      }
    } else if(what == 'EtaX') {
      if(par==1){
        Vine2$par[T+1+i,j] = Vine2$par[T+1+i,j] + h
      } else {
        Vine2$par2[T+1+i,j] = Vine2$par2[T+1+i,j] + h
      }
    } else if(what == 'EtaXX') {
      if(par==1){
        Vine2$par[T+1+i,j] = Vine2$par[T+1, j] + h
      } else {
        Vine2$par2[T+1+i,j] = Vine2$par2[T+1, j] + h
      }
    }
    newdep = VineSim(unifs[1:M,], Vine2, T)
    newsim = MarginalTransform(newdep, thetaDist, xDist, thetaParams, xParams, T)
    deriv = 0
    for(i in 1:M){
      deriv = deriv + (PLogDens(y, newsim[i,]) - PLogDens(y, sims[i,]) -
                           QLogDens(newdep[i,], newsim[i,], thetaDist, xDist, thetaParams, xParams, Vine2) +
                           QLogDens(unifsdep[i,], sims[i,], thetaDist, xDist, thetaParams, xParams, Vine))/h
    }
  }
  return(deriv/M)
}



CopulaVB = function(y, S, thetaDist, xDist, thetaParams, xParams, Vine, 
                    M, maxIter, threshold=0.01, alpha=0.01, beta1=0.9, beta2=0.999){
  require(VineCopula)
  require(pscl)
  T = length(y)
  if(S > T){
    print('Error: S must be less than or equal to the length of y')
    return(NULL)
  }
  e = 1E-8 # Adam optimsier parameter
  xDist = ifelse(xDist==1, 2, 3) # Convert from distribution marker to number of parameters in distribution
  xParams = cbind(matrix(runif(xDist*S), nrow=xDist), xParams) # Extend X parameter matrix to new latent variables
  
  MtLamT = VtLamT = matrix(0, 3, 4) # For theta Marginal Parameters
  MtLamX = VtLamX = matrix(0, xDist, S) # For new X Marginal Parameters
  MtEtaT1 = VtEtaT1 = MtEtaT2 = VtEtaT2 = matrix(0, 4, 4) # For theta copula parameters par1/par2
  MtEtaX1 = VtEtaX1 = MtEtaX2 = VtEtaX2 = matrix(0, 4, S) # for theta/X copula par1/par2
  MtEtaXX1 = VtEtaXX1 = MtEtaXX2 = VtEtaXX2 = rep(0, S) # for X_t/X_t+1 copula par1/par2
  
  # Create Vine including new data
  FullVine = list(Matrix = matrix(0, T+5, T+5),
                  family = matrix(0, T+5, T+5),
                  par = matrix(0, T+5, T+5),
                  par2 = matrix(0, T+5, T+5))
  
  diag(FullVine$Matrix)[1:(T+1)] = 5:(T+5)
  for(i in 2:(T+1)){
    for(j in 1:i-1){
      FullVine$Matrix[i, j] = i - j + 4
    }
  }
  
  # Copy old Vine entries across
  FullVine$Matrix[(S+1):(T+5),(S+1):(T+5)] = Vine$Matrix
  FullVine$family[(S+1):(T+5), (S+1):(T+5)] = Vine$family
  FullVine$par[(S+1):(T+5), (S+1):(T+5)] = Vine$par
  FullVine$par2[(S+1):(T+5), (S+1):(T+5)] = Vine$par2
  # Extend to new X variables
  FullVine$Matrix[(T+2):(T+5), 1:S] = FullVine$Matrix[(T+2):(T+5), S+1]
  FullVine$family[(T+1):(T+5), 1:S] = FullVine$family[(T+1):(T+5), S+1]
  FullVine$par[(T+1):(T+5), 1:S] = FullVine$par[(T+1):(T+5), S+1]
  FullVine$par2[(T+1):(T+5), 1:S] = FullVine$par2[(T+1):(T+5), S+1]
  

  
  # While Loop Parameters
  iter = 0
  unifs = matrix(runif(max(50, M)*(T+5)), ncol=T+5)
  unifsdep = VineSim(unifs, FullVine, T) 
  sims = MarginalTransform(unifsdep, thetaDist, xDist, thetaParams, xParams, T)
  LB = ELBO(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, FullVine, 50)
  lastDiff = threshold + 1
  meanDiff = 0
  
  print('Start Loop')
  while((lastDiff > 5*threshold) | (meanDiff > threshold) | (iter < 10)){
    iter = iter + 1
    if(iter > maxIter){
      break
    }
    # Reset Partial Derivatives
    PLamT = matrix(0, 3, 4)
    PLamX = matrix(0, xDist, S)
    PEtaT1 = PEtaT2 = matrix(0, 4, 4)
    PEtaX1 = PEtaX2 = matrix(0, 4, S)
    PEtaXX1 = PEtaXX2 = rep(0, S)
    
    # Calculate Partial Derivatives
    # j = variable, i = parameter

    for(j in 1:4){
      for(i in 1:ifelse(j>2, thetaDist[j]+1, 2)){
        PLamT[i, j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                              xParams, FullVine, M, 'LamT', i, j)
      }
    }
    for(i in 1:xDist){
      for(j in 1:S){
        PLamX[i, j] =  Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                               xParams, FullVine, M, 'LamX', i, j)
      }
    }
    for(i in 1:3){
      for(j in 1:i){
        if(FullVine$par[T+2+i,T+2+j] != 0){
          PEtaT1[i, j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                                 xParams, FullVine, M, 'EtaT', i, j, 1)
        }
        if(FullVine$par2[T+2+i,T+2+j] != 0){
          PEtaT2[i, j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                                 xParams, FullVine, M, 'EtaT', i, j, 2)
        }
      }
    }
    for(i in 1:4){
      for(j in 1:S){
        if(FullVine$par[T+1+i,j] != 0){
          PEtaX1[i, j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                                 xParams, FullVine, M, 'EtaX', i, j, 1)
        }
        if(FullVine$par[T+1+i,j] != 0){
          PEtaX2[i, j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                                 xParams, FullVine, M, 'EtaX', i, j, 2)
        }
      }
    }
    for(j in 1:S){
      if(FullVine$par[T+1, j] != 0){
        PEtaXX1[j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                               xParams, FullVine, M, 'EtaXX', 0, j, 1)
      }
      if(FullVine$par2[T+1, j] != 0){
        PEtaXX1[j] = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                                xParams, FullVine, M, 'EtaXX', 0, j, 2)
      }
    }
    
    # Update lambda and eta parameters with Adam
    # Update first moment estimators
    MtLamT = beta1*MtLamT + (1-beta1)*PLamT
    MtLamX = beta1*MtLamX + (1-beta1)*PLamX
    MtEtaT1 = beta1*MtEtaT1 + (1-beta1)*PEtaT1
    MtEtaT1 = beta1*MtEtaT2 + (1-beta1)*PEtaT2
    MtEtaX1 = beta1*MtEtaX1 + (1-beta1)*PEtaX1
    MtEtaX2 = beta1*MtEtaX2 + (1-beta1)*PEtaX2
    MtEtaXX1 = beta1*MtEtaXX1 + (1-beta1)*PEtaXX1
    MtEtaXX2 = beta1*MtEtaXX2 + (1-beta1)*PEtaXX2
    
    # Update second moment estimators
    VtLamT = beta2*VtLamT + (1-beta2)*PLamT^2
    VtLamX = beta2*VtLamX + (1-beta2)*PLamX^2
    VtEtaT1 = beta2*VtEtaT1 + (1-beta2)*PEtaT1^2
    VtEtaT1 = beta2*VtEtaT2 + (1-beta2)*PEtaT2^2
    VtEtaX1 = beta2*VtEtaX1 + (1-beta2)*PEtaX1^2
    VtEtaX2 = beta2*VtEtaX2 + (1-beta2)*PEtaX2^2
    VtEtaXX1 = beta2*VtEtaXX1 + (1-beta2)*PEtaXX1^2
    VtEtaXX2 = beta2*VtEtaXX2 + (1-beta2)*PEtaXX2^2
    
    # Correct Bias and Update Parameters
    thetaParams = thetaParams + alpha * (MtLamT/(1-beta1^iter)) / (sqrt(VtLamT/(1-beta2^iter))+e)
    xParams[,(T-S+1):T] = xParams[,(T-S+1):T] + alpha * (MtLamX/(1-beta1^iter)) / (sqrt(VtLamX/(1-beta2^iter))+e)
    FullVine$par[2:5+T,2:5+T] = FullVine$par[2:5+T,2:5+T] + alpha * (MtEtaT1/(1-beta1^iter)) / (sqrt(VtEtaT1/(1-beta2^iter))+e)
    FullVine$par2[2:5+T,2:5+T] = FullVine$par2[2:5+T,2:5+T] + alpha * (MtEtaT2/(1-beta1^iter)) / (sqrt(VtEtaT2/(1-beta2^iter))+e)
    FullVine$par[2:5+T, 1:S] = FullVine$par[2:5+T, 1:S] + alpha * (MtEtaX1/(1-beta1^iter)) / (sqrt(VtEtaX1/(1-beta2^iter))+e)
    FullVine$par2[2:5+T, 1:S] = FullVine$par2[2:5+T, 1:S] + alpha * (MtEtaX2/(1-beta1^iter)) / (sqrt(VtEtaX2/(1-beta2^iter))+e)
    FullVine$par[T+1, 1:S] = FullVine$par[T+1, 1:S] + alpha * (MtEtaXX1/(1-beta1^iter)) / (sqrt(VtEtaXX1/(1-beta2^iter))+e)
    FullVine$par2[T+1, 1:S] = FullVine$par2[T+1, 1:S] + alpha * (MtEtaXX2/(1-beta1^iter)) / (sqrt(VtEtaXX2/(1-beta2^iter))+e)
                                                  
    unifs = matrix(runif(max(50, M)*(T+5)), ncol=T+5)
    unifsdep = VineSim(unifs, FullVine, T) 
    sims = MarginalTransform(unifsdep, thetaDist, xDist, thetaParams, xParams, T)
    LB = c(LB, ELBO(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, FullVine, 50))
    if(iter > 10){
      lastDiff = abs(LB[iter+1]-LB[iter])
      meanDiff = mean(LB[1:5 + iter-4] - LB[1:5 + iter-5])
    }
    if(iter %% 10 == 0){
      print(paste0('Iteration: ', iter))
      print(paste0('ELBO: ', LB[iter+1]))
    }
  } # Close while loop
  return(list(thetaParams = thetaParams,
              xParams = xParams,
              Vine = FullVine,
              LB = LB,
              finalLB = tail(LB, 1),
              iter = iter))
}
