PLogDens = function(y, sims){
  priordens = log(densigamma(sims[1], 1, 1)) + log(densigamma(sims[2], 1, 1)) + 
    log(1/2) + dnorm(sims[4], 0, 10, log=TRUE)
  xdens = dnorm(sims[5], 0, sims[1]/(1-sims[3]), log=TRUE)
  for(t in 6:length(sims)){
    xdens = xdens + dnorm(sims[t], sims[3]*sims[t-1], sqrt(sims[2]),log=TRUE)
  }
  ydens = 0
  for(t in 1:length(y)){
    ydens = ydens + dnorm(y[t], sims[4] + sims[5+t], sqrt(sims[1]), log=TRUE)
  }
  return(priordens + xdens + ydens)
}

QLogDens = function(unifs, sims, thetaDist, xDist, thetaParams, xParams, Vine){
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
      margins = margins + dt(sims[i], thetaParams[1, i], thetaParams[2, i], thetaParams[3, i], log=TRUE)
    }
  }
  for(i in 5:length(sims)){
    if(xDist == 2){
      margins = margins + dnorm(sims[i], xParams[1, i-4], sqrt(xParams[2, i-4]), log=TRUE)
    } else {
      margins = margins + dt(sims[i], xParams[1, i-4], xParams[2, i-4], xParams[3, i-4], log=TRUE)
    }
  }
  copulas = log(RVineDensity(unifs, Vine))
  return(margins + copulas)
}

ELBO = function(unifs, sims, y, thetaDist, xDist, thetaParams, xParams, Vine, N){
  elbo = 0
  for(i in 1:N){
    elbo = elbo + PLogDens(y, sims[i,]) -
      QLogDens(unifs[i,], sims[i,], thetaDist, xDist, thetaParams, xParams, Vine)
  }
  return(elbo/N)
}

MarginalTransform = function(unifs, thetaDist, xDist, thetaParams, xParams){
  output = rep(0, length(unifs))
  for(i in 1:2){
    if(thetaDist[i]==1){
      output[i] = qweibull(unifs[i], thetaParams[1, i], thetaParams[2, i])
    } else if(thetaDist[i]==2){
      output[i] = qgamma(unifs[i], thetaParams[1, i], thetaParams[2, i])
    } else if(thetaDist[i]==3){
      output[i] = qigamma(unifs[i], thetaParams[1, i], thetaParams[2, i])
    } else {
      output[i] = qlnorm(unifs[i], thetaParams[1, i], sqrt(thetaParams[2, i]))
    }
  }
  for(i in 3:4){
    if(thetaDist[i]==1){
      output[i] = qnorm(unifs[i], tthetaParams[1, i], sqrt(thetaParams[2, i]))
    } else {
      output[i] = qt(unifs[i], thetaParams[1, i], thetaParams[2, i], thetaParams[3, i])
    }
  }
  for(i in 5:length(unifs)){
    if(xDist == 2){
      output[i] = qnorm(unifs[i], xParams[1, i-4], sqrt(xParams[2, i-4]))
    } else {
      output[i] = qt(unifs[i], xParams[1, i-4], xParams[2, i-4], xParams[3, i-4])
    }
  }
  return(output)
}

LambdaPartial = function(unifs, sims, y, thetaDist, xDist, thetaParams, xParams, Vine, par, i, j){
  if(par == 'theta'){
    if(j <= 2){ # SigmaSq parameters on R+
      if(thetaDist[j] == 1){
        if(i == 1){
          score = score of weibull wrt param 1
        } else {
          score = score of weibull wrt param 2
        }
      } else if(thetaDist[j] == 2){
        if(i == 1){
          score = score of gamma wrt alpha
        } else {
          score = score of gamma wrt beta
        }
      } else if(thetaDist[j] == 3){
        if(i == 1){
          score = score of IG wrt alpha
        } else {
          score = score of IG wrt beta
        }
      } else {
        if(i == 1){
          score = score of LN wrt mu
        } else {
          score = score of LN wrt var
        }
      }
    } else { # Phi and Gamma on R
      if(thetaDist[j] == 1){
        if(i == 1){
          score = score of norm wrt mu
        } else {
          score = score of norm wrt var
        }
      } else {
        if(i == 1){
          score = score of t wrt loc
        } else if(i == 2){
          score = score of t wrt scale
        } else {
          score = score of t wrt df
        }
      }
    }
  } else if(par == 'x'){ 
    if(xDist == 2){
      if(i == 1){
        score = score of normal wrt mu
      } else {
        score = score of normal wrt variance
      }
    } else {
      if(i == 1){
        score = score of t wrt location
      } else if(i ==2){
        score = score of t wrt scale
      } else {
        score = score of t wrt df
      }
    }
  }
  
  ELBO = PLogDens(y, sims) - QLogDens(unifs, sims, thetaDist, xDist, thetaParams, xParams, Vine)
  return(score*ELBO)
}

EtaPartial = function(){}

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
  xParams = cbind(xParams, matrix(runif(XDist*S), ncol=S)) # Extend X parameter matrix to new latent variables
  
  MtLamT = VtLamT = matrix(0, 3, 4) # For theta Marginal Parameters
  MtLamX = VtLamX = matrix(0, Xdist, S) # For new X Marginal Parameters
  
  # Create Vine including new data
  FullVine = list(Matrix = matrix(0, T+5, T+5),
                  family = matrix(0, T+5, T+5),
                  par = matrix(0, T+5, T+5),
                  par2 = matrix(0, T+5, T+5))
  
  diag(FullVine$Matrix)[1:(T+1)] = 5:(T+5)
  for(i in 2:(T+1)){
    for(j in 1:i-1){
      FullVine$Matrix[i, j] = T + j - i + 6
    }
  }
  
  # Copy old Vine entries across
  FullVine$Matrix[2:5+T,2:5+T] = Vine$Matrix[2:5+T-S, 2:5+T-S]
  FullVine$Matrix[2:5+T,1:(T+1)] = Vine$Matrix[2:5+T-S, 1]
  FullVine$family[2:5+T,2:5+T] = Vine$family[2:5+T-S, 2:5+T-S]
  FullVine$family[2:5+T,1:(T+1)] = Vine$family[2:5+T-S, 1]
  FullVine$family[T+1, 1:T] = Vine$family[T+1, 1]
  FullVine$par[2:5+T,2:5+T] = Vine$par[2:5+T-S, 2:5+T-S]
  FullVine$par2[2:5+T,2:5+T] = Vine$par2[2:5+T-S, 2:5+T-S]
  FullVine$par[2:5+T,1:(T-S+1)] = Vine$par[2:5+T-S, 1:(T-S+1)]
  FullVine$par2[2:5+T,1:(T-S+1)] = Vine$par2[2:5+T-S, 1:(T-S+1)]
  FullVine$par[T+1, 1:(T-S)] = Vine$par[T-S+1, 1:(T-S)]
  FullVine$par2[T+1, 1:(T-S)] = Vine$par2[T-S+1, 1:(T-S)]
  
  #Initialise new parameters at last X's value
  FullVine$par[2:5+T, (T-S+2):(T+1)] = FullVine$par[2:5+T, T-S+1]
  FullVine$par2[2:5+T, (T-S+2):(T+1)] = FullVine$par2[2:5+T, T-S+1]
  FullVine$par[T+1, (T-S+1):T] = FullVine$par[T+1, T-S]
  FullVine$par2[T+1, (T-S+1):T] = FullVine$par2[T+1, T-S]
  
  MtEtaT1 = VtEtaT1 = MtEtaT2 = VtEtaT2 = matrix(0, 4, 4) # For theta copula parameters par1/par2
  MtEtaX1 = VtEtaX1 = MtEtaX2 = VtEtaX2 = matrix(0, 4, S) # for theta/X copula par1/par2
  MtEtaXX1 = VtEtaXX1 = MtEtaXX2 = VtEtaXX2 = rep(0, S) # for X_t/X_t+1 copula par1/par2
  
  # While Loop Parameters
  iter = 0
  unifs = RVineSim(FullVine, min(50, M))
  sims = MarginalTransform(VineSim)
  LB = ELBO(unifs, sims, y, thetaDist, xDist, thetaParams, xParams, Vine, 50)
  lastDiff = threshold + 1
  meanDiff = 0
  
  while(lastDiff > 5*threshold & meanDiff > threshold){
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
  
    for(m in 1:M){
      for(j in 1:4){
        for(i in 1:ifelse(i>2 & thetaDist[i]==2, 3, 2)){
          PLamT[i, j] = LambdaPartial(unifs[m,], sims[m,], y, thetaDist, xDist, thetaParams, xParams, Vine, 'theta', i, j)/M
        }
      }
    
      for(i in 1:xDist){
        for(t in 1:S){
          PLamX = LambdaPartial(unifs[m,], sims[m,], y, thetaDist, xDist, thetaParams, xParams, Vine, 'x', i, j)/M
        }
      }
    
      for(i in 1:3){
        for(j in 1:i){
          if(FullVine$par[T+2+i,T+2+j] != 0){
        }
          if(FullVine$par2[T+2+i,T+2+j] != 0){
        }
        }
      }
    
      for(i in 1:4){
        for(j in 1:S){
          if(FullVine$par[T+1+i,T-S+1+j] != 0){
          }
          if(FullVine$par[T+1+i,T-S+1+j] != 0){
          }
        }
      }
      
      for(j in 1:S){
        if(FullVine$par[T+1, T-S+j] != 0){
        }
        if(FullVine$par2[T+1, T-S+j] != 0){
        }
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
    VtLaVt = beta2*VtLaVt + (1-beta2)*PLaVt^2
    VtLamX = beta2*VtLamX + (1-beta2)*PLamX^2
    VtEtaT1 = beta2*VtEtaT1 + (1-beta2)*PEtaT1^2
    VtEtaT1 = beta2*VtEtaT2 + (1-beta2)*PEtaT2^2
    VtEtaX1 = beta2*VtEtaX1 + (1-beta2)*PEtaX1^2
    VtEtaX2 = beta2*VtEtaX2 + (1-beta2)*PEtaX2^2
    VtEtaXX1 = beta2*VtEtaXX1 + (1-beta2)*PEtaXX1^2
    VtEtaXX2 = beta2*VtEtaXX2 + (1-beta2)*PEtaXX2^2
    
    # Update Parameters
    thetaParams = thetaParams + alpha * (MtLamT/(1-beta1^iter)) / (sqrt(VtLamT/(1-beta2^iter))+e)
    xParams[,(T-S+1):T] = xParams[,(T-S+1):T] + alpha * (MtLamX/(1-beta1^iter)) / (sqrt(VtLamX/(1-beta2^iter))+e)
    FullVine$par[2:5+T,2:5+T] = FullVine$par[2:5+T,2:5+T] + alpha * (MtEtaT1/(1-beta1^iter)) / (sqrt(VtEtaT1/(1-beta2^iter))+e)
    FullVine$par2[2:5+T,2:5+T] = FullVine$par2[2:5+T,2:5+T] + alpha * (MtEtaT2/(1-beta1^iter)) / (sqrt(VtEtaT2/(1-beta2^iter))+e)
    FullVine$par[2:5+T, (T-S+2):(T+1)] = FullVine$par[2:5+T, (T-S+2):(T+1)] + alpha * (MtEtaX1/(1-beta1^iter)) / (sqrt(VtEtaX1/(1-beta2^iter))+e)
    FullVine$par2[2:5+T, (T-S+2):(T+1)] = FullVine$par2[2:5+T, (T-S+2):(T+1)] + alpha * (MtEtaX2/(1-beta1^iter)) / (sqrt(VtEtaX2/(1-beta2^iter))+e)
    FullVine$par[T+1, (T-S+1):T] = FullVine$par[T+1, (T-S+1):T] + alpha * (MtEtaXX1/(1-beta1^iter)) / (sqrt(VtEtaXX1/(1-beta2^iter))+e)
    FullVine$par2[T+1, (T-S+1):T] = FullVine$par2[T+1, (T-S+1):T] + alpha * (MtEtaXX2/(1-beta1^iter)) / (sqrt(VtEtaXX2/(1-beta2^iter))+e)
                                                  
    unifs = RVineSim(FullVine, min(50, M))
    sims = MarginalTransform(VineSim)
    LB = ELBO(unifs, sims, y, thetaDist, xDist, thetaParams, xParams, Vine, 50)
    if(iter > 10){
      lastDiff = abs(LB[iter]-LB[iter-1])
      meanDiff = mean(LB[-4:0+iter] - LB[-4:0+iter-1])
    }
  } # Close while loop
  return(List(thetaParams = thetaParams,
              xParams = xParams,
              Vine = FullVine,
              LB = LB,
              finalLB = tail(LB, 1),
              iter = iter))
}