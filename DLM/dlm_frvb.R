##Functions for fitting a full rank Gaussian via gradient ascent
logigamma = function(x, alpha, beta){
  alpha * log(beta) - log(gamma(alpha)) - (alpha + 1)*log(x) - beta/x
}

logjoint = function(y, params){
  theta = params[1:4] #first two params are log(sigma2)
  x = params[5:(T+5)]
  T = length(y)
  logprior = log(1/2) + dnorm(theta[4], 0, sqrt(10), log = TRUE) -
    theta[1] - 1 / exp(theta[1]) - theta[2] - 1 / exp(theta[2])
  logstates = dnorm(x[1], 0, sqrt(exp(theta[2])/(1-theta[3]^2)), log = TRUE) + 
    sum(dnorm(x[2:(T+1)], theta[3]*x[1:T], sqrt(exp(theta[2])), log = TRUE))
  logy = sum(dnorm(y, x[2:(T+1)] + theta[4], sqrt(exp(theta[1])), log = TRUE))
  return(logprior + logstates + logy)
}

logq = function(x, mean, L){
  cons = -length(x)/2 * log(2*pi) + sum(log(diag(L)))
  kernel = -1/2 * t(x - mean) %*% solve(L %*% t(L)) %*% (x - mean)
  cons + kernel
}

simul = function(mean, L){
  standardNormal = rnorm(length(mean))
  normal = mean + L %*% standardNormal
  while(normal[3] > 0.95 | normal[3] < -0.95){
    standardNormal = rnorm(length(mean))
    normal = mean + L %*% standardNormal
  }
  output = list(z = standardNormal, normal = normal)
  return(output)
}

ELBO = function(y, mean, L, n = 250){
  eval = 0
  sigma = L %*% t(L)
  for(i in 1:n){
    sims = simul(mean, L)
    lj = logjoint(y, sims$normal)
    qcons = -length(x)/2 * log(2*pi) + sum(log(diag(L)))
    qkernel = -1/2 * t(x - mean) %*% solve(L %*% t(L)) %*% (x - mean)
    eval = eval + lj - (qcons + qkernel)
  }
  return(eval/n)
}

elboMeanDeriv = function(y, simulation, mean, L, i){
  #Calculates the analytical derivative of the ELBO for a given simulated value of the parameters
  #dlog(p(f(epsilon, lambda), y))/dlambda is calculated via the chain rule of dp/df * df/dl. 
  #dloq(q(f(epsilon, lambda))) is calcualted directly as the transformed standard normal q results in a location/scale normal 
  #derivatives of that location/scale can be easily taken with respect to lambda as dq/dl
  T = length(y)
  if(i == 1){ #corresponds to sigmaSqY
    dpdf = -(T/2 + 2)/exp(simulation$normal[i]) + 1/(exp(simulation$normal[i])^2) + 
      sum((y - simulation$normal[4] - simulation$normal[6:length(mean)])^2)/(2 * exp(simulation$normal[i])^2)
  } else if(i == 2) { #corresponds to sigmaSqX
    dpdf = -(T/2 + 5/2)/exp(simulation$normal[i]) + 1/(exp(simulation$normal[i])^2) +
      (1-simulation$normal[3]^2)*simulation$normal[5]^2 / (2 * exp(simulation$normal[i])^2) +
      sum((simulation$normal[6:length(mean)] - simulation$normal[3]*simulation$normal[5:(length(mean)-1)])^2) / (2 * exp(simulation$normal[i])^2)
  } else if(i == 3) { #phi
    dpdf = - (simulation$normal[i] * sum(simulation$normal[5:(length(mean)-1)]^2) - 
                sum(simulation$normal[5:(length(mean)-1)]*simulation$normal[6:length(mean)])) / exp(simulation$normal[2])
  } else if(i == 4) { #mu
    dpdf = - (simulation$normal[i] - muBar) / muVar - (T*simulation$normal[i] - sum(y) + sum(simulation$normal[6:length(mean)])) / exp(simulation$normal[1])
  } else if(i == 5) { #x0
    dpdf = - (simulation$normal[i]*(1 + simulation$normal[3]^2) - simulation$normal[3]*simulation$normal[6]) / exp(simulation$normal[2])
  } else if(i == length(mean)) { #xT
    dpdf = - (simulation$normal[length(mean)] - simulation$normal[3]*simulation$normal[length(mean)-1]) / exp(simulation$normal[2]) -
      (simulation$normal[i] + simulation$normal[4] - y[T]) / exp(simulation$normal[1])
  } else { #xt, t = 1, 2, ... ,T-1
    dpdf = - ((1 + simulation$normal[2]^2)*simulation$normal[i] - simulation$normal[2]*(simulation$normal[i-1] + simulation$normal[i+1])) / exp(simulation$normal[2]) -
      (simulation$normal[i] + simulation$normal[4] - y[i-5]) / exp(simulation$normal[1])
  }
  dfdl = 1
  if(i <= 2){
    dfdl = exp(simulation$normal[i])
  }
  dqdl = 0
  return(dpdf*dfdl - dqdl)
}

elboLDeriv = function(y, simulation, mean, L, i, j){
  #Calculates the analytical derivative of the ELBO for a given simulated value of the parameters
  #dlog(p(f(epsilon, lambda), y))/dlambda is calculated via the chain rule of dp/df * df/dl. 
  #dloq(q(f(epsilon, lambda))) is calcualted directly as the transformed standard normal q results in a location/scale normal 
  #derivatives of that location/scale can be easily taken with respect to lambda as dq/dl
  T = length(y)
  if(i == 1){ #corresponds to sigmaSqY
    dpdf = -(T/2 + 2)/exp(simulation$normal[i]) + 1/(exp(simulation$normal[i])^2) + 
      sum((y - simulation$normal[4] - simulation$normal[6:length(mean)])^2)/(2 * exp(simulation$normal[i])^2)
  } else if(i == 2) { #corresponds to sigmaSqX
    dpdf = -(T/2 + 5/2)/exp(simulation$normal[i]) + 1/(exp(simulation$normal[i])^2) +
      (1-simulation$normal[3]^2)*simulation$normal[5]^2 / (2 * exp(simulation$normal[i])^2) +
      sum((simulation$normal[6:length(mean)] - simulation$normal[3]*simulation$normal[5:(length(mean)-1)])^2) / (2 * exp(simulation$normal[i])^2)
  } else if(i == 3) { #phi
    dpdf = - (simulation$normal[i] * sum(simulation$normal[5:(length(mean)-1)]^2) - 
                sum(simulation$normal[5:(length(mean)-1)]*simulation$normal[6:length(mean)])) / exp(simulation$normal[2])
  } else if(i == 4) { #mu
    dpdf = - (simulation$normal[i] - muBar) / muVar - (T*simulation$normal[i] - sum(y) + sum(simulation$normal[6:length(mean)])) / exp(simulation$normal[1])
  } else if(i == 5) { #x0
    dpdf = - (simulation$normal[i]*(1 + simulation$normal[3]^2) - simulation$normal[3]*simulation$normal[6]) / exp(simulation$normal[2])
  } else if(i == length(mean)) { #xT
    dpdf = - (simulation$normal[length(mean)] - simulation$normal[3]*simulation$normal[length(mean)-1]) / exp(simulation$normal[2]) -
      (simulation$normal[i] + simulation$normal[4] - y[T]) / exp(simulation$normal[1])
  } else { #xt, t = 1, 2, ... ,T-1
    dpdf = - ((1 + simulation$normal[2]^2)*simulation$normal[i] - simulation$normal[2]*(simulation$normal[i-1] + simulation$normal[i+1])) / exp(simulation$normal[2]) -
      (simulation$normal[i] + simulation$normal[4] - y[i-5]) / exp(simulation$normal[1])
  }
  dfdl = simulation$z[j]
  if(i <= 2){
    dfdl = simulation$z[j]*exp(simulation$normal[i])
  }
  if(i == j){
    dqdl = -1 / L[i, i]
  } else {
    dqdl = 0
  }
  return(dpdf*dfdl - dqdl)
}

FRSGA = function(y, mean, L, S, threshold, maxIter, meanfield = FALSE, adagrad = FALSE){
  #ADAM tuning parameters
  stepsizeMean = 0.1
  stepsizeL = 0.1
  beta1 = 0.9
  beta2 = 0.999
  epsilon = 10^(-8)
  #Initialise moment vectors, alpha and beta have the same length, as does mean and sd
  mMean = rep(0, length(mean))
  mL = matrix(0, length(mean), length(mean))
  vMean = rep(0, length(mean))
  vL = matrix(0, length(mean), length(mean))
  #AdaGrad tuning parameters
  gMean = rep(0, length(mean))
  gL = matrix(0, length(mean), length(mean))
  #Step counter
  t = 0
  #Track last 10 ELBO values
  LB = rep(0, 10)
  #Break when mean change of past 10 values is lower than threshold, allow time for algorithm to get started
  while(mean(abs(LB[2:10]-LB[1:9])) > threshold | t < 50){
    t = t + 1
    if(t > maxIter){
      break
    }
    #Initialise derivatives to 0
    meanDeriv = rep(0, length(mean))
    LDeriv = matrix(0, length(mean), length(mean))
    #Take partial derivatives the mean of S monte carlo simulations
    for(s in 1:S){
      sims = simul(mean, L)
      for(i in 1:length(mean)){
        meanDeriv[i] = meanDeriv[i] + elboMeanDeriv(y, sims, mean, L, i)/S
        if(meanfield){
          LDeriv[i,i] = LDeriv[i,i] + elboLDeriv(y, sims, mean, L, i, i)/S
        } else {
          for(j in 1:i){
            LDeriv[i,j] = LDeriv[j] + elboLDeriv(y, sims, mean, L, i, j)/S
          }
        }
      }
    }
    if(adagrad){
      gMean = gMean + meanDeriv^2
      gL = gL + LDeriv^2
      pMean = stepsizeMean * gMean^(-1/2)
      pL = stepsizeL * gL^(-1/2)
      if(t >= 2){ #first step always has deriv and p cancel
        mean = mean + pMean * meanDeriv
        L = L + pL * LDeriv 
      }
    } else {
      #Update ADAM moments
      mMean = beta1 * mMean + (1 - beta1) * meanDeriv
      mL = beta1 * mL + (1 - beta1) * LDeriv
      vMean = beta2 * vMean + (1 - beta2) * meanDeriv^2
      vL = beta2 * vL + (1 - beta2) * LDeriv^2
      #Calculate bias corrections
      mMeanHat = mMean / (1 - beta1^t)
      mLHat = mL / (1 - beta1^t)
      vMeanHat = vMean / (1 - beta2^t)
      vLHat = vL / (1 - beta2^t)
      #Update parameters, first update always has m and v cancel
      if(t >= 2){
        mean = mean + stepsizeMean * mMeanHat / (sqrt(vMeanHat) + epsilon)
        L = L + stepsizeL * mLHat / (sqrt(vLHat) + epsilon)
      }
    }
    #Calculate new ELBO
    LB[1:9] = LB[2:10]
    LB[10] = ELBO(y, mean, L)
  } #End of while loop
  if(meanfield){
    output = list(LB = LB, iter = t, mean = mean, sd = diag(L)^2)
  } else {
    output = list(LB = LB, iter = t, mean = mean, L = L)
  }
  return(output)
}