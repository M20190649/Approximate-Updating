##Functions for fitting MFVB via gradient ascent
logigamma = function(x, alpha, beta){
  alpha * log(beta) - log(gamma(alpha)) - (alpha + 1)*log(x) - beta/x
}

logjoint = function(y, params){
  theta = params[1:4]
  x = params[5:(T+5)]
  T = length(y)
  logprior = log(1/2) + dnorm(theta[4], 0, sqrt(10), log = TRUE) +
    logigamma(theta[1], 1, 1) + logigamma(theta[2], 1, 1)
  logstates = dnorm(x[1], 0, sqrt(theta[2]/(1-theta[3]^2)), log = TRUE) + sum(dnorm(x[2:(T+1)], theta[3]*x[1:T], sqrt(theta[2])), log = TRUE)
  logy = sum(dnorm(y, x[2:(T+1)] + theta[4], sqrt(theta[1]), log = TRUE))
  return(logprior + logstates + logy)
}

logq = function(params, mean, sd){
  lognorms = sum(dlnorm(params[1:2], mean[1:2], sd[1:2], log = TRUE))
  norms = sum(dnorm(params[3:55], mean[3:55], sd[3:55], log = TRUE))
  lognorms + norms
}

simul = function(mean, sd){
  standardNormal = rnorm(length(mean))
  normal = mean + sd*standardNormal
  normal[1:2] = exp(normal[1:2]) #sigma2 is lognormal
  while(normal[3] < -0.95 | normal[3] > 0.95){
    standardNormal[3] = rnorm(1)
    normal[3] = mean[3] + sd[3]*standardNormal[3]
  }
  output = list(z = standardNormal, normal = normal)
  return(output)
}

elboDeriv = function(y, simulation, mean, sd, argument, j){
  #Calculates the analytical derivative of the ELBO for a given simulated value of the parameters
  #dlog(p(f(epsilon, lambda), y))/dlambda is calculated via the chain rule of dp/df * df/dl. 
  #dloq(q(f(epsilon, lambda))) is calcualted directly as the transformed standard normal q results in a location/scale normal 
  #derivatives of that location/scale can be easily taken with respect to lambda as dq/dl
  T = length(y)
  if(j == 1){ #corresponds to sigmaSqY
    dpdf = -(T/2 + 2)/simulation$normal[j] + 1/(simulation$normal[j]^2) + 
      sum((y - simulation$normal[4] - simulation$normal[6:length(mean)])^2)/(2 * simulation$normal[j]^2)
  } else if(j == 2) { #corresponds to sigmaSqX
    dpdf = -(T/2 + 5/2)/simulation$normal[j] + 1/(simulation$normal[j]^2) +
      (1-simulation$normal[3]^2)*simulation$normal[5]^2 / (2 * simulation$normal[j]^2) +
      sum((simulation$normal[6:length(mean)] - simulation$normal[3]*simulation$normal[5:(length(mean)-1)])^2) / (2 * simulation$normal[j]^2)
  } else if(j == 3) { #phi
    dpdf = - (simulation$normal[j] * sum(simulation$normal[5:(length(mean)-1)]^2) - 
                sum(simulation$normal[5:(length(mean)-1)]*simulation$normal[6:length(mean)])) / simulation$normal[2]
  } else if(j == 4) { #mu
    dpdf = - (simulation$normal[j] - muBar) / muVar - (T*simulation$normal[j] - sum(y) + sum(simulation$normal[6:length(mean)])) / simulation$normal[1]
  } else if(j == 5) { #x0
    dpdf = - (simulation$normal[j]*(1 + simulation$normal[3]^2) - simulation$normal[3]*simulation$normal[6]) / simulation$normal[2]
  } else if(j == length(mean)) { #xT
    dpdf = - (simulation$normal[length(mean)] - simulation$normal[3]*simulation$normal[length(mean)-1]) / simulation$normal[2] -
      (simulation$normal[j] + simulation$normal[4] - y[T]) / simulation$normal[1]
  } else { #xt, t = 1, 2, ... ,T-1
    dpdf = - ((1 + simulation$normal[2]^2)*simulation$normal[j] - simulation$normal[2]*(simulation$normal[j-1] + simulation$normal[j+1])) / simulation$normal[2] -
      (simulation$normal[j] + simulation$normal[4] - y[j-5]) / simulation$normal[1]
  }
  
  if(argument == "mean"){ #lambda is a mean parameter
    if(j == 1){ #sigma2 = exp(mean + sd*z), transform handled differently
      dfdl = exp(mean[1] + sd[1]*simulation$z[1]) 
      dqdf = -simulation$z[1]
      dqdl = dqdf * dfdl
    } else if(j == 2) {
      dfdl = simulation$z[1] * exp(mean[1] + sd[1]*simulation$z[1]) 
      dqdf = -simulation$z[2]
      dqdl = dqdf * dfdl
    } else {
      dfdl = 1
      dqdl = (mean[j] - simulation$normal[j]) / sd[j]
    }
  } else { #lambda is an sd parameter
    if(j == 1){
      dfdl = exp(mean[1] + sd[1]*simulation$z[1]) 
      dqdf = -simulation$z[1]
      dqdl = dqdf * dfdl
    } else if(j == 2) {
      dfdl = simulation$z[1] * exp(mean[1] + sd[1]*simulation$z[1]) 
      dqdf = -simulation$z[2]
      dqdl = dqdf * dfdl
    } else {
      dfdl = simulation$z[j]
      dqdl = -1/sd[j] + (mean[j] - simulation$normal[j])^2 / (sd[j]^3)
    }
  }
  
  return(dpdf*dfdl - dqdl)
}

scoreDeriv = function(y, simulation, mean, sd, argument, j){
  if(argument == "mean"){
    score = (simulation$normal[j] - mean[j]) / (sd[j]^2)
  } else {
    score = - 1 / sd[j] + (simulation$normal[j] - mean[j])^2 / (sd[j]^3)
  }
  deriv = score * (logjoint(y, simulation$normal) - 
                     logq(simulation$normal, mean, sd))
  return(deriv)
}

ELBO = function(y, mean, sd, n = 250){
  eval = 0
  for(i in 1:n){
    sims = simul(mean, sd)
    eval = eval + logjoint(y, sims$normal) - logq(sims$normal, mean, sd)
  }
  return(eval/n)
}

MFSGA = function(y, mean, sd, S, threshold, maxIter, score = TRUE){
  #ADAM tuning parameters
  stepsize = 0.15
  beta1 = 0.9
  beta2 = 0.999
  epsilon = 10^(-8)
  #Initialise moment vectors, alpha and beta have the same length, as does mean and sd
  mMean = rep(0, length(mean))
  mSd = rep(0, length(mean))
  vMean = rep(0, length(mean))
  vSd = rep(0, length(mean))
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
    sdDeriv = rep(0, length(mean))
    #Take partial derivatives the mean of S monte carlo simulations
    for(s in 1:S){
      sims = simul(mean, sd)
      for(j in 1:length(mean)){
        if(score) {
          meanDeriv[j] = meanDeriv[j] + scoreDeriv(y, sims, mean, sd, "mean", j)/S
          sdDeriv[j] = sdDeriv[j] + scoreDeriv(y, sims, mean, sd, "sd", j)/S
        } else {
          meanDeriv[j] = meanDeriv[j] + elboDeriv(y, sims, mean, sd, "mean", j)/S
          sdDeriv[j] = sdDeriv[j] + elboDeriv(y, sims, mean, sd, "sd", j)/S
        }
      }
    }
    #Update ADAM moments
    mMean = beta1 * mMean + (1 - beta1) * meanDeriv
    mSd = beta1 * mSd + (1 - beta1) * sdDeriv
    vMean = beta2 * vMean + (1 - beta2) * meanDeriv^2
    vSd = beta2 * vSd + (1 - beta2) * sdDeriv^2
    #Calculate bias corrections
    mMeanHat = mMean / (1 - beta1^t)
    mSdHat = mSd / (1 - beta1^t)
    vMeanHat = vMean / (1 - beta2^t)
    vSdHat = vSd / (1 - beta2^t)
    #Update parameters
    mean = mean + stepsize * mMeanHat / (sqrt(vMeanHat) + epsilon)
    sd = sd + stepsize * mSdHat / (sqrt(vSdHat) + epsilon)
    sd[sd < 0.01] = 0.01
    if(mean[3] < -0.95){
      mean[3] = -0.95
    } else if(mean[3] > 0.95){
      mean[3] = 0.95
    }
    #Calculate new ELBO
    LB[1:9] = LB[2:10]
    LB[10] = ELBO(y, mean, sd)
  } #End of while loop
  output = list(LB = LB, iter = t, mean = mean, sd = sd)
  return(output)
}


