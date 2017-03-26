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

logq = function(params, lambda){
  alpha = lambda[1:2]
  beta = lambda[3:4]
  mean = lambda[5:57]
  sd = lambda[58:110]
  invGamma = logigamma(params[1], alpha[1], beta[1]) + logigamma(params[2], alpha[2], beta[2])
  normal = sum(dnorm(params[3:(T+5)], mean, sd, log = TRUE))
  return(invGamma + normal)
}

simul = function(lambda){
  alpha = lambda[c(1, 3)]
  beta = lambda[c(2,4)]
  mean = lambda[5:57]
  sd = lambda[58:110]
  uniform = runif(2)
  standardNormal = rnorm(T+3)
  sigmaSq = qigamma(uniform[1], alpha[1], beta[1])
  sigmaSq = c(sigmaSq, qigamma(uniform[2], alpha[2], beta[2]))
  normal = mean + sd*standardNormal
  while(normal[1] < -0.95 | normal[1] > 0.95){
    standardNormal[1] = rnorm(1)
    normal[1] = mean[1] + sd[1]*standardNormal[1]
  }
  output = list(unif = uniform, z = standardNormal, sigmaSq = sigmaSq, normal = normal)
  return(output)
}

elboDeriv = function(y, simulation, lambda, j){
  T = length(y)
  h = 0.000001
  lambda2 = lambda
  lambda2[j] = lambda2[j] + h
  sigmaSq2 = simulation$sigmaSq
  normal2 = simulation$normal
  if(j <= 2){
    sigmaSq2[1] = qigamma(simulation$unif[1], lambda2[1], lambda2[2])
  } else if(j <= 4){
    sigmaSq2[2] = qigamma(simulation$unif[2], lambda2[3], lambda2[4])
  } else if(j <= 57){
    normal2[j-4] = normal2[j-4] + h
  } else {
    normal2[j-57] = normal2[j-57] + h * simulation$z[j-57]
  }
  params = c(simulation$sigmaSq, simulation$normal)
  params2 = c(sigmaSq2, normal2)
  fxh = logjoint(y, params2) - logq(params2, lambda2)
  fx = logjoint(y, params) - logq(params, lambda)
  return((fxh - fx)/h)
}

ELBO = function(y, lambda, n = 250){
  eval = 0
  for(i in 1:n){
    sims = simul(lambda)
    params = c(sims$sigmaSq, sims$normal)
    eval = eval + logjoint(y, params) - logq(params, lambda)
  }
  return(eval/n)
}

MFSGA = function(y, lambda, S, threshold, maxIter){
  #ADAM tuning parameters
  stepsize = 0.25
  beta1 = 0.9
  beta2 = 0.999
  epsilon = 10^(-8)
  #Initialise moment vectors, alpha and beta have the same length, as does mean and sd
  m = rep(0, 110)
  v = rep(0, 110)
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
    partialDeriv = rep(0, 110)
    #Take partial derivatives the mean of S monte carlo simulations
    for(s in 1:S){
      sims = simul(lambda)
      for(j in (1:110)[-c(1, 3)]){
        partialDeriv[j] = partialDeriv[j] + elboDeriv(y, sims, lambda, j)/S
      }
    }
    #Update ADAM moments
    m = beta1 * m + (1 - beta1) * partialDeriv
    v = beta2 * v + (1 - beta2) * partialDeriv^2
   #Calculate bias corrections
    mHat = m / (1 - beta1^t)
    vHat = v / (1 - beta2^t)
    #Update parameters
    lambda = lambda + stepsize * mHat / (sqrt(vHat) + epsilon)
    if(lambda[5] < -0.95){
      lambda[1] = -0.95
    } else if(lambda[1] > 0.95){
      lambda[5] = 0.95
    }
    for(j in 58:110){
      if(lambda[j] < 0.02){
        lambda[j] = 0.02
      }
    }
    #Calculate new ELBO
    LB[1:9] = LB[2:10]
    LB[10] = ELBO(y, lambda)
  } #End of while loop
  output = list(LB = LB, iter = t, sigy = lambda[1:2], sigx = lambda[3:4], mean = lambda[5:57], sd = lambda[58:110])
  return(output)
}
