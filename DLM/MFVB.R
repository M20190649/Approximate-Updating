##Functions for fitting MFVB via gradient ascent

logjoint = function(y, params){
  theta = params[1:4]
  x = params[5:(T+5)]
  T = length(y)
  logprior = log(1/2) + dnorm(theta[2], 0, sqrt(10), log = TRUE) +
    log(densigamma(theta[3], 1, 1)) + log(densigamma(theta[4], 1, 1))
  logstates = dnorm(x[1], 0, sqrt(theta[4]), log = TRUE) + sum(dnorm(x[2:(T+1)], theta[1]*x[1:T], sqrt(theta[4])), log = TRUE)
  logy = sum(dnorm(y, x[2:(T+1)] + theta[2], sqrt(theta[3]), log = TRUE))
  return(logprior + logstates + logy)
}

logq = function(params, alpha, beta, mean, sd){
  invGamma = log(densigamma(params[1], alpha[1], beta[1])) + log(densigamma(params[2], alpha[2], beta[2]))
  normal = sum(dnorm(params[3:(T+5)], mean, sd, log = TRUE))
  return(invg + norm)
}

simul = function(alpha, beta, mean, sd){
  uniform = runif(2)
  standardNormal = rnorm(T+3)
  sigmaSq = qigamma(uniform, alpha, beta) 
  normal = mean + sd*standardNormal
  output = list(unif = uniform, z = standardNormal, sigmaSq = sigmaSq, normal = normal)
  return(output)
}

elboderiv = function(y, simulation, alpha, beta, mean, sd, argument, j){
  T = length(y)
  if(argument == "alpha"){
    alpha2 = alpha; alpha2[j] = alpha2[j] + 0.000001
    dfdl = (qigamma(simulation$uniform[j], alpha2[j], beta[j]) -  qigamma(simulation$uniform[j], alpha[j], beta[j])) / 0.000001
    dqdl = log(beta[j]) - digamma(alpha[j]) - simulation$sigmaSq[j]
    if(j == 1){
      dpdf = - (alphaY + T/2 + 1) / simulation$sigmaSq[j] + (betaY + sum((y - simulation$normal[2] - simulation$normal[4:(T+3)])^2)/2 ) / simulation$sigmaSq[j]^2
    } else {
      dpdf = - (alphaX + T/2 + 3/2) / simulation$sigmaSq[j] + (betaX + simulation$normal[3]^2/2 +
                        sum((simulation$normal[4:(T+3)] - simulation$normal[1]*simulation$normal[3:(T+2)])^2)/2) / simulation$sigmaSq[j]^2 
    }
  } else if(argument == "beta"){
    beta2 = beta; beta2[j] = beta2[j] + 0.000001
    dfdl = (qigamma(simulation$uniform[j], alpha[j], beta2[j]) -  qigamma(simulation$uniform[j], alpha[j], beta[j])) / 0.000001
    dqdl = alpha[j]/beta[j] - 1 / simulation$sigmaSq[j]
    if(j == 1){
      dpdf = - (alphaY + T/2 + 1) / simulation$sigmaSq[j] + (betaY + sum((y - simulation$normal[2] - simulation$normal[4:(T+3)])^2) ) / simulation$sigmaSq[j]^2
    } else {
      dpdf = - (alphaX + T/2 + 3/2) / simulation$sigmaSq[j] + (betaX + simulation$normal[3]^2/2 +
                        sum((simulation$normal[4:(T+3)] - simulation$normal[1]*simulation$normal[3:(T+2)])^2)/2) / simulation$sigmaSq[j]^2 
    }
  } else if(argument == "mean"){
    dfdl = 1
    dqdl = (simulation$normal[j] - mean[j]) / (sd[j]^2)
    if(j == 1) {
      dpdf = (T*simulation$normal[1] + sum(simulation$normal[4:(T+3)]) - sum(y))/simulation$sigmaSq[1] - (simulation$normal[1]-muBar)/muVar
    } else if(j == 2){
      dpdf = -(sum(simulation$normal[3:(T+2)]^2)*simulation$normal[2] - sum(simulation$normal[3:(T+2)]*simulation$normal[4:(T+3)])) / simulation$sigmaSq[2]
    } else if(j == 3){
      dpdf = (simulation$normal[2] - (1 + simulation$normal[2]^2)*simulation$normal[3]) / simulation$sigmaSq[2]
    } else {
      dpdf = - (simulation$normal[j]*(1+simulation$normal[2]^2) - simulation$normal[2]*(simulation$normal[j-1] + simulation$normal[j+1])) / simulation$sigmaSq[2]
    }
  } else {
    dfdl = simulation$z[j]
    dqdl = (simulation$normal[j] - mean[j])^2 / (sd[j]^3) - 1 / sd[j]
    if(j == 1) {
      dpdf = (T*simulation$normal[1] + sum(simulation$normal[4:(T+3)]) - sum(y))/simulation$sigmaSq[1] - (simulation$normal[1]-muBar)/muVar
    } else if(j == 2){
      dpdf = -(sum(simulation$normal[3:(T+2)]^2)*simulation$normal[2] - sum(simulation$normal[3:(T+2)]*simulation$normal[4:(T+3)])) / simulation$sigmaSq[2]
    } else if(j == 3){
      dpdf = (simulation$normal[2] - (1 + simulation$normal[2]^2)*simulation$normal[3]) / simulation$sigmaSq[2]
    } else {
      dpdf = - (simulation$normal[j]*(1+simulation$normal[2]^2) - simulation$normal[2]*(simulation$normal[j-1] + simulation$normal[j+1])) / simulation$sigmaSq[2]
    }
  
  }
  return(dpdf*dfdl - dqdl)
}

ELBO = function(y, alpha, beta, mean, sd, n = 100){
  eval = 0
  for(i in 1:n){
    sims = simul(alpha, beta, mean, sd)
    params = c(sims$sigmaSq, sims$normal)
    eval = eval + logjoint(y, params) - logq(params, alpha, beta, mean, sd)
  }
  return(eval/n)
}

MFSGA = function(y, alpha, beta, mean, sd, s, eta, threshold, maxIter){
  T = length(y)
  #use Adam optimiser
  #converge when mean change of last ten iterations is below threshold
}