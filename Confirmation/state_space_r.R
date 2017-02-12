simtheta = function(Vine, lambda){
  #Steps to add lambda to Vine$par
  theta = RVineSim(Vine)
  return(theta)
}
  
ELBO = function(Vine, lambda, y, n = 1000){
 out = rep(0, n)
  for(i in 1:n){
    theta = simtheta(Vine, lambda)
    x = FFBS(y, theta)
    a = logjoint(y, x, theta)
    b = logq(theta, x, lambda)
    out[i] = a - b
  }
  return(mean(out))
}

logq = function(T, theta, lambda, x){
  phi = dnorm(theta[0], lambda[0], lambda[1], log = TRUE)
  mu = dnorm(theta[1], lambda[2], lambda[3], log = TRUE)
  sigy = log(densigamma(theta[2], lambda[4], lambda[5]))
  sigx = log(densigamma(theta[3], lambda[6], lambda[7]))
  copula = 0
  smoothed = FFBS(y, theta, FALSE)
  logqx = 0
  for(i in 1:(T+1)){
    logqx = logqx + dnorm(x[i], smoothed[i], sqrt(smoothed[T+1+i]), log = TRUE)
  }
  return(phi + mu + sigy + sigx + copula + logqx)
}

ELBOderiv = function(T, Vine, lambda, y, param){
  theta = simtheta(Vine, lambda)
  x = FFBS(y, theta)
  
  #if/else for params
  
  #return deriv(logjoint) * (logjoint - logq)
}

SGA = function(y, lambda, s, M, threshold, eta, Vine){
  k = length(initial)
  T = length(y)
  Gt = rep(0, k)
  iter = 0
  LBnew = ELBO(Vine, lambda, y)
  diff = threshold + 1
  
  while(diff > threshold | iter < 50){
    if(iter >= M){
      break
    }
    partials = rep(0, k)
    for(i in 1:s){
      for(j in 1:k){
        partials[j] = partials[j] + ELBOderiv(T, Vine, lambda, y, j)/s
      }
    }
    Gt = Gt + partials^2
    pt = eta * Gt^(-1/2)
    if(iter > 0){
      lambda = lambda + pt * partials
    }
    LBold = LBnew
    LBnew = ELBO(Vine, lambda, y)
    diff = abs(LBold - LBnew)
    iter = iter + 1
  }
  return(list(lambda = lambda, iterations = iter, ELBO = LBnew))
}