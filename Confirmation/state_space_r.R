simtheta = function(Vine, lambda){
  Vine$par[4, 3] = lambda[10]
  Vine$par[4, 2] = lambda[11]
  Vine$par2[4, 2] = lambda[12]
  theta = RVineSim(1, Vine)
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

dent = function(x, df, mean, var, log = FALSE){
  front = gamma((df + 1)/ 2) / (gamma(1/2) * gamma(df/2))
  middle = (df * var)^(-1/2)
  end = (1 + 1/df * ((x - mean)/sqrt(var))^2)^(-(v+1)/2)
  if(log){
    return(log(front) + log(middle) + log(end))
  } else {
    return(front*middle*end)
  }
}

logq = function(T, theta, lambda, x){
  phi = log(dtruncnorm(theta[1], -1, 1, lambda[1], lambda[2]))
  mu = dent(theta[2], lambda[3], lambda[4], lambda[5], log = TRUE)
  sigy = log(densigamma(theta[3], lambda[6], lambda[7]))
  sigx = log(densigamma(theta[4], lambda[8], lambda[9]))
  copula1 = log(BiCopPDF(pigamma(theta[3], lambda[6], lambda[7]), 
                         pigamma(theta[4], lambda[8], lambda[9]), family = 1, par = lambda[10]))
  copula2 = log(BiCopPDF(ptruncnorm(theta[1], lambda[1], lambda[2]), 
                         pigamma(theta[4], lambda[8], lambda[9]), 
                         family = 30, par = lambda[11], par2 = lambda[12]))
  smoothed = FFBS(y, theta, FALSE)
  logqx = rep(0, T+1)
  for(i in 1:(T+1)){
    logqx[i]= dnorm(x[i], smoothed[i], sqrt(smoothed[T+1+i]), log = TRUE)
  }
  return(phi + mu + sigy + sigx + copula1 + copula2 + sum(logqx))
}

ELBOderiv = function(T, Vine, lambda, y, param){
  theta = simtheta(Vine, lambda)
  x = FFBS(y, theta)
  
  #if/else for params
  
  #return deriv(logq) * (logjoint - logq)
}

#lambda = 1. Phi Mean, 2. Phi Variance, 3. Mu df, 4. Mu mean, 5. Mu Variance, 6. Sigy Alpha, 7. Sigy Beta, 8. Sigx Alpha, 9. Sigx Beta
#10 Sigx/Sigy Rho, #11 Sigx/Phi p1, #12 Sigx/Phi p2
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
  return(list(lambda = lambda, iter = iter, ELBO = LBnew))
}