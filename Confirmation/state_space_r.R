simtheta = function(Vine, lambda, meanX, varX){
  #simulate cdf via Vine
  epsilon = RVineSim(1, Vine)
  #Transform (0, 1) eps to theta and states
  theta = c(qtruncnorm(epsilon[1], -1, 1, lambda[1], sqrt(lambda[2])),
            QRM::qst(epsilon[2], lambda[3], sqrt(lambda[4]), lambda[5]),
            qigamma(epsilon[3], lambda[6], lambda[7]),
            qigamma(epsilon[4], lambda[8], lambda[9]))
  x = rep(0, 101)
  for(i in 1:101){
    x[i] = qnorm(epsilon[4+i], meanX[i], sqrt(varX[i]))
  }
  
  return(list(theta = theta, x = x))
}
  
ELBO = function(Vine, lambda, meanX, varX, y, epsilon, n = 1000){
  out = rep(0, n)
  for(i in 1:n){
    theta = c(qtruncnorm(epsilon[i, 1], -1, 1, lambda[1], sqrt(lambda[2])),
              QRM::qst(epsilon[i, 2], lambda[3], sqrt(lambda[4]), lambda[5]),
              qigamma(epsilon[i, 3], lambda[6], lambda[7]),
              qigamma(epsilon[i, 4], lambda[8], lambda[9]))
    x = rep(0, 101)
    for(j in 1:101){
      x[i] = qnorm(epsilon[i, 4+j], meanX[j], sqrt(varX[j]))
    }
    a = logjoint(y, x, theta)
    b = logq(T, theta, lambda, meanX, varX, Vine, x)
    out[i] = a - b
  }
  return(mean(out))
}

dent = function(x, mean, var, df, log = FALSE){
  front = gamma((df + 1)/ 2) / (gamma(1/2) * gamma(df/2))
  middle = (df * var)^(-1/2)
  end = (1 + 1/df * ((x - mean)/sqrt(var))^2)^(-(v+1)/2)
  if(log){
    return(log(front) + log(middle) + log(end))
  } else {
    return(front*middle*end)
  }
}

logq = function(T, theta, lambda, meanX, varX, Vine, x, epsilon){
  phi = log(dtruncnorm(theta[1], -1, 1, lambda[1], sqrt(lambda[2])))
  mu = dent(theta[2], lambda[3], lambda[4], lambda[5], log = TRUE)
  sigy = log(densigamma(theta[3], lambda[6], lambda[7]))
  sigx = log(densigamma(theta[4], lambda[8], lambda[9]))
  states = 0
  for(i in 1:(T+1)){
    states = states + dnorm(x[i], meanX[i], sqrt(varX[i]), log = TRUE)
  }
  
  copulas = 0
  for(j in 1:104){
    for(i in 101:105){
      if(Vine$family[i, j] != 0){
        copulas = copulas + log(BiCopPDF(epsilon[Vine$Matrix[j,j]], epsilon[Vine$Matrix[i, j]], 
                                    Vine$family[i, j], Vine$par[i, j], Vine$par2[i, j]))
      }
    }
  }
  
  return(phi + mu + sigy + sigx + states + copulas)
}

ELBOderiv = function(y, Vine, lambda, meanX, varX, epsilon, param){
  T = length(y)
  theta = c(qtruncnorm(epsilon[1], -1, 1, lambda[1], sqrt(lambda[2])),
            QRM::qst(epsilon[2], lambda[3], sqrt(lambda[4]), lambda[5]),
            qigamma(epsilon[3], lambda[6], lambda[7]),
            qigamma(epsilon[4], lambda[8], lambda[9]))
  x = rep(0, 101)
  for(j in 1:101){
    x[j] = qnorm(epsilon[4+j], meanX[j], sqrt(varX[j]))
  }
  
  h = 0.000000001
  if(param == 1 | param == 2) {
    lambda2 = lambda 
    lambda2[param] = lambda[param] + h
    marginalDeriv = (log(dtruncnorm(theta[1], -1, 1, lambda2[1], sqrt(lambda2[2]))) - 
                       log(dtruncnorm(theta[1], -1, 1, lambda[1], sqrt(lambda[2]))))/h
  } else if(param == 3){
    marginalDeriv = -(2*lambda[5] + 1)*(lambda[3] - theta[2]) / (lambda[4]*lambda[5] + (lambda[3] - theta[2])^2)
  } else if(param == 4){
    marginalDeriv = lambda[5]*(lambda[4] + 2*(theta[2] - lambda[3])^2) / (2*lambda[4] * (lambda[5]*lambda[4] + (theta[2] - lambda[3])^2))
  } else if(param == 5){
    a = (theta[2] - lambda[3])^2
    marginalDeriv = a * (lambda[5] - 1/2) / (lambda[5] * (a + lambda[4]*lambda[5])) -
      log(a / (lambda[4]*lambda[5]) + 1) - 1 / (2*lambda[5]) +
      1/2 * (digamma((lambda[5] + 1)/2) - digamma(lambda[5]/2))
  }
    
  } else if(param >= 10 & param <= 110){
    marginalDeriv = (meanX[param - 9] - x[param - 9]) / varX[param - 9]
  } else if(param >= 111 & param <= 211){
    marginalDeriv = ((x[param - 110] - meanX[param - 110])^2 - varX[param - 110]) / (2*varX[param - 110]^2)
  }
  
  #if/else for params
  
  #return deriv(logq) * (logjoint - logq)
}

#lambda = 1. Phi Mean, 2. Phi Variance, 3. Mu mean, 4. Mu Variance, 5. Mu df, 6. Sigy Alpha, 7. Sigy Beta, 8. Sigx Alpha, 9. Sigx Beta
SGA = function(y, lambda, meanX, varX, s, M, threshold, eta, Vine, n){
  k = length(lambda)
  T = length(y)
  Gt = rep(0, 786)
  iter = 0
  
  epsilon = RVineSim(n, Vine)
  LBnew = ELBO(Vine, lambda, meanX, varX, y, epsilon, n)
  diff = threshold + 1
  Par1 = matrix(0, T + 5, T + 5)
  Par2 = matrix(0, T + 5, T + 5)
  
  while(diff > threshold | iter < 10){
    if(iter >= M){
      break
    }
    
    epsilon = RVineSim(n, Vine)
    epsSample = sample(1:n, s, replace = TRUE)
    partials = rep(0, 786)
    for(i in 1:s){
      for(j in 1:786){
        partials[j] = partials[j] + ELBOderiv(y, Vine, lambda, meanX, varX, epsilon[epsSample[i]], j)/s
      }
    }
    
    Gt = Gt + partials^2
    pt = eta * Gt^(-1/2)
    if(iter > 0){
      c(lambda, meanX, varX) = c(lambda, meanX, varX) + pt[1:211] * partials[1:211]
      Par1[Vine$par != 0] = pt[212:631] * partials[212:631]
      Par2[Vine$par2 != 0] = pt[632:786] * partials[632:786]
      Vine$par = Vine$par + Par1
      Vine$par2 = Vine$par2 + Par2
    }
    
    LBold = LBnew
    LBnew = ELBO(Vine, lambda, meanX, varX, y, n)
    diff = abs(LBold - LBnew)
    iter = iter + 1
  }
  return(list(lambda = lambda, meanX = meanX, varX = varX, Vine = Vine, iter = iter, ELBO = LBnew))
}