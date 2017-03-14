ELBO = function(Vine, lambda, meanX, varX, y, epsilon, n = 1000){
  out = rep(0, n)
  for(i in 1:n){
    theta = c(qtruncnorm(epsilon[i, 1], -1, 1, lambda[1], sqrt(lambda[2])),
              QRM::qst(epsilon[i, 2], lambda[3], sqrt(lambda[4]), lambda[5]),
              qigamma(epsilon[i, 3], lambda[6], lambda[7]),
              qigamma(epsilon[i, 4], lambda[8], lambda[9]))
    x = rep(0, 101)
    for(j in 1:101){
      x[j] = qnorm(epsilon[i, 4+j], meanX[j], sqrt(varX[j]))
    }
    a = logjoint(y, x, theta)
    b = logq(T, theta, lambda, meanX, varX, Vine, x, epsilon[i,])
    out[i] = a - b
  }
  return(mean(out))
}

dent = function(x, mean, var, df, log = FALSE){
  front = gamma((df + 1)/ 2) / (gamma(1/2) * gamma(df/2))
  middle = (df * var)^(-1/2)
  end = (1 + 1/df * (x - mean)^2/var)^(-(df+1)/2)
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

logCopulaDeriv = function(eps1, eps2, family, par1, par2){
  h = 0.000001
  fxh = log(BiCopPDF(eps1 + h, eps2, family, par1, par2))
  fx = log(BiCopPDF(eps1, eps2, family, par1, par2))
  return((fxh - fx)/h)
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
  
  h = 0.0000001
  if(param == 1 | param == 2) {
    var = 1
    lambda2 = lambda 
    lambda2[param] = lambda[param] + h
    marginalDeriv = (log(dtruncnorm(theta[1], -1, 1, lambda2[1], sqrt(lambda2[2]))) - 
                       log(dtruncnorm(theta[1], -1, 1, lambda[1], sqrt(lambda[2]))))/h
    CdfDeriv = (log(ptruncnorm(theta[1], -1, 1, lambda2[1], sqrt(lambda2[2]))) - 
                  log(ptruncnorm(theta[1], -1, 1, lambda[1], sqrt(lambda[2]))))/h
  } else if(param == 3){
    var = 2
    marginalDeriv = -(2*lambda[5] + 1)*(lambda[3] - theta[2]) / (lambda[4]*lambda[5] + (lambda[3] - theta[2])^2)
    CdfDeriv = - dent(theta[2], lambda[3], lambda[4], lambda[5])
  } else if(param == 4){
    var = 2
    marginalDeriv = lambda[5]*(lambda[4] + 2*(theta[2] - lambda[3])^2) / (2*lambda[4] * (lambda[5]*lambda[4] + (theta[2] - lambda[3])^2))
    CdfDeriv = -1/2 *  (theta[2] - lambda[3]) / lambda[4] * dent(theta[2], lambda[3], lambda[4], lambda[5])
  } else if(param == 5){
    var = 2
    a = (theta[2] - lambda[3])^2
    marginalDeriv = a * (lambda[5] - 1/2) / (lambda[5] * (a + lambda[4]*lambda[5])) -
      log(a / (lambda[4]*lambda[5]) + 1) - 1 / (2*lambda[5]) +
      1/2 * (digamma((lambda[5] + 1)/2) - digamma(lambda[5]/2))
    CdfDeriv = (pt((theta[2] - lambda[3]) / sqrt(lambda[4]), lambda[5] + h) - 
                  pt((theta[2] - lambda[3]) / sqrt(lambda[4]), lambda[5]))/h
  } else if(param == 6 | param == 8){
    var = param/2
    marginalDeriv = log(lambda[param + 1]) - digamma(lambda[param]) -log(theta[var])
    CdfDeriv = (pigamma(theta[var], lambda[param] + h, lambda[param + 1]) - 
                  pigamma(theta[var], lambda[param], lambda[param + 1]))/h
  } else if(param == 7 | param == 9){
    var = (param - 1)/2
    marginalDeriv = lambda[param - 1]/lambda[param] - 1/theta[var]
    CdfDeriv = (pigamma(theta[var], lambda[param - 1], lambda[param]+ h) - 
                  pigamma(theta[var], lambda[param - 1], lambda[param]))/h
  } else if(param >= 10 & param <= 110){
    var = param - 9
    marginalDeriv = (meanX[var] - x[var]) / varX[var]
    CdfDeriv = - dnorm(x[var], meanX[var], sqrt(varX[var]))
  } else {
    var = param - 110
    marginalDeriv = ((x[var] - meanX[var])^2 - varX[var]) / (2*varX[var]^2)
    CdfDeriv = - 1/2 * (x[var] - meanX[var]) / varX[var] * dnorm(x[var], meanX[var], sqrt(varX[var]))
  }
  
  if(var <= 4){
    CopDeriv = 0
    for(i in (1:105)[-(106 - var)]){
      family = Vine$family[(106 - var),i]
      if(family != 0){
        CopDeriv = CopDeriv + logCopulaDeriv(epsilon[var], epsilon[Vine$Matrix[i, i]], family, 
                                             Vine$par[(106 - var), i], Vine$par2[(106 - var), i])
      }
    }
  } else if(var == 5 | var == 105){
    k = 3 - max(1, var - 103)
    family = Vine$family[101, k]
    CopDeriv = logCopulaDeriv(epsilon[var], epsilon[ifelse(k == 1, var - 1, var + 1)], family, Vine$par[101, k], Vine$par2[101, k])
    for(i in 102:105){
      family = Vine$family[i, k]
      if(family != 0){
        CopDeriv = CopDeriv + logCopulaDeriv(epsilon[var], epsilon[106 - var], family, Vine$par[i, k], Vine$par2[i, k])
      }
    }
  } else {
    family = Vine$family[101, var - 3]
    CopDeriv = logCopulaDeriv(epsilon[var], epsilon[var+1], family, Vine$par[101, var - 3], Vine$par2[101, var - 3])
    family = Vine$family[101, var - 4]
    CopDeriv = CopDeriv + logCopulaDeriv(epsilon[var], epsilon[var - 1], family, Vine$par[101, var - 4], Vine$par2[101, var - 4])
    for(i in 102:105){
      family = Vine$family[i, var - 3]
      if(family != 0){
        CopDeriv = CopDeriv + logCopulaDeriv(epsilon[var], epsilon[106 - var], family, Vine$par[i, var - 3], Vine$par2[i, var - 3])
      }
    }
  }
  
  logqDeriv = marginalDeriv + CdfDeriv * CopDeriv
  logqEval = logq(T, theta, lambda, meanX, varX, Vine, x, epsilon)
  logJointEval = logjoint(y, x, theta)
  
  return(logqDeriv * (logqEval - logJointEval))
  
}

copulaDeriv = function(y, Vine, lambda, meanX, varX, epsilon, i, j, whichpar){
  T = length(y)
  theta = c(qtruncnorm(epsilon[1], -1, 1, lambda[1], sqrt(lambda[2])),
            QRM::qst(epsilon[2], lambda[3], sqrt(lambda[4]), lambda[5]),
            qigamma(epsilon[3], lambda[6], lambda[7]),
            qigamma(epsilon[4], lambda[8], lambda[9]))
  x = rep(0, 101)
  for(m in 1:101){
    x[m] = qnorm(epsilon[4+m], meanX[m], sqrt(varX[m]))
  }
  
  family = Vine$family[i, j]
  u1 = Vine$Matrix[j, j]
  u2 = Vine$Matrix[i, j]
  wrt = "par"
  if(whichpar == 2){
    wrt = "par 2"
  }
  par = Vine$par[i, j]
  par2 = Vine$par2[i, j]
  
  if(family %in% c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)){
    CopDeriv = BiCopDeriv(epsilon[u1], epsilon[u2], family, par, par2, wrt, log = TRUE)
  } else {
    h = 0.000001
    if(whichpar == 1){
      parh = par + h
      CopDeriv = (log(BiCopPDF(epsilon[u1], epsilon[u2], family, parh, par2)) - 
                    log(BiCopPDF(epsilon[u1], epsilon[u2], family, par, par2)))/h
    } else {
      par2h = par2 + h
      CopDeriv = (log(BiCopPDF(epsilon[u1], epsilon[u2], family, par, par2h)) - 
                    log(BiCopPDF(epsilon[u1], epsilon[u2], family, par, par2)))/h
    }
  }
  logqEval = logq(T, theta, lambda, meanX, varX, Vine, x, epsilon)
  logJointEval = logjoint(y, x, theta)
  
  return(CopDeriv * (logqEval - logJointEval))
}

#lambda = 1. Phi Mean, 2. Phi Variance, 3. Mu mean, 4. Mu Variance, 5. Mu df, 6. Sigy Alpha, 7. Sigy Beta, 8. Sigx Alpha, 9. Sigx Beta
SGA = function(y, lambda, meanX, varX, S, M, threshold, eta, Vine, n){
  k = length(lambda)
  T = length(y)
  Gt = rep(0, 209)
  GtCop = matrix(0, ncol = 104, nrow = 5)
  GtCop2 = matrix(0, ncol = 104, nrow = 5)
  iter = 0
  
  epsilon = RVineSim(n, Vine)
  LBnew = ELBO(Vine, lambda, meanX, varX, y, epsilon, n)
  diff = threshold + 1
  
  while(diff > threshold | iter < 10){
    if(iter >= M){
      break
    }
    
    epsSample = sample(1:n, S, replace = TRUE)
    partials = rep(0, 209)
    partialsCop = matrix(0, ncol = 104, nrow = 5)
    partialsCop2 = matrix(0, ncol = 104, nrow = 5)
    for(s in 1:S){
      for(j in 1:209){
        partials[j] = partials[j] + ELBOderiv(y, Vine, lambda, meanX, varX, epsilon[epsSample[s],], j)/S
      }
      
      for(j in 1:104){
        for(i in 101:105){
          if(Vine$par[i, j] != 0){
            partialsCop[i-100, j] = partialsCop[i-100, j] + copulaDeriv(y, Vine, lambda, meanX, varX, epsilon[epsSample[s],], i, j, 1)/S
          }
          if(Vine$par2[i, j] != 0){
            partialsCop2[i-100, j] = partialsCop2[i-100, j] + copulaDeriv(y, Vine, lambda, meanX, varX, epsilon[epsSample[s],], i, j, 2)/S
          }
        }
      }
    }
    
    Gt = Gt + partials^2
    GtCop = GtCop + partialsCop^2
    GtCop2 = GtCop2 + partialsCop2^2
    pt = eta * Gt^(-1/2)
    ptCop = eta * GtCop^(-1/2)
    ptCop2 = eta * GtCop2^(-1/2)
    if(iter > 0){
      c(lambda, meanX, varX) = c(lambda, meanX, varX) + pt * partials
      Vine$par[101:105, 1:104] = Vine$par[101:105, 1:104] + ptCop * partialsCop 
      Vine$par2[101:105, 1:104] = Vine$par2[101:105, 1:104] + ptCop2 * partialsCop2
    }
    
    LBold = LBnew
    epsilon = RVineSim(n, Vine)
    LBnew = ELBO(Vine, lambda, meanX, varX, y, n)
    diff = abs(LBold - LBnew)
    iter = iter + 1
  }
  return(list(lambda = lambda, meanX = meanX, varX = varX, Vine = Vine, iter = iter, ELBO = LBnew))
}

