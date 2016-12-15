invert_cdf = function(u, lambda, i){
  if(i == 1 | i == 2){
    support = seq(0, 1, 0.001)
    dens = dnormmix(support, lambda[(5*i-4):(5*i)])
    CDF = cumsum(dens)
    u = u * max(CDF)
    draws = rep(0, length(u))
    for(j in 1:length(u)){
      draws[i] = support[min(which(CDF > u))]
    }
  } else if(i == 3) {
    draws = qigamma(u, lambda[11], lambda[12])
  } else {
    draws = qnorm(u, lambda[5+2*i], lambda[6+2*i])
  }
  return(draws)
}

dnormmix = function(theta, lambda){
  dens = lambda[5] * dnorm(theta, lambda[1], sqrt(lambda[2])) + 
    (1-lambda[5]) * dnorm(theta, lambda[3], sqrt(lambda[4]))
  return(dens)
}

simtheta = function(s, RVM, lambda){
  #Simulate theta - R only
  #Require RVM Object created from Vine Structure (prespecified) and dependence parameters (eta)
  simvine = RVineSim(s, RVM)
  #Require Inverse CDF transformation for all marginals from mean field parameters (lambda)
  sim = matrix(0, s, 10)
  for(i in 1:10){
    sim[,i] = invert_cdf(simvine[,s], lambda[[s]], i)
  }
  return(sim)
}

ELBO = function(y, lambda, eta, RVM, reps = 10){
  theta = simtheta(reps, RVM, lambda)
  output = 0;
  for(i in 1:reps){
    output = output +  (log_joint_dens(y, theta[i,]) - logq(theta[i,], lambda, eta, RVM))/reps;
  }
  return(output)
}

qnormmix = function(theta, lambda){
  support = seq(0, 1, 0.001)
  dens = lambda[5]*dnorm(support, lambda[1], sqrt(lambda[2])) + (1-lambda[5])*dnorm(support, lambda[3], sqrt(lambda[4]))
  cdf = cumsum(dens)/1000
  u = cdf[min(which(theta < support))]
  return(u)
}

logq = function(theta, lambda, eta, Vine){
  require(pscl)
  require(VineCopula)
  marginals = rep(0, 10)
  copulas = rep(0, 37)
  for(i in 1:10){
    if(i == 1 | i == 2){
      marginals[i] =  log(lambda[5*i] * dnorm(theta[i], lambda[5*(i-1) + 1], lambda[5*(i-1) + 2]) + 
        (1-lambda[5*i]) * dnorm(theta[i], lambda[5*(i-1) + 3], lambda[5*(i-1) + 4]))
    } else if (i == 3) {
      marginals[i] = log(densigamma(theta[i], lambda[11], lambda[12]))
    } else {
      marginals[i] = dnorm(theta[i], lambda[5 + 2*i], lambda[6 + 2*i], log = TRUE)
    }
  }
  counter = 1
  for(i in 1:10){
    for(j in 1:10){
      if(Vine$family[i, j] != 0){
        t1 = Vine$Matrix[j, j]
        t2 = Vine$Matrix[i, j]
        if(t1 <= 2){
          u1 = qnormmix(theta[t1], lambda[(t1*5-4):(t1*5)])
        } else if(t1 == 3){
          u1 = qigamma(theta[t1], lambda[11], lambda[12])
        } else {
          u1 = qnorm(theta[t1], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])) 
        }
        if(t2 <= 2){
          u2 = qnormmix(theta[t2], lambda[(t1*5-4):(t1*5)])
        } else if(t1 == 3){
          u2 = qigamma(theta[t2], lambda[11], lambda[12])
        } else {
          u2 = qnorm(theta[t2], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])) 
        }
        copulas[counter] = log(BiCopPDF(u1, u2, Vine$family[i, j], eta[[1]][i, j], eta[[2]][i, j]))
        counter = counter + 1
      }
    }
  }
  return(sum(marginals) + sum(copulas))
}

deriv_lambda = function(theta, lambda, eta, y, j, Vine){
  h = 0.000001
  if(j <= 10){
    dens = ceiling(j/5)
  } else {
    dens = ceiling(j/2) -3
  }
  lambda2 = lambda
  lambda2[j] = lambda2[j] + h
  if(dens <= 2){
    deriv_marginal = (dnormmix(theta[dens], lambda2[(dens*5-4):(dens*5)]) - dnormmix(theta[dens], lambda[(dens*5-4):(dens*5)]))/h
    deriv_cdf = (qnormmix(theta[dens], lambda2[(dens*5-4):(dens*5)]) - qnormmix(theta[dens], lambda[(dens*5-4):(dens*5)]))/h
  } else if(dens = 3){
    deriv_marginal = (densigamma(theta[dens], lambda2[11], lambda2[12]) - densigamma(theta[dens], lambda[11], lambda[12]))/h
    deriv_cdf = (qigamma(theta[dens], lambda2[11], lambda2[12]) - qigamma(theta[dens], lambda[11], lambda[12]))/h
  } else {
    deriv_marginal = (dnorm(theta[dens], lambda2[5 + 2*i], sqrt(lambda2[6 + 2*i])) - dnorm(theta[dens], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])))/h
    deriv_cdf = (qnorm(theta[dens], lambda2[5 + 2*i], sqrt(lambda2[6 + 2*i])) - qnorm(theta[dens], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])))/h
  }
  cop_derivs = NULL
  for(i in 1:10){
    if(Vine$Matrix[i, i] == dens){
      t1 = dens
      for(j in (i+1):10){
        if(Vine$family[i, j] != 0){
          t2 = Vine$Matrix[i, j]
          if(t1 <= 2){
            u1 = qnormmix(theta[t1], lambda[(t1*5-4):(t1*5)])
          } else if(t1 == 3){
            u1 = qigamma(theta[t1], lambda[11], lambda[12])
          } else {
            u1 = qnorm(theta[t1], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])) 
          }
          if(t2 <= 2){
            u2 = qnormmix(theta[t2], lambda[(t1*5-4):(t1*5)])
          } else if(t1 == 3){
            u2 = qigamma(theta[t2], lambda[11], lambda[12])
          } else {
            u2 = qnorm(theta[t2], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])) 
          }
          cop_derivs = c(cop_derivs, BiCopDeriv(u1, u2, Vine$family[i, j], eta[[1]][i, j], eta[[2]][i, j], deriv = "u1"))
        }
      }
    }
  }
  for(j in 1:10){
    for(i in (j+1):10){
      if(Vine$Matrix[i, j] == dens){
        if(Vine$family[i, j] != 0){
          t1 = Vine$Matrix[j, j]
          t2 = dens
          if(t1 <= 2){
            u1 = qnormmix(theta[t1], lambda[(t1*5-4):(t1*5)])
          } else if(t1 == 3){
            u1 = qigamma(theta[t1], lambda[11], lambda[12])
          } else {
            u1 = qnorm(theta[t1], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])) 
          }
          if(t2 <= 2){
            u2 = qnormmix(theta[t2], lambda[(t1*5-4):(t1*5)])
          } else if(t1 == 3){
            u2 = qigamma(theta[t2], lambda[11], lambda[12])
          } else {
            u2 = qnorm(theta[t2], lambda[5 + 2*i], sqrt(lambda[6 + 2*i])) 
          }
          cop_derivs = c(cop_derivs, BiCopDeriv(u1, u2, Vine$family[i, j], eta[[1]][i, j], eta[[2]][i, j], deriv = "u2"))
        }
      }
    }
  }
  return(deriv_marginal + deriv_cdf*sum(cop_derivs))
}

deriv_eta = function(theta, lambda, eta, y, which, j, k, Vine){}

SGA_lambda = function(y, theta, lambda, eta, tuning, s, LB, Vine, k, M){
  diff = 1
  RVM = RVineMatrix(Vine$Matrix, Vine$family, eta[[1]], eta[[2]])
  r = 0
  Gt = rep(0, 26)
  while(diff > threshold | r < k){
    if(r > M){
      break
    }
    r = r + 1
    theta = simtheta(s, RVM, lambda)
    partials = rep(0, 26)
    for(i in 1:s){
      for(j in 1:26){
        partials[j] = deriv_lambda(theta[i,], lambda, eta, y, j)
      }
    }
    for(j in 1:26){
      Gt[j] = Gt[j] + partials[j]^2
      pt = tuning * Gt[j]^(-0.5)
      lambda[j] = lambda[j] + pt * partials[j]
    }
    LBnew = ELBO(y, lambda, eta, RVM)
    diff = abs(LBnew - LB)
    LB = LBnew
  }
  return(c(lambda, LB))
}

SGA_eta = function(y, theta, lambda, eta, tuning, s, LB, Vine, k, M){
  diff = 1
  r = 0
  Gt = array(0, dim = c(10, 10, 2))
  while(diff > threshold | r < k){
    if(r > M){
      break
    }
    r = r + 1
    RVM = RVineMatrix(Vine$Matrix, Vine$family, eta[[1]], eta[[2]])
    theta = simtheta(s, RVM, lambda)
    partials = array(0, dim = c(10, 10, 2))
    for(i in 1:s){
      for(j in 1:10){
        for(k in 1:10){
          for(l in 1:2){
            if(eta[[l]][j, k] != 0){
              partials[j, k, l] = deriv_eta(theta[i,], lambda, eta, y, l, j, k, Vine)
            }
          }
        }
      }
    }
    for(j in 1:10){
      for(k in 1:10){
        for(l in 1:2){
          if(eta[[l]][j, k] != 0){
            Gt[j, k, l] = Gt[j, k, l] + partials[j, k, l]^2
            pt = tuning * Gt[j, k, l]^(-0.5)
            eta[[l]][j, k] = eta[[l]][j, k] + pt * partials[j, k, l]
          }
        }
      }
    }
    LBnew = ELBO(y, lambda, eta, RVM)
    diff = abs(LBnew - LB)
    LB = LBnew
  }
  eta[[3]] = LB
  return(eta)
}


