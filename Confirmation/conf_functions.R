#MCMC
posterior.statistics = function(x){
  l95 = quantile(x, probs = 0.025)
  mean = mean(x)
  median = median(x)
  u95 = quantile(x, probs = 0.975)
  return(c(l95, mean, median, u95))
}

lphidens = function(y, phi1, phi2, sigma2){
  T = length(y)
  rho0 = (1-phi2) / ((1+phi2) * ((1-phi2)^2 - phi1^2))
  rho1 = phi1 / ((1+phi2) * ((1-phi2)^2 - phi1^2))
  p1 = -1/2 * log(rho0^2 - rho1^2)
  p2 = -1 / (2*sigma2) * (rho0*(y[1]^2 + y[2]^2) - 2*rho1*y[1]*y[2]) / (rho0^2 - rho1^2)
  p3 = -1 / (2*sigma2) * sum((y[3:T] - phi1 * y[2:(T-1)] - phi2*y[1:(T-2)])^2)
  return(p1 + p2 + p3)
}

lqdens = function(x, mean, sd, otherphi, whichphi){
  density = dnorm(x, mean, sd)
  if(whichphi == 1){
    range = pnorm(1 - otherphi, mean, sd) - pnorm(otherphi - 1, mean, sd)
  } else {
    range = pnorm(1 - abs(otherphi), mean, sd) - pnorm(-1, mean, sd)
  }
  return(log(density/range))
}

#SVB
simr = function(n, lambda, transform = FALSE){
  eps = matrix(0, nrow = n, ncol = 3)
  theta = matrix(0, nrow = n, ncol = 3)
  for(i in 1:n) {
    flag = FALSE
    while(!flag) {
      eps[i, 1:2] = rnorm(2)
      eps[i, 3] = runif(1)
      theta[i, 1] = lambda[1] + lambda[3]*eps[i, 1]
      theta[i, 2] = lambda[2] + lambda[4]*eps[i, 1] + lambda[5]*eps[i, 2]
      theta[i, 3] = qigamma(eps[i, 3], lambda[6], lambda[7])
      if(theta[i, 2] > -1 & theta[i, 2]  < 1 + theta[i, 1] & theta[i, 2] < 1 - theta[i, 1]){
        flag = TRUE
      }
    }
  }
  if(transform){
    return(theta)
  } else {
    return(eps)
  }
}

ELBOr = function(lambda, y, n = 1000){
  theta = simr(n, lambda, TRUE)
  out = rep(0, n)
  for(i in 1:n){
    a = logjointc(y, theta[i, 1], theta[i, 2], theta[i, 3])
    b = - logqc(theta[i, ], lambda)
    out[i] = a + b
  }
  return(mean(out))
}

allderivr = function(lambda, y, param){
  epsilon = simr(1, lambda)
  theta = rep(0, 3)
  theta[1] = lambda[3]*epsilon[1] + lambda[1]
  theta[2] = lambda[4]*epsilon[1] + lambda[5]*epsilon[2] + lambda[2]
  theta[3] = qigamma(epsilon[3], lambda[6], lambda[7])

  if(param == 1){
    deriv1 = phideriv(y, theta[1], theta[2], theta[3], 1)
  } else if(param == 2) {
    deriv1 = phideriv(y, theta[1], theta[2], theta[3], 2)
  } else if(param == 3) {
    deriv1 = epsilon[1] * phideriv(y, theta[1], theta[2], theta[3], 1)
  } else if(param == 4) {
    deriv1 = epsilon[1] * phideriv(y, theta[1], theta[2], theta[3], 2)
  } else if(param == 5) {
    deriv1 = epsilon[2] * phideriv(y, theta[1], theta[2], theta[3], 2)
  }

  h = 0.000001
  lambda2 = lambda
  lambda2[param] = lambda2[param] + h
  
  if(param == 6 | param == 7){
    sigmasq2 = qigamma(epsilon[3], lambda2[6], lambda2[7])
    dSigmasqDLambda = (sigmasq2 - theta[3])/h
    deriv1 = dSigmasqDLambda *  sigderiv(y, theta[1], theta[2], theta[3])
  }
  
  theta2 = rep(0, 3)
  theta2[1] = lambda2[3]*epsilon[1] + lambda2[1]
  theta2[2] = lambda2[4]*epsilon[1] + lambda2[5]*epsilon[2] + lambda2[2]
  theta2[3] = qigamma(epsilon[3], lambda2[6], lambda2[7])
  deriv2 = (logqc(theta2, lambda2)-logqc(theta, lambda))/h
  return(deriv1 - deriv2)
}


SGAr = function(y, lambda, s, threshold, M, eta1, eta2){
  Gt = rep(0, 7)
  eta = c(rep(eta1, 5), rep(eta2, 2))
  k = 0
  LBnew = ELBOr(lambda, y)
  diff = threshold + 1
  pt = rep(0, 7)
  
  while(diff > threshold | k < 50){ 
    if(k >= M){
      break
    }
    partials = rep(0, 7)
    for(i in 1:s){
      for(j in 1:7){
        partials[j] = partials[j] + allderivr(lambda, y, j)/s
      }
    }
    for(j in 1:7){
      Gt[j] = Gt[j] + partials[j]^2
      pt[j] = eta[j] * Gt[j]^(-0.5)
      if(k > 0){
        lambda[j] = lambda[j] + pt[j] * partials[j] ##in step one the product of pt and partials is eta. so wait two steps for pt to settle a little
      }
    }
    LBold = LBnew
    LBnew = ELBOr(lambda, y)
    diff = abs(LBnew - LBold)
    k = k + 1
  }
  out = c(lambda, k, LBnew)
  return(out)
}