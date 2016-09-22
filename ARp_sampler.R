library(GGally)
library(coda)
library(ggplot2)
library(gridExtra)
library(reshape)

T = 100    #Length of y
N = 5000   #MCMC draws

#True parameters
sigma2 <- 2
phi1 <- 0.4
phi2 <- -0.2

a = 1
b = 1
shape = T/2 -1 + a


#Generate data
y = c(0, 0.5)
for(i in 3:T){
  y = c(y, phi1*y[(i-1)] + phi2*y[(i-2)] + rnorm(1, 0, sqrt(sigma2)))
}

genRho = function(invroots) {
  roots = 1/invroots
  lr = roots[Im(roots) == 0]
  pr = length(lr)
  lc = roots[Im(roots) != 0]
  pc = length(lc)
  if(pc > 0) {
    lc = matrix(c(Re(lc[1]), Im(lc[1])), ncol = 2)
  }
  p = pr + pc
  Psi = matrix(0, ncol = p, nrow = p+1)
  Psi[1, ] = 1
  if(pr > 0){
    Psi[2, 1] = -lr[1]
    if(pr > 1) {
      for(i in 2:pr){
        Psi[2:(i+1), i] = Psi[2:(i+1), i-1] - lr[i]*Psi[1:i, i-1]
      }
    }
  }
  if(pc > 0) {
    if(pr > 0) {
      Psi[2, pr + 2] = -2*lc[1] + Psa[2, pr]
      Psi[3:(pr+3), pr+2] = (lc[1]^2 + lc[2]^2)*Psi[1:(pr+1), pr]-
        2*lc[1]*Psi[2:(pr+2), pr] + Psi[3:(pr+3), pr]
    } else {
      Psi[2, 2] = -2*lc[1]
      Psi[3, 2] = (lc[1]^2 + lc[2]^2)
      }
    if(pc > 2) {
      for(i in seq(4, pc, 2)) {
        pri = pr + i
        prim = pri - 2
        Psi[2, pri] = -2*lc[i-1] + Psi[2, prim]
        Psi[3:(pri+1), pri] = (lc[i-1]^2 + lc[i]^2)*Psi[1:(pri-1), prim]-
          2*lc[i-1]*Psi[2:pri, prim] + Psi[3:(pri+1), prim]
      }
    }
  }
  Rho = Psi[1:(p+1), p]
  Rho[(p+1):1] = Rho[1:(p+1)]/Rho[p+1]
  return(c(-Rho[2], -Rho[3]))
}

loglike = function(y, phi1, phi2, sig2) {
  loglike = -(T-p)/2*log(sig2) - sum((y[3:T] - phi1*y[2:(T-1)] - phi2*y[1:(T-2)])^2)/(2*sig2)
  return(loglike)
}

gen.invroots = function(invroots) {
  if(Im(invroots[1]) == 0) {
    a = runif(1)
    if(a < 0) {
      invroots[1] = runif(1, -1, 1)
    } else {
      invroots[2] = runif(1, -1, 1)
    }
  } else {
    theta = 2*pi*runif(1)
    rho = sqrt(runif(1))
    a = rho*cos(theta)
    b = rho*sin(theta)
    invroots[1] = complex(real = a, imaginary = b)
    invroots[2] = complex(real = a, imaginary = - b)
  }
  return(invroots)
}

gen.invroots2 = function(invroots) {
  if(Im(invroots[1]) != 0) {
    invroots[1] = runif(1, -1, 1)
    invroots[2] = runif(1, -1, 1)
  } else {
    theta = 2*pi*runif(1)
    rho = sqrt(runif(1))
    a = rho*cos(theta)
    b = rho*sin(theta)
    invroots[1] = complex(real = a, imaginary = b)
    invroots[2] = complex(real = a, imaginary = - b)
  }
  return(invroots)
}

theta = matrix(0, ncol = 3, nrow = N)
theta[1, 3] = 1
theta[1, 1:2] = 0.3

accept1 = 0
accept2 = 0

invroots = 1/polyroot(c(1, -theta[1, 1], - theta[1, 2]))

for(i in 2:N){
  #First step
  newroots = gen.invroots(invroots)
  candidate1 = genRho(newroots)
  loglike.candidate = loglike(y, candidate1[1], candidate1[2], theta[i-1, 3])
  loglike.old = loglike(y, theta[i-1, 1], theta[i-1, 2], theta[i-1, 3])
  alpha = min(1, exp(loglike.candidate-loglike.old))
  u = runif(1)
  if(u <= alpha) {             #accept draw
    eps = newroots
    accept1 = accept1 + 1
  } else {                    #reject draw
    eps = invroots
  }
  
  #Second step
  newroots = gen.invroots(eps)
  candidate2 = genRho(newroots)
  loglike.candidate2 = loglike(y, candidate2[1], candidate2[2], theta[i-1, 3])
  alpha = min(1, exp(loglike.candidate2-loglike.old))
  u = runif(1)
  if(u <= alpha) {             #accept draw
    invroots = newroots
    theta[i, 1:2] = candidate2
    accept2 = accept2 + 1
  } else {                    #reject draw
    eps = invroots
    theta[i, 1:2] = theta[i-1, 1:2]
  }
  
  inv.scale = b + sum(((y[3:T] - theta[i, 1]*y[2:(T-1)]- theta[i, 2]*y[1:(T-2)])^2))/2
  theta[i, 3] = 1/rgamma(1, shape = shape, rate = inv.scale)
}