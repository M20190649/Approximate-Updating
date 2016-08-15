library(GGally)
library(coda)

T = 100    #Length of y
N = 50000   #MCMC draws

#True parameters
sigma2 <- 2
mu <- 2
phi <- 0.5

#Generate data
y = rnorm(1, mu, sqrt(sigma2/(1-phi^2)))
for(i in 2:T){
  y = c(y, mu + phi*(y[(i-1)]-mu) + rnorm(1, 0, sqrt(sigma2)))
}

#Prior parameters
mubar = 0
lamb = 10
a = 1
b = 1
p1 = 5
p2 = 5

#constant
shape = T/2 + a

#log phi kernel for the MH step
lphiden = function(draw, mu, sig2){
  t1 = 1/2*log(1-draw^2) -(1-draw^2)*(y[1]-mu)^2/(2*sig2)  #The normal exponential for y1
  t2 = -(sum(((y[2:T]-mu)-(draw*(y[1:(T-1)]-mu)))^2))/(2*sig2)    #The normal exponentials for the rest of the y's
  t3 = log(1+draw)*(p1-1) + log(1-draw)*(p2-1)   #The stretched beta prior
  return(t1+t2+t3)
}

#put draws here
theta = matrix(ncol=3, nrow=N)
theta[1,-2] = 1       #Starting parameters
theta[1, 2] = 0.5

#count accepted draws
accept = 0

for(i in 2:N){
  
  #mu conditional
  denom = lamb^2 * ((T-1)*(1-theta[(i-1),2])^2+(1-theta[(i-1),2]^2)) + theta[(i-1), 3]
  sigsqmu = theta[(i-1), 3]*lamb^2 / denom
  muhat = (lamb^2*((1-theta[(i-1),2]^2)*y[1] + (1-theta[(i-1),2])*sum(y[2:T]-theta[(i-1),2]*y[1:(T-1)])) + mubar*theta[(i-1), 3])/denom
  theta[i,1] = rnorm(1, muhat, sqrt(sigsqmu))
  
  #phi conditional
  phihat = sum((y[2:T]-theta[i, 1])*(y[1:(T-1)]-theta[i,1]))/sum((y[1:(T-1)]-theta[i,1])^2)   #from the Catherine's notes for ETC4541, pg 84
  vphi = theta[(i-1),3]/sum((y[1:(T-1)]-theta[i,1])^2)    #acceptance rate for this candidate ~75%, trying random walk instead
  phidraw = rnorm(1, phihat, sqrt(vphi))
  #phidraw = rnorm(1, theta[(i-1),2], 0.4)               #random walk MH candidate draw to tune acceptance ratios to ~25%
  
  #phi MH step
  if(phidraw > -1 & phidraw < 1) {                       #within stationary conditions
    mh.candidate = lphiden(phidraw, theta[i, 1], theta[(i-1), 3])     #symmetrical candidate so only need to evaluate the true density components
    mh.prev = lphiden(theta[(i-1), 2], theta[i, 1], theta[(i-1), 3])  #using log density
    alpha = min(1, exp(mh.candidate-mh.prev))
    u = runif(1)
    if(u <= alpha) {             #accept draw
      theta[i, 2] = phidraw
      accept = accept + 1
    } else {                    #reject draw
      theta[i, 2] = theta[(i-1), 2]
    }
  } else {                      #non stationary
    theta[i, 2] = theta[(i-1), 2]
  }
  
  #sig2 conditional
  inv.scale = b + ((y[1] - theta[i, 1])^2*(1-theta[i, 2]^2) + sum(((y[2:T] - theta[i, 1]) - theta[i, 2]*(y[1:(T-1)]-theta[i, 1]))^2))/2
  theta[i, 3] = 1/rgamma(1, shape = shape, rate = inv.scale)
}
accept/N
keep = theta[(N/5+1):N,] #Drop first 20% of draws
colnames(keep) = c("Mu", "Phi", "Sigma2")
effectiveSize(keep)
#ggpairs(keep)

thin = keep[seq(1, 4*N/5, 10),] #Keep every 10th
effectiveSize(thin)
ggpairs(thin)


### Copula VI
library(pscl)
library(mvtnorm)
library(copula)

#Copula Simulation function

cop.sim = function(lambda, eta, s) {
  P <- matrix(c(1, eta[1], eta[2],               # Correlation matrix
              eta[1], 1, eta[3],
              eta[2], eta[3], 1), nrow = 3)

## Simulation (compact vectorized version) 
  U <- pnorm(matrix(rnorm(s*3), ncol = 3) %*% chol(P))

  mu = qnorm(U[,1], mean = lambda[1], sd = sqrt(lambda[2]))
  phi = qnorm(U[,2], mean = lambda[3], sd = sqrt(lambda[4]))
  sig2 = qigamma(U[,3], alpha = lambda[5], beta = lambda[6])
  return(cbind(mu, phi, sig2))
}

logjoint = function(y, mu, phi, sig2) {
  t1 = 1/2*log(1-phi^2) - (T/2+a+1)*log(sig2) - (1-phi^2)*(y[1]-mu)^2/(2*sig2)
  t2 = - sum((y[2:T] - mu - phi*(y[1:(T-1)] - mu))^2)/(2*sig2)
  t3 = -b/sig2 - (mu-mubar)^2/(2*lamb^2)
  t4 = (p1-1)*log((1+phi)/2) + (p2-1)*log((1-phi)/2)
  return(t1 + t2 + t3 + t4)
}

logq = function(theta, lambda, eta) {
  t1 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  t2 = densigamma(theta[3], lambda[5], lambda[6])
  c.ms = normalCopula(eta[2])
  c.ps = normalCopula(eta[3])
  c1 = dCopula(c(pnorm(theta[1], mean = lambda[1], sd = sqrt(lambda[2])), pigamma(theta[3], alpha = lambda[5], beta = lambda[6])), c.ms, log = TRUE)
  c2 = dCopula(c(pnorm(theta[2], mean = lambda[3], sd = sqrt(lambda[4])), pigamma(theta[3], alpha = lambda[5], beta = lambda[6])), c.ms, log = TRUE)
  return(t1 + t2 + c1 + c2)
}

BVNderiv.muhat= function(theta, lambda, eta) {
  h = 0.000001
  p1 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1]+h, lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  p2 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

BVNderiv.phihat= function(theta, lambda, eta) {
  h = 0.00001
  p1 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]+h), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  p2 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  out = (p1 - p2)/h
  return(out)
}


BVNderiv.lmu= function(theta, lambda, eta) {
  h = 0.000001
  p1 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2]+h, eta[1], eta[1], lambda[4]), 2), log = TRUE)
  p2 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

BVNderiv.lphi= function(theta, lambda, eta) {
  h = 0.000001
  p1 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]+h), 2), log = TRUE)
  p2 = dmvnorm(c(theta[1], theta[2]), mean = c(lambda[1], lambda[3]), sigma = matrix(c(lambda[2], eta[1], eta[1], lambda[4]), 2), log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

IGderiv.alpha = function(theta, lambda, eta) {
  h = 0.000001
  p1 = densigamma(theta[3], lambda[5] + h, lambda[6])
  p2 = densigamma(theta[3], lambda[5], lambda[6])
  out = (p1 - p2)/h
  return(out)
}

IGderiv.beta = function(theta, lambda, eta) {
  h = 0.000001
  p1 = densigamma(theta[3], lambda[5], lambda[6] + h)
  p2 = densigamma(theta[3], lambda[5], lambda[6])
  out = (p1 - p2)/h
  return(out)
}

Qnormderiv.mean = function(x, mean, var){
  h = 0.000001
  p1 = pnorm(x, mean+h, sqrt(var))
  p2 = pnorm(x, mean, sqrt(var))
  out = (p1 - p2)/h
  return(out)
}

Qnormderiv.var = function(x, mean, var){
  h = 0.000001
  p1 = pnorm(x, mean, sqrt(var+h))
  p2 = pnorm(x, mean, sqrt(var))
  out = (p1 - p2)/h
  return(out)
}

QIGderiv.alpha = function(theta, lambda, eta) {
  h = 0.000001
  p1 = pigamma(theta[3], lambda[5] + h, lambda[6])
  p2 = pigamma(theta[3], lambda[5], lambda[6])
  out = (p1 - p2)/h
  return(out)
}

QIGderiv.beta = function(theta, lambda, eta) {
  h = 0.000001
  p1 = pigamma(theta[3], lambda[5], lambda[6] + h)
  p2 = pigamma(theta[3], lambda[5], lambda[6])
  out = (p1 - p2)/h
  return(out)
}

copderiv.mu = function(theta, lambda, rho) {
  h=0.000001
  copula = normalCopula(rho)
  p1 = dCopula(c(pnorm(theta[1], mean = lambda[1], sd = sqrt(lambda[2])) + h, pigamma(theta[3], lambda[5], lambda[6])), copula, log = TRUE)
  p2 = dCopula(c(pnorm(theta[1], mean = lambda[1], sd = sqrt(lambda[2])), pigamma(theta[3], lambda[5], lambda[6])), copula, log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

copderiv.phi = function(theta, lambda, rho) {
  h=0.000001
  copula = normalCopula(rho)
  p1 = dCopula(c(pnorm(theta[2], mean = lambda[3], sd = sqrt(lambda[4])) + h, pigamma(theta[3], lambda[5], lambda[6])), copula, log = TRUE)
  p2 = dCopula(c(pnorm(theta[2], mean = lambda[3], sd = sqrt(lambda[4])), pigamma(theta[3], lambda[5], lambda[6])), copula, log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

copderiv.sigmu = function(theta, lambda, rho) {
  h=0.000001
  copula = normalCopula(rho)
  p1 = dCopula(c(pnorm(theta[1], mean = lambda[1], sd = sqrt(lambda[2])), pigamma(theta[3], lambda[5], lambda[6]) + h), copula, log = TRUE)
  p2 = dCopula(c(pnorm(theta[1], mean = lambda[1], sd = sqrt(lambda[2])), pigamma(theta[3], lambda[5], lambda[6])), copula, log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

copderiv.sigphi = function(theta, lambda, rho) {
  h=0.000001
  copula = normalCopula(rho)
  p1 = dCopula(c(pnorm(theta[2], mean = lambda[3], sd = sqrt(lambda[4])), pigamma(theta[3], lambda[5], lambda[6])  + h), copula, log = TRUE)
  p2 = dCopula(c(pnorm(theta[2], mean = lambda[3], sd = sqrt(lambda[4])), pigamma(theta[3], lambda[5], lambda[6])), copula, log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

muhat.deriv= function(theta, lambda, eta) {
  p1 = BVNderiv.muhat(theta, lambda, eta) + Qnormderiv.mean(theta[1], lambda[1], lambda[2])*copderiv.mu(theta, lambda, eta[2])
  return(p1)
}

lmu.deriv = function(theta, lambda, eta) {
  p1 = BVNderiv.lmu(theta, lambda, eta) + Qnormderiv.var(theta[1], lambda[1], lambda[2])*copderiv.mu(theta, lambda, eta[2])
  return(p1)
}

phihat.deriv= function(theta, lambda, eta) {
  p1 = BVNderiv.phihat(theta, lambda, eta) + Qnormderiv.mean(theta[2], lambda[3], lambda[4])*copderiv.phi(theta, lambda, eta[3])
  return(p1)
}

lphi.deriv = function(theta, lambda, eta) {
  p1 = BVNderiv.lphi(theta, lambda, eta) + Qnormderiv.var(theta[2], lambda[3], lambda[4])*copderiv.phi(theta, lambda, eta[3])
  return(p1)
}

alpha.deriv = function(theta, lambda, eta) {
  p1 = IGderiv.alpha(theta, lambda, eta) + QIGderiv.alpha(theta, lambda, eta) * 
    (copderiv.sigmu(theta, lambda, eta[2]) + copderiv.sigphi(theta, lambda, eta[3]))
  return(p1)
}

beta.deriv = function(theta, lambda, eta) {
  p1 = IGderiv.beta(theta, lambda, eta) + QIGderiv.beta(theta, lambda, eta) * 
    (copderiv.sigmu(theta, lambda, eta[2]) + copderiv.sigphi(theta, lambda, eta[3]))                                                                                 
  return(p1)
}

part2 = function(y, theta, lambda, eta) {
  logjoint(y, theta[1], theta[2], theta[3]) - logq(theta, lambda, eta)
}

#VB iteration storage
lambda = matrix(0, nrow = 100, ncol = 6)
eta = matrix(0, nrow = 100, ncol = 3)
i.lambda = c(2, 0.5, 0.5, 0.02, 10, 10)
i.eta = c(0, 0, 0)

p = 0.01/(1:10)
pl = 0.001/(1:10)
pt = cbind(p, pl, p, pl, p, p)

#First iteration: Draw simulations
s = 2000
theta = cop.sim(i.lambda, i.eta, s)

#First iteration: Calculate partial derivatives
muhat.partial = vector(length = s)
lmu.partial = vector(length = s)
phihat.partial = vector(length = s)
lphi.partial = vector(length = s)
alpha.partial = vector(length = s)
beta.partial = vector(length = s)
for(j in 1:s) {
  logdens = part2(y, theta[j,], i.lambda, i.eta)
  muhat.partial[j] = logdens * muhat.deriv(theta[j,], i.lambda, i.eta)
  lmu.partial[j] = logdens * lmu.deriv(theta[j,], i.lambda, i.eta)
  phihat.partial[j] = logdens * phihat.deriv(theta[j,], i.lambda, i.eta)
  lphi.partial[j] = logdens * lphi.deriv(theta[j,], i.lambda, i.eta)
  alpha.partial[j] = logdens * alpha.deriv(theta[j,], i.lambda, i.eta)
  beta.partial[j] = logdens * beta.deriv(theta[j,], i.lambda, i.eta)
}
mustep <- mean(muhat.partial)
lmustep <- mean(lmu.partial)qq
phistep <- mean(phihat.partial)
lphistep <- mean(lphi.partial)
alphastep <- mean(alpha.partial)
betastep <- mean(beta.partial)

steps = c(mustep, lmustep, phistep, lphistep, alphastep, betastep)
lambda[1,] = i.lambda + pt[1,]*steps
