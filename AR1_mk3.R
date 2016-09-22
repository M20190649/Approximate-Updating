library(GGally)
library(coda)
library(ggplot2)
library(gridExtra)
library(reshape)

T = 100    #Length of y
N = 5000   #MCMC draws

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
psi = 10
sigbar = 0
gamma = 10
p1 = 1
p2 = 1
a = 1
b = 1


x = seq(-0.99, 0.99, 0.03)
sbeta <- function(x, p1, p2) {
  ((1+x)/2)^(p1-1)*((1-x)/2)^(p2-1)
}

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
  denom = psi^2 * ((T-1)*(1-theta[(i-1),2])^2+(1-theta[(i-1),2]^2)) + theta[(i-1), 3]
  sigsqmu = theta[(i-1), 3]*psi^2 / denom
  muhat = (psi^2*((1-theta[(i-1),2]^2)*y[1] + (1-theta[(i-1),2])*sum(y[2:T]-theta[(i-1),2]*y[1:(T-1)])) + mubar*theta[(i-1), 3])/denom
  theta[i,1] = rnorm(1, muhat, sqrt(sigsqmu))
  
  #phi conditional
  #phihat = sum((y[2:T]-theta[i, 1])*(y[1:(T-1)]-theta[i,1]))/sum((y[1:(T-1)]-theta[i,1])^2)   #from the Catherine's notes for ETC4541, pg 84
  #vphi = theta[(i-1),3]/sum((y[1:(T-1)]-theta[i,1])^2)    #acceptance rate for this candidate ~75%, trying random walk instead
  #phidraw = rnorm(1, phihat, sqrt(vphi))
  phidraw = rnorm(1, theta[(i-1),2], 0.18)               #random walk MH candidate draw to tune acceptance ratios to ~25%
  
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
thin[,3] = log(thin[,3])
thin = data.frame(thin)
MCMCmu = ggplot(data=thin, aes(x=Mu)) + geom_density() + labs(y=NULL)
MCMCphi = ggplot(data=thin, aes(x=Phi)) + geom_density() + labs(y=NULL)
MCMCsig = ggplot(data=thin, aes(x=Sigma2)) + geom_density() + labs(x = "Log Sigma Squared", y=NULL)
MCMCmp = ggplot(data=thin, aes(x=Mu, y = Phi)) + geom_point()
MCMCms = ggplot(data=thin, aes(x=Mu, y = Sigma2)) + geom_point() + labs(y = "Log Sigma Squared")
MCMCps = ggplot(data=thin, aes(x=Phi, y = Sigma2)) + geom_point() + labs(y = "Log Sigma Squared")

#Bayesian Essentials MCMC#

p = 1
Psi = matrix(0, ncol = p, nrow = p+1)
Psi[1,1] = 1
Psi[2,1] = -0.5
Rho = Psi[1:2, 1]

loglike = function(x, mu, Rho, sig2) {
  T = length(x)
  x = x - mu
  loglike = 0
  for(i in 2:T) {
    loglike = loglike - (t(Rho)%*%x[i:(i-1)])^2
  }
  loglike = (loglike/(2*sig2)-(T/2)*log(sig2))
  return(loglike)
}

theta.be = matrix(0, nrow = N, ncol = 3)
theta.be[1, ] = c(2, 0.2, 1)
accept = c(0, 0, 0)

for(i in 2:N){
  phi.candidate = runif(1, 0, 1)
  alpha = min(1, exp(loglike(y, theta.be[(i-1), 1], c(1, -phi.candidate), theta.be[(i-1), 3]^2)) /
                exp(loglike(y, theta.be[(i-1), 1], c(1, -theta.be[(i-1), 2]), theta.be[(i-1), 3]^2)))
  u = runif(1)
  if(u < alpha) {
    theta.be[i, 2] = phi.candidate
    accept[2] = accept[2] + 1
  } else {
    theta.be[i, 2] = theta.be[(i-1), 2]
  }
  mu.candidate = rnorm(1, theta.be[(i-1), 1], 0.5)
  alpha = min(1, exp(loglike(y, mu.candidate, c(1, -theta.be[i, 2]), theta.be[(i-1), 3]^2)) /
                exp(loglike(y, theta.be[(i-1), 1],  c(1, -theta.be[i, 2]), theta.be[(i-1), 3]^2)))
  u = runif(1)
  if(u < alpha) {
    theta.be[i, 1] = mu.candidate
    accept[1] = accept[1] + 1
  } else {
    theta.be[i, 1] = theta.be[(i-1), 1]
  }
  
  sig.candidate = runif(1, 0.5, 2.5)
  alpha = min(1, exp(loglike(y, theta.be[i, 1],  c(1, -theta.be[i, 2]), sig.candidate^2)) /
                exp(loglike(y, theta.be[i, 1],  c(1, -theta.be[i, 2]), theta.be[(i-1), 3]^2)))
  u = runif(1)
  if(u < alpha) {
    theta.be[i, 3] = sig.candidate
    accept[3] = accept[3] + 1
  } else {
    theta.be[i, 3] = theta.be[(i-1), 3]
  }
}
accept
effectiveSize(theta.be)
theta.be = theta.be[seq(1001, 5000, 10),]
effectiveSize(theta.be)
theta.be[,3] = theta.be[,3]^2
colnames(theta.be) = c("Mu", "Phi", "Sigma2")
ggpairs(theta.be)
ggpairs(thin)


### Black Box VI
library(mvtnorm)

#Starting components

#no longer needed
simul = function(n, lambda, L) {
  out = matrix(-5, ncol = 3, nrow = n)
  for(i in 1:n){
    test = c(-5, -5, -5)
    while(test[2] > 1 | test[2] < -1) {
      out[i, ] = rmvnorm(1, rep(0, 3), diag(1, 3))
      test = L %*% out[i, ] + lambda
    }
  }
  return(out)
}

logjoint = function(y, mu, phi, sig2) {
  t1 = 1/2*log(1-phi^2) -(T/2)*sig2 - (1-phi^2)*(y[1]-mu)^2/(2*exp(sig2))
  t2 = - sum((y[2:T] - mu - phi*(y[1:(T-1)] - mu))^2)/(2*exp(sig2))
  t3 = -(mu-mubar)^2/(2*psi^2) - (sig2 - sigbar)^2/gamma
  t4 = (p1-1)*log((1+phi)/2) + (p2-1)*log((1-phi)/2)
  return(t1 + t2 + t3 + t4)
}

logq = function(theta, lambda, Sigma) {
  dmvnorm(theta, lambda, Sigma, log = TRUE)
}

#no longer needed
MVNderiv = function(theta, lambda, Sigma, param) {
  h = 0.000001
  theta2 = theta
  lambda2 = lambda
  Sigma2 = Sigma
  if(param == 1 | param == 2 | param == 3) {
    theta2[param] = theta[param] + h
  } else if(param == 4 | param == 5 | param ==6) {
    Sigma2[(param - 3), (param - 3)] = Sigma[(param - 3), (param - 3)] + h
  } else if(param == 7) {
    Sigma2[1, 2] = Sigma2[2, 1] = Sigma[1, 2] + h
  } else if(param == 8) {
    Sigma2[1, 3] = Sigma2[3, 1] = Sigma[1, 3] + h
  } else {
    Sigma2[2, 3] = Sigma2[3, 2] = Sigma[2, 3] + h
  }
  p1 = dmvnorm(theta2, lambda2, Sigma2, log = TRUE)
  p2 = dmvnorm(theta, lambda, Sigma, log = TRUE)
  out = (p1 - p2)/h
  return(out)
}

jointderiv = function(theta, L, lambda, y, param) {
  h =  0.000001
  L2 = L
  lambda2 = lambda
  if(param %in% c(1, 2, 3)) {
    lambda2[param] = lambda[param] + h
  } else if(param %in% c(4, 5, 6)) {
    L2[(param - 3), (param - 3)] = L[(param - 3), (param - 3)] + h
  } else if(param == 7) {
    L2[2, 1] = L[2, 1] + h
  } else if(param == 8) {
    L2[3, 1] = L[3, 1] + h
  } else {
    L2[3, 2] = L[3, 2] + h
  }
  params = L %*% theta + lambda
  params2 = L2 %*% theta + lambda2
  p1 = logjoint(y, params[1], params[2], params[3])
  p2 = logjoint(y, params2[1], params2[2], params2[3])
  return((p2-p1)/h)
}

jakobderiv = function(theta, L, lambda, param) {
  if(param %in% c(1, 2, 4, 5, 7, 8, 9)) {
    return(0)
  } else if(param == 3) {
    return(1)
  } else {
    return(theta[3])
  }
}

#no longer needed
part2 = function(y, theta, lambda, Sigma) {
  logjoint(y, theta[1], theta[2], theta[3]) - logq(theta, lambda, Sigma)
}

#number of iterations
M = 1000
#VB iteration storage
lambda = matrix(0, nrow = M, ncol = 3)
#variances = matrix(0, nrow = 100, ncol = 3)
#covariances = matrix(0, nrow = 100, ncol = 3)
i.lambda = c(mean(y), cor(y[2:T], y[1:(T-1)]), log(var(y)))
#i.Sigma = diag(c(0.15, 0.15, 0.15))
Lp = matrix(0, nrow = M, ncol = 6)
i.L = diag(0.4, 3)

p = 0.01/(1:M)
#pl = 0.001/(1:M)

#First iteration: Draw simulations
s = 10
theta = simul(s, i.lambda, i.L)
#First iteration: Calculate partial derivatives
partials = matrix(0, nrow = s, ncol = 9)

for(j in 1:s) {
  #logdens = part2(y, theta[j,], i.lambda, i.Sigma)
  
  #partial derivs
  for(i in c(1:3, 7:9)) {
    partials[j, i] = jointderiv(theta[j, ], i.L, i.lambda, y, i) + jakobderiv(theta[j, ], i.L, i.lambda, i)
  }
  for(i in 4:6) {
    partials[j, i] = jointderiv(theta[j, ], i.L, i.lambda, y, i) + jakobderiv(theta[j, ], i.L, i.lambda, i) + 1/(i.L[(i - 3), (i - 3)])
  }
}

steps = colMeans(partials, na.rm = TRUE)
lambda[1,] = i.lambda + p[1]*steps[1:3]
Lp[1, ] = c(0.4, 0.4, 0.4, 0, 0, 0)+ p[1]*steps[4:9]


#Repeat iterations
for(k in 2:M) {
  
  L = diag(Lp[(k-1), 1:3])
  L[lower.tri(L)] = Lp[(k-1), 4:6]
  theta = simul(s, lambda[(k-1),], L)
  
  partials = matrix(0, nrow = s, ncol = 9)
  
  for(j in 1:s) {
    #logdens = part2(y, theta[j,], i.lambda, i.Sigma)
    
    #partial derivs
    for(i in c(1:3, 7:9)) {
      partials[j, i] = jointderiv(theta[j, ], L, lambda[(k-1), ], y, i) + jakobderiv(theta[j, ], L, lambda[(k-1), ], i)
    }
    for(i in 4:6) {
      partials[j, i] = jointderiv(theta[j, ], L, lambda[(k-1), ], y, i) + jakobderiv(theta[j, ], L, lambda[(k-1), ], i) + 1/(L[(i - 3), (i - 3)])
    }
  }
  steps = colMeans(partials, na.rm = TRUE)
  lambda[k,] = lambda[(k-1), ] + p[k]*steps[1:3]
  Lp[k, ] = Lp[(k-1), ] + p[k]*steps[4:9]
}

Lf = diag(Lp[M, 1:3])
Lf[lower.tri(Lf)] = Lp[M, 4:6]
Sigma = Lf %*% t(Lf)

#noise in the draws
VBdraws = rmvnorm(2000, lambda[M,], Sigma)
VBdraws[,3] = exp(VBdraws[,3])
ggpairs(VBdraws, lower = list(continuous = "density"))

#univariate densities
xm = seq(0.5, 4, length.out=1000)
mdens = dnorm(xm, lambda[M, 1], sqrt(Sigma[1, 1]))
mudat = data.frame(cbind(xm, mdens))
muplot = ggplot(mudat, aes(x=xm, y=mdens)) + geom_line() + labs(x="Mu", y=NULL)

xp = seq(0.2, 1, length.out=1000)
pdens = dnorm(xp, lambda[M, 2], sqrt(Sigma[2, 2]))
phidat = data.frame(cbind(xp, pdens))
phiplot = ggplot(phidat, aes(x=xp, y=pdens)) + geom_line() + labs(x="Phi", y=NULL)

xs = seq(0, 2, length.out = 1000)
sdens = dnorm(xs, lambda[M, 3], sqrt(Sigma[3, 3]))
sigdat = data.frame(cbind(xs, sdens))
sigplot = ggplot(sigdat, aes(x=xs, y=sdens)) + geom_line() + labs(x="Log Sigma Squared", y=NULL)

#bivariate densities

x1 = seq(0, 3, length.out=1000)
x2 = seq(0, 3, length.out=1000)
z = matrix(0, length(x1), length(x2))
for (i in 1:length(x1)) {
  a = x1
  b = x2[i]
  z[,i] = dmvnorm(cbind(a,b), lambda[M, 1:2], Sigma[1:2, 1:2])
}
colnames(z) = x1
rownames(z) = x2
# reshape the data
dat.mp <- melt(z)
colnames(dat.mp) = c("Mu", "Phi", "value")

muphi <- ggplot(dat.mp, aes(x=Mu, y=Phi, z = value)) + geom_contour()

for (i in 1:length(x1)) {
  a = x1
  b = x2[i]
  z[,i] = dmvnorm(cbind(a,b), lambda[M, -2], Sigma[-2, -2])
}
colnames(z) = x1
rownames(z) = x2
# reshape the data
dat.ms <- melt(z)
colnames(dat.ms) = c("Mu", "Ln_Sigma2", "value")

musig <- ggplot(dat.ms, aes(x=Mu, y=Ln_Sigma2, z = value)) + geom_contour() + labs(y = "Log Sigma Squared")

for (i in 1:length(x1)) {
  a = x1
  b = x2[i]
  z[,i] = dmvnorm(cbind(a,b), lambda[M, -1], Sigma[-1, -1])
}
colnames(z) = x1
rownames(z) = x2
# reshape the data
dat.ps <- melt(z)
colnames(dat.ps) = c("Phi", "Ln_Sigma2", "value")
phisig <- ggplot(dat.ps, aes(x=Phi, y=Ln_Sigma2, z = value)) + geom_contour() + labs(y = "Log Sigma Squared")

grid.arrange(muplot, phiplot, sigplot, muphi, musig, phisig, ncol = 3)
grid.arrange(MCMCmu, MCMCphi, MCMCsig, MCMCmp, MCMCms, MCMCps, ncol = 3)

