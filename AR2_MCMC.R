library(Rcpp)
library(ggplot2)
library(dplyr)
library(RcppArmadillo)
library(GGally)
library(coda)
library(gridExtra)
library(reshape)
library(mvtnorm)
library(pscl)
library(microbenchmark)
library(copula)

sourceCpp("AR2.cpp")


T = 100    #Length of y
N = 120000   #MCMC draws

#True parameters
sigma2 <- 2
phi1 <- 0.7
phi2 <- 0.15

phi1b = 0
phi2b = 0
phi1l = 10
phi2l = 10

a = 1
b = 1
shape = T/2 -1 + a


#Generate data
y = c(0, 0.5)
for(i in 3:(T+50)){
  y = c(y, phi1*y[(i-1)] + phi2*y[(i-2)] + rnorm(1, 0, sqrt(sigma2)))
}
y = y[51:(T+50)]

phistar = function(y, phi1, phi2, sig2){
  loglike = -(T-2)/2*log(sig2) - sum((y[3:T] - phi1*y[2:(T-1)] - phi2*y[1:(T-2)])^2)/(2*sig2)
  logprior = -(phi1-phi1b)^2/(2*phi1l) -(phi2-phi2b)^2/(2*phi2l)
  return(loglike + logprior)
}

ym22 = sum(y[1:(T-2)]^2)
y2 = sum(y[3:T]^2)
ym12 = sum(y[2:(T-1)]^2)
yym1 = sum(y[3:T]*y[2:(T-1)])
yym2 = sum(y[3:T]*y[1:(T-2)])
ym1ym2 = sum(y[2:(T-1)]*y[1:(T-2)])

theta = matrix(0, ncol = 3, nrow = N)
theta[1, 3] = 1
theta[1, 1:2] = 0.3
for(i in 2:N){
  #phi 1 conditional
  mean1 = (phi1l*(yym1 - theta[i-1, 2]*ym1ym2) + theta[i-1, 3]*phi1b)/(phi1l*ym12 + theta[i-1, 3])
  var1 = theta[i-1, 3]*phi1l/(phi1l*ym12+theta[i-1, 3]) 
  theta[i, 1] = rnorm(1, mean1, sqrt(var1))
  
  #phi2 conditional
  mean2 = (phi2l*(yym2 - theta[i, 1]*ym1ym2) + theta[i-1, 3]*phi2b)/(phi2l*ym22 + theta[i-1, 3])
  var2 = theta[i-1, 3]*phi2l/(phi2l*ym22+theta[i-1, 3]) 
  theta[i, 2] = rnorm(1, mean2, sqrt(var2))
  
  #sigma2 conditional
  inv.scale = b + sum(((y[3:T] - theta[i, 1]*y[2:(T-1)]- theta[i, 2]*y[1:(T-2)])^2))/2
  theta[i, 3] = 1/rgamma(1, shape = shape, rate = inv.scale)
}


keep = theta[(N/5+1):N, ]
effectiveSize(keep)
thin = keep[seq(1, 4*N/5, 20),]
effectiveSize(thin)
colnames(thin) = c("Phi1", "Phi2", "SigmaSquared")

ggpairs(thin)

#Fit copulas
#Normal copula is the best by BIC
psuedo = pobs(thin)
normcop = normalCopula(dim = 2)
fitn = fitCopula(normcop, psuedo[,1:2])
BICn = -2*fitn@loglik+log(4000)
tcop = tCopula(dim = 2)
fitt = fitCopula(tcop, psuedo[,1:2])
BICt = -2*fitt@loglik+2*log(4000)

#BVN Normal + IG
#Phi 1 - Phi 2 BVN Normal(mu_1, mu_2, taul, taul, rho)
#convert phi1l, phi2l, rho to L
#Sigma2 Inv.Gamma(alpha, beta)

phi1b = 0
phi2b = 0
phi1l = 10
phi2l = 10

a = 1
b = 1

sim2 = function(n, lambda) {
  L = matrix(c(lambda[3:4], 0, lambda[5]), nrow = 2)
  draws = rmvnorm(n, c(0,0))
  phi = t(apply(draws, 1, function(x) lambda[1:2] + L %*% x))
  sig2 = 1/rgamma(n, exp(lambda[6]), exp(lambda[7]))
  out = cbind(phi, sig2)
  colnames(out) = c("Phi1", "Phi2", "Sig2")
  return(out)
}

logj2 = function(y, phi1, phi2, sig2) {
  loglike = -(T-2)/2*log(sig2) - sum((y[3:T] - phi1*y[2:(T-1)] - phi2*y[1:(T-2)])^2)/(2*sig2)
  z = (phi1-phi1b)^2/phi1l + (phi2-phi2b)^2/phi2l
  p1 = -z/2
  p2 = -(a+1)*log(sig2) - b/sig2
  return(loglike + p1 + p2)
}

logq2 = function(theta, lambda) {
  L = matrix(c(lambda[3:4], 0, lambda[5]), nrow = 2)
  Sigma = L %*% t(L)
  bvn = dmvnorm(theta[1:2], lambda[1:2], Sigma, log = TRUE)
  ig = log(densigamma(theta[3], exp(lambda[6]), exp(lambda[7])))
  return(bvn + ig)
}

logqderiv = function(theta, lambda, param){
  h = 0.000001
  lambda2 = lambda
  lambda2[param] = lambda2[param] + h
  p1 = logq2(theta, lambda2)
  p2 = logq2(theta, lambda)
  if(param < 6) {
   return((p1-p2)/h)
  } else {
    return(exp(lambda[param])*((p1-p2)/h))
  }
}

ELBO = function(lambda,n=10){
  L = matrix(c(lambda[3:4], 0, lambda[5]), 2)
  Sigma = L %*% t(L)
  theta = cbind(rmvnorm(n, c(lambda[1], lambda[2]), Sigma), 1/rgamma(n, exp(lambda[6]), exp(lambda[7])))
  out = apply(theta, 1, function(x) logj2(y, x[1], x[2], x[3]) - logq2(x, lambda))
  return(mean(out))
}

M = 2000
lambda2 = matrix(0, ncol = 7, nrow = M)
Sigma = cov(thin)[1:2, 1:2]
i.L = t(chol(Sigma))

mu = mean(thin[,3])
v = var(thin[,3])

i.a = mu^2/v + 2
i.b = mu*(mu^2/v + 1)

i.lambda = c(mean(thin[,1]), mean(thin[,2]), i.L[1,1], i.L[2,1], i.L[2,2], log(i.a), log(i.b))
s = 50
LB = vector(length = M)
var = vector(length = M)
#first iteration
theta = sim2(s, i.lambda)
partials2 = matrix(0, ncol = 7, nrow = s)
h = matrix(0, ncol = 7, nrow = s)
for(i in 1:s){
  part2 = logj2(y, theta[i,1], theta[i, 2], theta[i, 3]) - logq2(theta[i,], i.lambda)
  for(j in 1:7){
    h[i, j] = logqderiv(theta[i,], i.lambda, j)
    partials2[i, j] = h[i, j] * part2
  }
}
rows = sample(1:s, 20)
astar = vector(length = 7)
for(j in 1:7){
  astar[j] = cov(partials2[rows, j], h[rows, j])/var(h[rows, j])
}
var[1] = var(partials2[,7])
fd = colMeans(partials2)
hd = colMeans(h)
steps2 = fd - astar*hd
G_1 = steps2 %*% t(steps2)
Gt = G_1
p_1 = 0.1* diag(G_1^(-1/2))
lambda2[1, ] = i.lambda + p_1*steps2
L.i = ELBO(i.lambda, 100)
LB[1] = ELBO(lambda2[1, ], 100)
diff = LB[1] - L.i

#loop
k = 2

while(abs(diff) > 0.001) {
  if(k > M){
    break
  }
  theta = sim2(s, lambda2[k-1,])
  partials2 = matrix(0, ncol = 7, nrow = s)
  h = matrix(0, ncol = 7, nrow = s)
  for(i in 1:s){
    part2 = logj2(y, theta[i,1], theta[i, 2], theta[i, 3]) - logq2(theta[i,], lambda2[k-1,])
    for(j in 1:7){
      h[i, j] = logqderiv(theta[i,], lambda2[k-1,], j)
      partials2[i, j] = h[i, j] * part2
    }
  }
  rows = sample(1:s, 20)
  astar = vector(length = 7)
  for(j in 1:7){
    astar[j] = cov(partials2[rows, j], h[rows, j])/var(h[rows, j])
  }
  var[k] = var(partials2[,7])
  fd = colMeans(partials2)
  hd = colMeans(h)
  steps2 = fd- astar*hd
  
  Gt = Gt + steps2 %*% t(steps2)
  p_t = 0.1 * diag(Gt^(-1/2))
  lambda2[k, ] = lambda2[k-1,] + p_t*steps2
  LB[k] = ELBO(lambda2[k,], 100)
  diff = LB[k] - LB[k-1]
  k = k + 1
}


L = matrix(c(lambda2[k-1, 3:4], 0, lambda2[k-1, 5]), 2)
Sigma = L %*% t(L)

xp1 = seq(0.3, 1.1, length.out=1000)
p1dens = dnorm(xp1, lambda2[k-1, 1], sqrt(Sigma[1, 1]))
p1dat = data.frame(cbind(xp1, p1dens))
p1plot = ggplot(p1dat, aes(x=xp1, y=p1dens)) + geom_line() + labs(x="Phi 1", y=NULL)

xp2 = seq(-0.3, 0.7, length.out=1000)
p2dens = dnorm(xp2, lambda2[k-1, 2], sqrt(Sigma[2, 2]))
p2dat = data.frame(cbind(xp2, p2dens))
p2plot = ggplot(p2dat, aes(x=xp2, y=p2dens)) + geom_line() + labs(x="Phi 2", y=NULL)

xs = seq(0.5, 3.5, length.out = 1000)
sdens = densigamma(xs, exp(lambda2[k-1, 6]), exp(lambda2[k-1, 7]))
sigdat = data.frame(cbind(xs, sdens))
sigplot = ggplot(sigdat, aes(x=xs, y=sdens)) + geom_line() + labs(x="Sigma Squared", y=NULL)

x1 = seq(-1, 1, length.out=1000)
x2 = seq(-1, 1, length.out=1000)
z = matrix(0, length(x1), length(x2))
for (i in 1:length(x1)) {
  d1 = x1
  d2= x2[i]
  z[,i] = dmvnorm(cbind(d1,d2), lambda2[k-1, 1:2], Sigma)
}
colnames(z) = x1
rownames(z) = x2
# reshape the data
dat.p1p2 <- melt(z)
colnames(dat.p1p2) = c("Phi1", "Phi2", "value")

p1p2 <- ggplot(dat.p1p2, aes(x=Phi1, y=Phi2, z = value)) + geom_contour()+ coord_equal()
grid.arrange(p1plot, p2plot, sigplot, p1p2)

thin = data.frame(thin)
MCMCp1 = ggplot(data=thin, aes(x=Phi1)) + geom_density() + labs(y=NULL) 
MCMCp2 = ggplot(data=thin, aes(x=Phi2)) + geom_density() + labs(y=NULL)
MCMCsig = ggplot(data=thin, aes(x=SigmaSquared)) + geom_density() + labs(y=NULL)
MCMCpp = ggplot(data=thin, aes(x=Phi1, y = Phi2)) + geom_point()+ coord_equal()
grid.arrange(MCMCp1, MCMCp2, MCMCsig, MCMCpp)
#does not capture var matrix well

#Mean Field
phi1b = 0
phi2b = 0
phi1l = 10
phi2l = 10
a = 1
b = 1
astar = a + T/2 -1
mf = matrix(0, ncol = 5, nrow = 10)
colnames(mf) = c("bstar", "phi1h", "l1", "phi2h", "l2")
mf[1,] = c(1, 0.4, 0.1, -0.2, 0.1)

#sum stats
ym22 = sum(y[1:(T-2)]^2)
y2 = sum(y[3:T]^2)
ym12 = sum(y[2:(T-1)]^2)
yym1 = sum(y[3:T]*y[2:(T-1)])
yym2 = sum(y[3:T]*y[1:(T-2)])
ym1ym2 = sum(y[2:(T-1)]*y[1:(T-2)])

for(i in 2:10){
  #beta
  mf[i, 1] = b + 1/2*(y2 + ym12*(mf[i-1, 3] + mf[i-1, 2]^2) + ym22*(mf[i-1, 5] + mf[i-1, 4]^2)) + 
    ym1ym2*mf[i-1, 2]*mf[i-1, 4] - mf[i-1, 2]*yym1 - mf[i-1, 4]*yym2
  
  #phi1
  mf[i, 3] = 1/(ym12*astar/mf[i, 1] + 1/phi1l)
  mf[i, 2] = (phi1l*(yym1 - mf[i-1, 4]*ym1ym2) + mf[i, 1]/(astar -1)*phi1b)/(ym12*phi1l + mf[i, 1]/(astar-1))
  
  #phi2
  mf[i, 5] = 1/(ym22*astar/mf[i, 1] + 1/phi2l)
  mf[i, 4] = (phi2l*(yym2 - mf[i, 2]*ym1ym2) + mf[i, 1]/(astar -1)*phi2b)/(ym22*phi2l + mf[i, 1]/(astar-1))
}
ELBOc(c(mf[10, 2], mf[10, 4], chol(diag(c(mf[10, 3], mf[10, 5])))[1,1], 0, chol(diag(c(mf[10, 3], mf[10, 5])))[2,2], astar, mf[10, 1]), 1000, y)

xp1 = seq(-0.1, 1, length.out=1000)
p1dens = dnorm(xp1, mf[10, 2], sqrt(mf[10, 3]))
p1dat = data.frame(cbind(xp1, p1dens))
p1plot = ggplot(p1dat, aes(x=xp1, y=p1dens)) + geom_line() + labs(x="Phi 1", y=NULL)

xp2 = seq(-0.8, 0.4, length.out=1000)
p2dens = dnorm(xp2, mf[10, 4], sqrt(mf[10, 5]))
p2dat = data.frame(cbind(xp2, p2dens))
p2plot = ggplot(p2dat, aes(x=xp2, y=p2dens)) + geom_line() + labs(x="Phi 2", y=NULL)

xs = seq(0.5, 3.5, length.out = 1000)
sdens = densigamma(xs, astar, mf[10, 1])
sigdat = data.frame(cbind(xs, sdens))
sigplot = ggplot(sigdat, aes(x=xs, y=sdens)) + geom_line() + labs(x="Sigma Squared", y=NULL)

grid.arrange(p1plot, p2plot, sigplot, ncol = 2)

#Inverse Transform VB Copula

sim3 = function(n) {
  out = rmvnorm(n, c(0,0))
  out = cbind(out, runif(n))
  return(out)
}

logj2 = function(y, phi1, phi2, sig2) {
  loglike = -(T-2)/2*log(sig2) - sum((y[3:T] - phi1*y[2:(T-1)] - phi2*y[1:(T-2)])^2)/(2*sig2)
  z = (phi1-phi1b)^2/phi1l + (phi2-phi2b)^2/phi2l
  p1 = -z/2
  p2 = -(a+1)*log(sig2) - b/sig2
  return(loglike + p1 + p2)
}

logq3 = function(theta, lambda) {
  L = matrix(c(lambda[3:4], 0, lambda[5]), nrow = 2)
  Sigma = L %*% t(L)
  bvn = dmvnorm(theta[1:2], lambda[1:2], Sigma, log = TRUE)
  ig = log(densigamma(theta[3], lambda[6], lambda[7]))
  return(bvn + ig)
}

allderiv = function(theta, lambda, y, param) {
  h =  0.000001
  lambda2 = lambda
  lambda2[param] = lambda[param] + h
  L = matrix(c(lambda[3], lambda[4], 0, lambda[5]), 2)
  L2 = matrix(c(lambda2[3], lambda2[4], 0, lambda2[5]), 2)
  
  params = c(lambda[1:2] + L %*% theta[1:2], qigamma(theta[3], lambda[6], lambda[7]))
  params2 =  c(lambda2[1:2] + L2 %*% theta[1:2], qigamma(theta[3], lambda2[6], lambda2[7]))
  
  p1 = logj2(y, params[1], params[2], params[3])
  p2 = logj2(y, params2[1], params2[2], params2[3])
  return((p2-p1)/h)
}

ELBO2 = function(lambda,n=10){
  L = matrix(c(lambda[3:4], 0, lambda[5]), 2)
  Sigma = L %*% t(L)
  theta = cbind(rmvnorm(n, c(lambda[1], lambda[2]), Sigma), 1/rgamma(n, lambda[6], lambda[7]))
  out = apply(theta, 1, function(x) logj2(y, x[1], x[2], x[3]) - logq3(x, lambda))
  return(mean(out))
}

M = 10000
lambda3 = matrix(0, ncol = 7, nrow = M)
Sigma = cov(thin)[1:2, 1:2]
i.L = t(chol(Sigma))

mu = mean(thin[,3])
v = var(thin[,3])

i.a = mu^2/v + 2
i.b = mu*(mu^2/v + 1)

LB2 = vector(length=M)
var = vector(length = M)
i.lambda = c(mean(thin[,1]), mean(thin[,2]), i.L[1,1], i.L[2,1], i.L[2,2], i.a, i.b)
#i.lambda = c(0,0,0.1,0,0.1,3,3)
s = 5
#p = 0.01/(1:M)
#p2 = 0.0001/(1:M)
#first iteration
theta = sim3(s)
partials3 = matrix(0, ncol = 7, nrow = s)
for(i in 1:s){
  for(j in 1:7){
    partials3[i, j] = allderiv(theta[i,], i.lambda, y, j)
  }
}
var[1] = var(partials3[,1])
steps3 = colMeans(partials3)
G_1 = steps3 %*% t(steps3)
Gt = G_1
p_1 = 0.1* diag(G_1^(-1/2))
lambda3[1, ] = i.lambda + p_1*steps3
L.i = ELBO2(i.lambda, 5)
LB2[1] = ELBO2(lambda3[1, ], 5)
diff = LB2[1] - L.i
#lambda3[1, 6:7] = i.lambda[6:7] + p2[1]*steps3[6:7]

#loop
k = 2
while(abs(diff) > 0.001) {
  if(k > M) {
    break
  }
  theta = sim3(s)
  for(i in 1:s){
    for(j in 1:7){
      partials3[i, j] = allderiv(theta[i,], lambda3[k-1,], y, j)
    }
  }
  var[k] = var(partials3[,1])
  steps3 = colMeans(partials3)
  Gt = Gt + steps3 %*% t(steps3)
  p_t = 0.1* diag(Gt^(-1/2))
  lambda3[k, ] = lambda3[k-1, ] + p_t*steps3
  LB2[k] = ELBO2(lambda3[k, ], 5)
  diff = LB2[k] - LB2[k-1]
  k = k+1
}

#qplot(1:395, lambda3[1:395, 1])

L = matrix(c(lambda3[k-1, 3:4], 0, lambda3[k-1, 5]), 2)
Sigma = L %*% t(L)
cov2cor(Sigma)
#write.csv(mf, "AR2mf.csv")
write.csv(thin, "AR2thin.csv")
#write.csv(lambda2, "AR2copula.csv")

der = function(x, theta, lambda, y, param) {
  lambda[7] = x
  mean(apply(theta, 1, allderiv, lambda, y = y , param = param))
}

test = sapply(seq(1, 100, 0.1), der, theta = theta, lambda = lambda2[500,], y = y, param = 7)
seq(1, 500, 1)[which(test < 0.5 & test > -0.5)]

a1 = seq(0.1, 150, length.out = 250)
b1 = seq(0.1, 150, length.out = 250)
z = matrix(0, length(a1), length(b1))
for(i in 1:length(a1)){
  for(j in 1:length(b1)){
    z[i,j] = ELBO(c(lambda3[M, 1:5], a1[i], b1[j]))
  }
}
z[which.max(z)]

gridsearch = function(phis, sigma2, N, T){
  phi1 = phis[1]
  phi2 = phis[2]
  phi1b = 0
  phi2b = 0
  phi1l = 10
  phi2l = 10
  a = 1
  b = 1

  y = c(0, 0.5)
  for(i in 3:(T+50)){
    y = c(y, phi1*y[(i-1)] + phi2*y[(i-2)] + rnorm(1, 0, sqrt(sigma2)))
  }
  y = y[51:(T+50)]
  
  draws = MCMC(N)
  return(cor(draws[,1], draws[,2]))
}
  
phi1 = seq(-1.99, 1.99, length.out = 20)
phi2 = seq(-0.99, 0.99, length.out = 20)
phigrid = expand.grid(phi1, phi2)
phigrid$c = apply(phigrid, 1, gridsearch, 2, 2000, 100)
ggplot(data=phigrid, aes(x=Var1, y = Var2, fill = c)) + geom_tile()