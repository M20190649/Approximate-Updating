library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
sourceCpp('ADPF.cpp')
sourceCpp('PMMH.cpp')

sigSq = 0.02
phi = 0.8
mu = -0.5
T = 1000

x0 = rnorm(1, mu, sqrt(sigSq/(1-phi^2)))
x = rep(0, T)
y = rep(0, T)
for(t in 1:T){
  if(t == 1){
    x[t] = mu + phi*(x0-mu) + rnorm(1, 0, sqrt(sigSq))
  } else {
    x[t] = mu + phi*(x[t-1]-mu) + rnorm(1, 0, sqrt(sigSq))
  }
  y[t] = rnorm(1, 0, exp(x[t]/2))
}
yM = matrix(y, T)


MCMC = PMMH(y, 50000, 10000, 200, c(2.5, 0.025, 0, 10, 20, 1.5), c(0.05, 0.05, 0.05))

meam = c(-3, 0, 0.5, 0)
var = diag(c(1, 0.3, 0.05, 0.3))
lambda = cbind(mean, var)
VBfit = VBIL_PF(yM, lambda, S=3, alpha=0.05, maxIter=5000, threshold=0.01, thresholdIS=0)

U = VBfit$U
U[is.na(U)] = 0
Sigma = t(U) %*% U
supPos = seq(0.001, 0.1, length.out=500)
supR = seq(-3, 3, length.out=500)
qplot(supPos, dlnorm(supPos, VBfit$Mu[1], sqrt(Sigma[1, 1])), geom='line')
qplot(supR, dnorm(supR, VBfit$Mu[2], sqrt(Sigma[2, 2])), geom='line')
qplot(supR, dnorm(supR, VBfit$Mu[3], sqrt(Sigma[3, 3])), geom='line')
qplot(supR, dnorm(supR, VBfit$Mu[4], sqrt(Sigma[4, 4])), geom='line')
