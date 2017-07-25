library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
sourceCpp('PMMH_toy.cpp')
sourceCpp('ADToy.cpp')

sigSqV = 10
sigSqW = 10
T = 100

x = rep(0, T)
y = rep(0, T)

for(t in 1:100){
  if(t==1){
    x[t] = rnorm(1, 0, sqrt(5))
  } else {
    x[t] = x[t-1]/2 + 25*x[t-1] / (1+x[t-1]^2) + 8*cos(1.2*t) + rnorm(1, 0, sqrt(sigSqV))
  }
  y[t] = x[t]^2 / 20 + rnorm(1, 0, sqrt(sigSqW))
}

MCMC = PMMHToy(y, 20000, 0, 50, c(1, 1, 1, 1), c(0.01, 0.01))

mean = c(2.3, 2.3)
var = c(0.1, 0, 0.1)
lambda = c(mean, var)
yM = matrix(y, T)
lM = matrix(lambda, 5)
VBfit = VBIL_Toy(yM, lM, S=3, alpha=0.05, maxIter=0, threshold=0.01, thresholdIS=0)

U = VBfit$U
U[is.na(U)] = 0
Sigma = t(U) %*% U
supPos = seq(0.001, 15, length.out=500)
qplot(supPos, dlnorm(supPos, VBfit$Mu[1], sqrt(Sigma[1, 1])), geom='line')
qplot(supPos, dlnorm(supPos, VBfit$Mu[2], sqrt(Sigma[2, 2])), geom='line')
