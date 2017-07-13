library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
sourceCpp('TNK2016.cpp')

#forex = readxl::read_excel("2010-2013.xls", skip=10)
T = 1001
#y = rep(0, T)
#for(i in 1:T){
#  y[i] = 100*log(forex$FXRUSD[i+1]/forex$FXRUSD[i]- 1/T * sum(log(forex$FXRUSD[i+1]/forex$FXRUSD[i])))
#}

x0 = rnorm(1, 0.1, sqrt(0.1/(1-0.9^2)))
x = rep(0, T)
y = rep(0, T)
for(t in 1:T){
  if(t == 1){
    x[t] = 0.1 + 0.9 * (x[1]-0.1) + rnorm(1, 0, sqrt(0.1))
  } else {
    x[t] = 0.1 + 0.9 * (x[t-1]-0.1) + rnorm(1, 0, sqrt(0.1))
  }
  y[t] = exp(x[t]/2) * rnorm(1)
}

lambda = c(0, log(0.3), log(95), log(5), log(11), log(1))
N = 100
S = 256

supportMu = seq(-2, 2, length.out=500)
supportTau = seq(0.3, 0.9999, length.out=500)
supportSigSq = seq(0.0001, 0.3, length.out=500)

dMuSt = dnorm(supportMu, lambda[1], sqrt(exp(lambda[2])))
dTauSt = dbeta(supportTau, exp(lambda[3]), exp(lambda[4]))
dSigSt = pscl::densigamma(supportSigSq, exp(lambda[5]), exp(lambda[6])) 

df = data.frame(support = c(supportMu, supportTau, supportSigSq), density = c(dMuSt, dTauSt, dSigSt), 
                variable = rep(c('mu', 'tau', 'sigma^2'), rep(500, 3)), version='starting')

yM = matrix(y, T)
lM = matrix(lambda, 6)

VBfit = VBIL(yM, lM, 256, 100, maxIter=200, alpha=0.1)
lambda = c(VBfit$lambda)

dMuCon = dnorm(supportMu, lambda[1], sqrt(exp(lambda[2])))
dTauCon = dbeta(supportTau, exp(lambda[3]), exp(lambda[4]))
dSigCon = pscl::densigamma(supportSigSq, exp(lambda[5]), exp(lambda[6]))

df2 = data.frame(support = c(supportMu, supportTau, supportSigSq), density = c(dMuCon, dTauCon, dSigCon), 
                 variable = rep(c('mu', 'tau', 'sigma^2'), rep(500, 3)), version='converged')

df = rbind(df, df2)
ggplot(df) + geom_line(aes(x=support, y=density, colour=version)) + facet_wrap(~variable, labeller=label_parsed, scales='free')





# particle filter in R

phi = 0.8
mu = 0.5
sigSq = 0.2
T = 500
x = rep(T+1)
y = rep(T)
x[1] = rnorm(1, mu, sqrt(sigSq / (1-phi^2)))
for(t in 1:T){
  x[t+1] = rnorm(1, mu + phi*(x[t]-mu), sqrt(sigSq))
  y[t] = rnorm(1, 0, exp(x[t+1]))
}

N = 100
xOld = rnorm(N, mu, sqrt(sigSq / (1-phi^2)))
xNew = rep(0, N)
xRS = rep(0, N)
pi = rep(1/N, N)
pisum = pi
omega = rep(0, N)
yDens = 0
for(t in 1:T){
  for(k in 2:N){
    pisum[k] = pisum[k] + pisum[k-1]
  }
  for(k in 1:N){
    u = runif(1)
    i = 1
    flag = TRUE
    while(flag){
      if(u < pisum[i]){
        xRS[k] = xOld[i]
        flag = FALSE
      }
      if(u > pisum[N]){
        xRS[k] = xOld[N]
      }
      i = i + 1
    }
  }
  for(k in 1:N){
    xNew[k] = rnorm(1, mu + phi *(xRS[k] - mu), sqrt(sigSq))
  }
  for(k in 1:N){
    omega[k] = 1.0 / sqrt(2*3.141598) * exp(-xNew[k]/2) * exp(-y[t]^2 / (2*exp(xNew[k])))
  }
  pi = omega / sum(omega)
  yDens = yDens + log(sum(omega) / N)
}


