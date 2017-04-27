library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
sourceCpp("SVM_SGA.cpp")

T = 50
phi = 0.7
gamma = 1
sigmaSq = 1

a0 = rnorm(1, 0, sqrt(sigmaSq / (1-phi^2)))
alpha = rep(0, T)
y = rep(0, T)
for(t in 1:T){
  if(t == 1){
    alpha[t] = gamma + phi*a0 + sqrt(sigmaSq)*rnorm(1)
  } else {
    alpha[t] = gamma + phi*alpha[t-1] + sqrt(sigmaSq)*rnorm(1)
  }
  y[t] = alpha[t] + log(rnorm(1)^2)
}

fit = SVM_SGA(y, 10, 25, 5, rep(0, T+4), diag(0.1, T+4))
