library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
sourceCpp("SVM_SGA.cpp")

set.seed(12)
T = 50
phi = 0.7
gamma = 1
sigmaSq = 1

a0 = rnorm(1, gamma/(1-phi), sqrt(sigmaSq / (1-phi^2)))
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

SVMmcmc = SVM_MCMC(y, 10000)


SVM = list()
SVMElbo = rep(0, 10)
for(i in 1:10){
  SVM[[i]] = SVM_SGA(y, T, 25, 5000, c(0, mean(y)-1.27*2, 0.5, rep(0, T+1)), diag(0.1, T+4))
  SVMElbo[i] = tail(SVM[[i]]$ELBO, 1)
}
SVMfit = SVM[[which.max(SVMElbo)]]


ysupport = seq(min(y)-2, max(y) + 2, length.out=1000)
draws = mvtnorm::rmvnorm(10000, SVMfit$Mu, SVMfit$L %*% t(SVMfit$L))
alphaT1 = rep(0, 10000)

yden = function(y, alpha){
  1 / sqrt(2*pi) * exp(1/2 * (y - alpha) - 1/2 * exp(y - alpha))
}

ydensity = rep(0, 1000)
for(i in 1:10000){
  u = sample(10000, 1)
  alphaT1[i] = rnorm(1, draws[u, 2] + draws[u, 3]*draws[u, T+4], sqrt(exp(draws[u, 1])))
  ydensity = ydensity + yden(ysupport, alphaT1[i])/10000
}

alphaT1t = rep(0, 10000)
truedens = rep(0, 1000)
for(i in 1:10000){
  alphaT1t[i] = rnorm(1, 1 + 0.7*alpha[T], 1)
  truedens = truedens + yden(ysupport, alphaT1t[i])/10000
}
forecasts = data.frame(support = ysupport, true = truedens, VB = ydensity)
forecasts = tidyr::gather(forecasts, method, density, -support)
ggplot(forecasts) + geom_line(aes(support, density, colour=method))
