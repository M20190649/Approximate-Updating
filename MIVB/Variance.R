library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(pscl)
sourceCpp('Variance.cpp')

N = 500
sigmaSq = 1
reps = 9
fit = matrix(0, reps, 4)

ldgamma = function(x, alpha, beta){
  alpha * log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
}

results = data.frame()
stats = data.frame()
for(i in 1:reps){
  y = rnorm(N, 0, sigmaSq)
  posObs = 1/rgamma(10000, 2+N/2, 1+0.5*sum(y^2))
  lnmu = mean(log(posObs))
  lnvar = var(log(posObs))
  lndelta = log(sqrt(lnvar))
  mu = mean(posObs)
  v = var(posObs)
  iga = mu^2/v +2
  igb = mu*(mu^2/v + 1)
  VBln = SGA_Var_LN(y, 25, 5000, 0.25, lnmu, lndelta)
  VBig = SGA_Var_IG(y, 25, 5000, 0.25, iga, igb)
  support = seq(0.01, 3, length.out=500)
  dln = dlnorm(support, VBln$Params[1], exp(VBln$Params[2]))
  dig = exp(ldgamma(support, VBig$Params[1], VBig$Params[2]))
  dtrue = exp(ldgamma(support, 2+N/2, 1 + 0.5*sum(y^2)))
  df = data.frame(support, c(dln, dig, dtrue), rep(c('Lognormal', 'Inverse-Gamma', 'True'), rep(500, 3)), i)
  results = rbind(results, df)
  dfstats = data.frame(c(VBln$Iter, VBig$Iter), c(tail(VBln$ELBO, 1), tail(VBig$ELBO, 1)), c('ln', 'ig'))
  stats = rbind(stats, dfstats)
}

colnames(results) = c('Support', 'Density', 'Distribution', 'Dataset')
ggplot(results) + geom_line(aes(Support, Density, colour=Distribution)) + facet_wrap(~Dataset, scales='free')
