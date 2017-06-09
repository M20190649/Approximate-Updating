library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
sourceCpp('Variance.cpp')


N = 500
sigmaSq = 1
reps = 6

# Log Density of Inverse Gamma Distribution. Using lgamma directly avoids gamma(alpha) = inf for high values of alpha
ldgamma = function(x, alpha, beta){
  alpha * log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
}



VB = data.frame()
True = data.frame()
stats = data.frame()
for(i in 1:reps){
  # Generate Data
  y = rnorm(N, 0, sigmaSq)

  # Fit from arbitary starting values
  lnmu = -1
  lndelta = log(0.5)
  iga = 0.5
  igb = 0.25
  
  # Fit the VB approx for both distributions
  VBln = SGA_Var(y=y, M=50, lognormal=TRUE, initPar1=lnmu, initPar2=lndelta, alpha=0.13)
  VBig = SGA_Var(y=y, M=50, lognormal=FALSE, initPar1=iga, initPar2=igb, alpha=0.5)
  
  # Evaluate density over a grid for fitted densities, starting densities and true density
  support = seq(0.01, 2, length.out=500)
  dln = dlnorm(support, VBln$Params[1], exp(VBln$Params[2]))
  dig = exp(ldgamma(support, exp(VBig$Params[1]), exp(VBig$Params[2])))
  dtrue = exp(ldgamma(support, 2+N/2, 1 + 0.5*sum(y^2)))
  dlnStart = dlnorm(support, lnmu, exp(lndelta))
  digStart = exp(ldgamma(support, exp(iga), exp(igb)))
  
  # Collect results
  df = data.frame(support, c(dln, dlnStart, dig, digStart),
                  rep(c('Lognormal', 'Inverse-Gamma'), rep(1000, 2)),
                  rep(c('Converged', 'Starting'), rep(500, 2)), i)
  VB = rbind(VB, df)
  df = data.frame(support, dtrue, rep(c('Lognormal', 'Inverse-Gamma'), rep(500, 2)), i)
  True = rbind(True, df)
  dfstats = data.frame(c(VBln$Iter, VBig$Iter), c(tail(VBln$ELBO, 1), tail(VBig$ELBO, 1)), c('ln', 'ig'), i)
  stats = rbind(stats, dfstats)
}

colnames(VB) = c('Support', 'Density', 'Distribution', 'Version', 'Dataset')
colnames(True) = c('Support', 'Density', 'Distribution', 'Dataset')
colnames(stats) = c('Iterations', 'ELBO', 'Distribution', 'Dataset')
ggplot() + geom_line(data=VB, aes(Support, Density, colour=Version)) + 
  geom_line(data=True, aes(Support, Density)) +
  facet_grid(Dataset~Distribution, scales='free')
