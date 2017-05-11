library(tidyr)
library(ggplot2)

lbsigy = function(mean){
  elbo = 0
  for(i in 1:100){
    sigx = exp(rnorm(1, mean, 0.1))
    dens = -2*log(sigx) - 1/sigx #- x0^2*(1-0.9^2)/(2*sigx) 
    #x = c(x0, x)
    for(t in 1:T){
      dens = dens + dnorm(y[t], 1+x[t], sqrt(sigx), log=TRUE)
    }
    dens = dens - dnorm(log(sigx), mean, 0.1, log = TRUE)
    elbo = elbo + dens/100
  }
  return(elbo)
}
lbsigx = function(mean){
  elbo = 0
  for(i in 1:100){
    sigx = exp(rnorm(1, mean, 0.1))
    dens = -2.5*log(sigx) - 1/sigx - x0^2*(1-0.9^2)/(2*sigx) 
    x = c(x0, x)
    for(t in 1:T){
      dens = dens + dnorm(x[t+1], 0.9*x[t], sqrt(sigx), log=TRUE)
    }
    dens = dens - dnorm(log(sigx), mean, 0.1, log = TRUE)
    elbo = elbo + dens/100
  }
  return(elbo)
}
lbphi = function(mean){
  elbo = 0
  for(i in 1:100){
    phi = rtruncnorm(1, -1, 1, mean, 0.1)
    dens = log(sqrt(1-phi^2)) - x0^2*(1-phi^2)/2 
    x = c(x0, x)
    for(t in 1:T){
      dens = dens + dnorm(x[t+1], phi*x[t], 1, log=TRUE)
    }
    dens = dens - dnorm(phi, mean, 0.1, log = TRUE)
    elbo = elbo + dens/100
  }
  return(elbo)
}
lbgamma = function(mean){
  elbo = 0
  
  for(i in 1:100){
    dens = 0
    gamma = rnorm(1, mean, 0.1)
    for(t in 1:T){
      dens = dens + dnorm(y[t], gamma+x[t], 1, log=TRUE)
    }
    dens = dens - dnorm(gamma, mean, 0.1, log = TRUE)
    elbo = elbo + dens/100
  }
  return(elbo)
}

mean = seq(-0.8, 3, 0.02)
dens = rep(0, length(mean))
meanphi = seq(-0.99, 0.99, 0.02)
densphi = rep(0, length(meanphi))

max = matrix(0, 50, 4)
for(k in 1:32){
  T=250
  x0 = rnorm(1, 0, sqrt(sigmaSqX/(1-phi^2)))
  x = rep(0, T)
  y = rep(0, T)
  for(t in 1:T){
    if(t == 1){
      x[t] = phi*x0 + rnorm(1, 0, sqrt(sigmaSqX)) 
    } else {
      x[t] = phi*x[t-1] + rnorm(1, 0, sqrt(sigmaSqX))
    }
    y[t] = mu + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
  }
  dens = rep(0, length(mean))
  for(i in 1:length(mean)){dens[i] = lbsigy(mean[i])}
  max[k,1] = mean[which.max(dens)]
  dens = rep(0, length(mean))
  for(i in 1:length(mean)){dens[i] = lbsigx(mean[i])}
  max[k,2] = mean[which.max(dens)]
  densphi = rep(0, length(meanphi))
  for(i in 1:length(mean)){densphi[i] = lbphi(meanphi[i])}
  max[k,3] = meanphi[which.max(densphi)]
  dens = rep(0, length(mean))
  for(i in 1:length(mean)){dens[i] = lbgamma(mean[i])}
  max[k,4] = mean[which.max(dens)]
  print(k)
}
maxdf = as.data.frame(max)
colnames(maxdf) = c('LSigySq', 'LSigxSq', 'Phi', 'Gamma')
maxdf = gather(maxdf, variable, maximum)
ggplot(maxdf, aes(x=maximum)) + facet_wrap(~variable, scales='free') + geom_histogram() + labs(y='ELBO Maximising Mean')