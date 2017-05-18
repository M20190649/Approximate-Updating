library(tidyr)
library(ggplot2)

# Inputs: Range of mean for approx. distribution, variance of MCMC posterior, latent variables x0 & x, data y
lbsigy = function(meanrange, var, x0, x, y){
  # Calculate the ELBO at each value in meanrange
  elbo = rep(0, length(meanrange))
  # Use the same random numbers for each mean
  eps = rnorm(100)
  # For each mean value
  for(j in seq_along(meanrange)){
    # The ELBO is an average across 100 draws of sigma^2_y
    for(i in 1:100){
      sigy = exp(meanrange[j] + sqrt(var)*eps[i])
      # Prior density
      logp =  -2*log(sigy) - 1/sigy
      for(t in 1:T){
        # From p(y | theta)
        logp = logp + dnorm(y[t], 1+x[t], sqrt(sigy), log=TRUE)
      }
      # Approx density q(log(sigma^2_y))
      logq = dnorm(log(sigy), meanrange[j], sqrt(var), log=TRUE)
      # All other components of log(p) and log(q) are constant with respect to sigma^2_y
      elbo[j] = elbo[j] + (logp-logq)/100
    }
  }
  return(elbo)
}

# Other three parameters work the same way (except for truncated normal sampling of phi)
lbsigx = function(meanrange, var, x0, x, y){
  elbo = rep(0, length(meanrange))
  eps = rnorm(100)
  for(j in seq_along(meanrange)){
    for(i in 1:100){
      sigx = exp(meanrange[j] + sqrt(var)*eps[i])
      logp = -2.5*log(sigx) - 1/sigx - x0^2*(1-0.95^2)/(2*sigx) 
      x = c(x0, x)
      for(t in 1:T){
        logp = logp + dnorm(x[t+1], 0.95*x[t], sqrt(sigx), log=TRUE)
      }
    logq = dnorm(log(sigx), meanrange[i], sqrt(var), log=TRUE)
    elbo[j] = elbo[j] + (logp-logq)/100
    }
  }
  return(elbo)
}

lbphi = function(meanrange, var, x0, x, y){
  elbo = rep(0, length(meanrange))
  eps = rnorm(100)
  for(j in seq_along(meanrange)){
    for(i in 1:100){
    phi = meanrange[j] + sqrt(var)*eps[i]
    while(phi > 1 | phi < -1){
      phi = rnorm(1, meanrange[j], sqrt(var))
    }
    logp = log(sqrt(1-phi^2)) - x0^2*(1-phi^2)/2 
    x = c(x0, x)
    for(t in 1:T){
      logp = logp + dnorm(x[t+1], phi*x[t], 1, log=TRUE)
    }
    logq = dnorm(phi, meanrange[j], sqrt(var), log = TRUE)
    elbo[j] = elbo[j] + (logp-logq)/100
    }
  }
  return(elbo)
}
  
lbgamma = function(meanrange, var, x0, x, y){
  elbo = rep(0, length(meanrange))
  eps = rnorm(100)
  for(j in seq_along(meanrange)){
    for(i in 1:100){
      gamma = meanrange[j] + sqrt(var)*eps[i]
      logp = 0
      for(t in 1:T){
        logp = logp + dnorm(y[t], gamma+x[t], 1, log=TRUE)
      }
    logq = dnorm(gamma, meanrange[j], sqrt(var), log = TRUE)
    elbo[j] = elbo[j] + (logp-logq)/100
    }
  }
  return(elbo)
}

# All parameters except phi have this grid
mean = seq(-0.8, 3, 0.02)
# Phi uses this grid
meanphi = seq(0.5, 0.99, 0.005)

set.seed(5)
# Max[k, i] = ELBO maximising mean for theta_i in data replication k
# Function is fairly slow
max = matrix(0, 20, 4)
for(k in 1:20){
  T = 250
  x0 = rnorm(1, 0, sqrt(1/(1-0.95^2)))
  x = rep(0, T)
  y = rep(0, T)
  for(t in 1:T){
    if(t == 1){
      x[t] = 0.95*x0 + rnorm(1, 0, sqrt(1)) 
    } else {
      x[t] = 0.95*x[t-1] + rnorm(1, 0, sqrt(1))
    }
    y[t] = 2 + x[t] + rnorm(1, 0, sqrt(1))
  }
  MCMC = DLM_MCMC(y, 5000)
  
  elbosigy = lbsigy(mean, var(log(MCMC$theta[2501:5000,1])), x0, x, y)
  max[k, 1] = mean[which.max(elbosigy)]
  elbosigx = lbsigx(mean, var(log(MCMC$theta[2501:5000,2])), x0, x, y)
  max[k, 2] = mean[which.max(elbosigx)]
  elbophi = lbphi(meanphi, var(MCMC$theta[2501:5000,3]), x0, x, y)
  max[k, 3] = meanphi[which.max(elbophi)]
  elbogamma = lbgamma(mean, var(MCMC$theta[2501:5000,4]), x0, x, y)
  max[k, 4] = mean[which.max(elbogamma)]
  print(k)
}

maxdf = as.data.frame(max)
colnames(maxdf) = c('log(sigma[y]^2)', 'log(sigma[x]^2)', 'phi', 'gamma')
maxdf = gather(maxdf, variable, maximum)

truev = data.frame(true=c(0, 0, 0.95, 2), variable=c('log(sigma[y]^2)', 'log(sigma[x]^2)', 'phi', 'gamma'))
ggplot() + geom_histogram(data=maxdf, aes(x=maximum)) + 
  geom_vline(data=truev, aes(xintercept=true), colour='red') +
  facet_wrap(~variable, scales='free', labeller=label_parsed) +
  labs(x='ELBO Maximising Mean', y='Number of Occurences')
