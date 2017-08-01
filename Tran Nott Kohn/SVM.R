library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
sourceCpp('ADPF.cpp')
sourceCpp('SVM_MCMC.cpp')
sourceCpp('ADDPF.cpp')

sigSq = 0.02
phi = 0.9
mu = -0.5
T = 500
S = 50

x = rep(0, T+S)
y = rep(0, T+S)
for(t in 1:(T+S)){
  if(t == 1){
    x[t] = rnorm(1, mu, sqrt(sigSq/(1-phi^2)))
  } else {
    x[t] = mu + phi*(x[t-1]-mu) + sqrt(sigSq) * rnorm(1)
  }
  y[t] = exp(x[t]/2) * rnorm(1)
}



yStar = log(y^2)
MCMC = SVM_MCMC(yStar[1:T], 50000, 0.05)


yM = matrix(y[1:T], T)
mean = c(-4.4, 0, 0.86, mu)
var = diag(c(0.7, 0.5, 0.101, sqrt(sigSq/(1-phi^2))))

lambda = cbind(mean, var)
fitBPF = VBIL_PF(yM, lambda, S=50, alpha=0.25, maxIter=5000, threshold=0.05, thresholdIS=0.9)

lambda = cbind(mean, var)
fitDPF = VBIL_DPF(matrix(yStar[1:T], T), lambda, S=10, N=100, alpha=0.25, maxIter=1, threshold=0.05, thresholdIS=0.9)






lambda = cbind(mean, var)
fitPartial = VBIL_PF(matrix(y[1:T], T), lambda, S=10, alpha=0.1, maxIter=5000, threshold=0.05, thresholdIS=0.9)

U = fitPartial$U
Sigma = t(U) %*% U
SigInv = solve(Sigma[1:3, 1:3])
muV = fitPartial$Mu
vComp = Sigma[4, 1:3] %*% SigInv
conVar = Sigma[4, 4] - vComp %*% Sigma[1:3, 4]

fitUpdate = VBIL_PF_Update(matrix(y[(T+1):(T+S)], S), lambda, matrix(diag(U)[1:3], 3), muV, SigInv, vComp, conVar, S=10, alpha=0.1, maxIter=5000, threshold=0.05, thresholdIS=0.9)






bestFit = which(finalELBO == max(finalELBO, na.rm=TRUE))
VBfit[[bestFit]]$Mu
VBfit[[bestFit]]$U




VBfit = VBfit[[bestFit]]

U = VBfit$U
U[is.na(U)] = 0
Sigma = t(U) %*% U
supSig = seq(0.001, 0.1, length.out=500)
supPhi = seq(0.001, 0.999, length.out=500)
supMu = seq(-1.5, 0, length.out=500)
supX = seq(x[T]-1.5, x[T]+1.5, length.out=500)
qplot(supSig, dlnorm(supPos, VBfit$Mu[1], sqrt(Sigma[1, 1])), geom='line')
qplot(supMu, dnorm(supR, VBfit$Mu[2], sqrt(Sigma[2, 2])), geom='line')
qplot(supPhi, dnorm(supR, VBfit$Mu[3], sqrt(Sigma[3, 3])), geom='line')
qplot(supX, dnorm(supR, VBfit$Mu[4], sqrt(Sigma[4, 4])), geom='line')


sigSq = 0.02
phi = 0.8
mu = -0.5
T = 500

set.seed(5)
reps = 500
results = data.frame()
for(i in 1:400){
  
  x0 = rnorm(1, mu, sqrt(sigSq/(1-phi^2)))
  x = rep(0, T+1)
  y = rep(0, T+1)
  for(t in 1:(T+1)){
    if(t == 1){
      x[t] = mu + phi*(x0-mu) + rnorm(1, 0, sqrt(sigSq))
    } else {
      x[t] = mu + phi*(x[t-1]-mu) + rnorm(1, 0, sqrt(sigSq))
    }
    y[t] = rnorm(1, 0, exp(x[t]/2))
  }
  yM = matrix(y[1:T], T)
  MCMC = SVM_MCMC(log(y^2)[1:T], 20000, 0.05)
  xsupport = seq(x[T]-1.5, x[T]+1.5, length.out=1000)
  ysupport = seq(min(y)-1.5, max(y)+1.5, length.out=1000)
  xDensity = rep(0, 1000)
  yDensity = rep(0, 1000)
  for(j in 2001:10000){
    mean = MCMC$theta[j,2] + MCMC$theta[j,3] * (MCMC$x[j, T+1] - MCMC$theta[j, 2])
    xDensity = xDensity + dnorm(xsupport, mean, sqrt(MCMC$theta[j, 1]))/8000
  }
  sumXT1 = cumsum(xDensity)
 
  for(j in 1:1000){
    u = runif(1)
    flag = TRUE
    k = 1
    while(flag){
      if(u < sumXT1[k]) {
        xT1draw = xsupport[k]
        flag = FALSE
      }
      k = k + 1
    }
    yDensity = yDensity + dnorm(ysupport, 0, exp(xT1draw/2))/1000
  }
  results = rbind(results, data.frame(Method = 'MCMC', 
                                     xLogScore = log(xDensity[min(which(x[T+1] < xsupport))]),
                                     yLogScore = log(yDensity[min(which(y[T+1] < ysupport))])))
  
  
  #EXTp1_1 = 0
  #EXTp1_2 = 0  
  #for(j in 1:1000){
  #  EXTp1_1 = EXTp1_1  + xDensity[j] * support[j] * (support[2] - support[1])
  #  EXTp1_2 = EXTp1_2  + xDensity[j] * support[j]^2 * (support[2] - support[1])
  #}
 
  #print(paste0('MCMC fininshed on iteration ', i))
  
  #mcmc = data.frame(method='MCMC', 
  #                  variable=c('sigma^2', 'mu', 'phi', 'X[T+1]'), 
  #                  statistic = c(rep('Mean', 4), rep('Variance', 4)),
  #                  value = c(colMeans(MCMC$theta), EXTp1_1, var(MCMC$theta[,1]), var(MCMC$theta[,2]), var(MCMC$theta[,3]), EXTp1_2 - EXTp1_1^2),
  #                  iteration = i)
  
  yM = matrix(y[1:T], T)
  mean = c(-4.4, 0, 0.86, mu)
  var = diag(c(0.7, 0.5, 0.101, sqrt(sigSq/(1-phi^2))))
  
  VBfit = list()
  finalELBO = vector(length=0)
  for(j in 1:3){
    lambda = cbind(mean, var)
    fit = VBIL_PF(yM, lambda, S=10, alpha=0.1, maxIter=5000, threshold=0.01, thresholdIS=0.9)
    finalELBO = c(finalELBO, fit$ELBO[fit$Iter])
    VBfit[[j]] = fit
  }
  finalELBO[finalELBO > 0] = -1000
  bestFit = which(finalELBO == max(finalELBO, na.rm=TRUE))
  VBfit = VBfit[[bestFit]]
  
  yDensity = rep(0, 1000)
  xT1draws = rnorm(1000, VBfit$Mu[4], sqrt(sum(VBfit$U[,4]^2)))
  for(j in 1:1000){
    yDensity = yDensity + dnorm(ysupport, 0, exp(xT1draws[j]/2))/1000
  }
  results = rbind(results, data.frame(Method = 'VB', 
                                      xLogScore = dnorm(x[T+1], VBfit$Mu[4], sqrt(sum(VBfit$U[,4]^2)), log=TRUE),
                                      yLogScore = log(yDensity[min(which(y[T+1] < ysupport))])))
  if(i %% 5 == 0){
    print(i)
  }
}

results = gather(results, variable, logscore, -Method)
results = rbind(results, results2)
write.csv(results, 'sim500.csv', row.names=FALSE)


results = read.csv('simulResults.csv')
resMCMC = filter(results, method=='MCMC')
results = rbind(resMCMC, resultsvb)

results %>% filter(statistic == "Mean") %>%
  ggplot() + geom_point(aes(method, value)) + 
  facet_wrap(~variable, scales='free', labeller = label_parsed) +
  labs(x=NULL, y=NULL, title='Posterior Mean') +
  theme(plot.title = element_text(hjust = 0.5))

results %>% filter(statistic == "Variance") %>% 
  mutate(value = sqrt(value)) %>% 
  ggplot() + geom_point(aes(method, value)) + 
  facet_wrap(~variable, scales='free', labeller = label_parsed) + 
  labs(x=NULL, y=NULL, title='Posterior Standard Deviation') +
  theme(plot.title = element_text(hjust = 0.5))


set.seed(12)
for(i in 1:reps){
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
  mean = c(-4, 0, 0.5, rep(0, T+2))
  var = c(0.1, 0.15, 0.15, rep(0.15, T+2))
  VBfit = list()
  finalELBO = vector(length=0)
  for(j in 1:10){
    lambda = matrix(c(mean, var), 2*T+10)
    fit = ADVI(yM, lambda, S=25, alpha=0.1, maxIter=5000, threshold=0.05, thresholdIS=0)
    finalELBO = c(finalELBO, fit$ELBO[fit$Iter])
    VBfit[[j]] = fit
  }
  bestFit = which(finalELBO == max(finalELBO, na.rm=TRUE))
  muV = VBfit[[bestFit]]$Mu
  sd = VBfit[[bestFit]]$Sd
  print(paste0('VB fininshed on iteration ', i))
  
  vbfull = data.frame(method='VB-Full',
                  variable = c('sigma^2', 'mu', 'phi', 'X[T+1]'),
                  statistic = c(rep('Mean', 4), rep('Variance', 4)),
                  value = c(exp(muV[1] + 0.5 * sd[1]^2), muV[2], 2*muV[3]-1, muV[T+5],(exp(sd[1]^2) -1) * exp(2*muV[1] + sd[1]^2), sd[2]^2, 4*sd[3]^2, sd[T+5]^2), 
                  iteration = i)
  results = rbind(results, vbfull)
  
}


# Conditional Distribution Calculation
theta = rnorm(3)
muCond = muV[4] + Sigma[4, 1:3] %*% solve(Sigma[1:3, 1:3]) %*% (theta - muV[1:3])
varCond = Sigma[4, 4] - Sigma[4, 1:3] %*% solve(Sigma[1:3, 1:3]) %*% Sigma[1:3, 4]



