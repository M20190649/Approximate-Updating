library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
library(mvtnorm)
sourceCpp('ADPF.cpp')
sourceCpp('SVM_MCMC.cpp')
#sourceCpp('ADDPF.cpp')
sourceCpp('ADPFUpdater.cpp')
#sourceCpp('APF.cpp')

sigSq = 0.02
phi = 0.9
mu = -0.5
T = 500

x = rep(0, T)
y = rep(0, T)
for(t in 1:T){
  if(t == 1){
    x[t] = rnorm(1, mu, sqrt(sigSq/(1-phi^2)))
  } else {
    x[t] = mu + phi*(x[t-1]-mu) + sqrt(sigSq) * rnorm(1)
  }
  y[t] = exp(x[t]/2) * rnorm(1)
}

yStar = log(y^2)
MCMC = SVM_MCMC(yStar, 50000, 0.05)
yM = matrix(y, T)
mean = c(-4.4, 0, 0.86, mu)
var = diag(c(0.7, 0.5, 0.101, sqrt(sigSq/(1-phi^2))))
lambda = cbind(mean, var)
fitBPF = VBIL_PF(yM, lambda, S=10, alpha=0.25, maxIter=5000, threshold=0.05, thresholdIS=0.9)
lambda = cbind(mean, var)
fitAPF = VBIL_APF(yM, lambda, S=10, alpha=0.25, maxIter=5000, threshold=0.05, thresholdIS=0.9)

VBfit = fitBPF
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

set.seed(16)
reps = 200
results = data.frame()
for(i in 1:reps){
  x0 = rnorm(1, mu, sqrt(sigSq/(1-phi^2)))
  x = rep(0, T+10)
  y = rep(0, T+10)
  for(t in 1:(T+10)){
    if(t == 1){
      x[t] = mu + phi*(x0-mu) + rnorm(1, 0, sqrt(sigSq))
    } else {
      x[t] = mu + phi*(x[t-1]-mu) + rnorm(1, 0, sqrt(sigSq))
    }
    y[t] = rnorm(1, 0, exp(x[t]/2))
  }
  MCMC = SVM_MCMC(log(y^2)[1:T], 15000, 0.05)
  ysupport = seq(min(y)-1.5, max(y)+1.5, length.out=1000)
  phivec = rep(1, 10)
  for(h in 1:10){
    yDensity = rep(0, 1000)
    xTh = rep(0, 10000)
    for(j in 5001:15000){
      xTh[j] = MCMC$theta[j, 3] + MCMC$theta[j, 2]^(h-1) * (MCMC$x[j, T+1] - MCMC$theta[j, 3])
      if(h > 1){
        for(k in 2:h){
          xTh[j] = xTh[j] + MCMC$theta[j, 2]^(k-2) * rnorm(1, 0, sqrt(MCMC$theta[j, 1]))
        }
      }
      yDensity = yDensity + dnorm(ysupport, 0, exp(xTh[j]/2))/10000
    }
    logscore = log(yDensity[min(which(y[T+h] < ysupport))])
    yCDF = cumsum(yDensity) / sum(yDensity)
    CRPS = -sum((yCDF - (ysupport > y[T+h]))^2) * (ysupport[2] - ysupport[1]) 
    results = rbind(results,
                    data.frame(Method = 'MCMC',
                               LogScore = logscore,
                               CRPS = CRPS,
                               MSE = (log(y[T+1]^2) - mean(xTh) + 1.27)^2,
                               yTh = y[T+h],
                               h = h,
                               iter = i))
    
  }

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
  Sigma = t(VBfit$U) %*% VBfit$U
  yDensity = rep(0, 1000)
  qDraws = rmvnorm(1000, VBfit$Mu, Sigma)
  qDraws[1, ] = exp(qDraws[1,])
  for(h in 1:10){
    yDensity = rep(0, 1000)
    xTh = rep(0, 1000)
    for(j in 1:1000){
      xTh[j] = qDraws[j, 3] + qDraws[j, 2]^(h-1) * (qDraws[j, 4] - qDraws[j, 3])
      if(h > 1){
        for(k in 2:h){
          xTh[j] = xTh[j] + qDraws[j, 2]^(k-2) * rnorm(1, 0, sqrt(qDraws[j, 1]))
        }
      }
      yDensity = yDensity + dnorm(ysupport, 0, exp(xTh[j]/2))/10000
    }
    logscore = log(yDensity[min(which(y[T+h] < ysupport))])
    yCDF = cumsum(yDensity) / sum(yDensity)
    CRPS = -sum((yCDF - (ysupport > y[T+h]))^2) * (ysupport[2] - ysupport[1]) 
    results = rbind(results,
                    data.frame(Method = 'VB',
                               LogScore = logscore,
                               CRPS = CRPS,
                               MSE = (log(y[T+1]^2) - mean(xTh) + 1.27)^2,
                               yTh = y[T+h],
                               h = h,
                               iter = i))
  }
  print(i)
  if(i %% 25 == 0){
    write.csv(results, 'sim500Extra.csv', row.names=FALSE)
  }
}
results = read.csv('sim500Fixed.csv')
resultsL = gather(results, statistic, value, -Method)
ggplot(resultsL) + geom_boxplot(aes(Method, value)) + facet_wrap(~statistic, scales='free')



VB = filter(results, Method=='VB-BPF')
MCMC = filter(results, Method == 'MCMC')
resJoined = cbind(select(VB, -Method), select(MCMC, -Method))
colnames(resJoined) = c('LS_VB', 'CRPS_VB', 'MSE_VB', 'LS_MCMC', 'CRPS_MCMC', 'MSE_MCMC')
resJoined %>% mutate(Logscore = LS_VB - LS_MCMC,
                     CRPS = CRPS_VB - CRPS_MCMC, 
                     MSE = MSE_VB - MSE_MCMC) %>%
              select(Logscore, CRPS, MSE) %>%
              mutate(iter = 1:500) %>%
              gather(statistic, difference, -iter) %>%
              ggplot() + geom_density(aes(difference)) + facet_wrap(~statistic, scales='free')

# Conditional Distribution Calculation
theta = rnorm(3)
muCond = muV[4] + Sigma[4, 1:3] %*% solve(Sigma[1:3, 1:3]) %*% (theta - muV[1:3])
varCond = Sigma[4, 4] - Sigma[4, 1:3] %*% solve(Sigma[1:3, 1:3]) %*% Sigma[1:3, 4]

# Average time to converge
set.seed(128)
N = c(25, 50, 100)
time = data.frame()
for(j in 1:3){
  for(k in 1:500){
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
    mean = c(-4.4, 0, 0.86, mu)
    var = diag(c(0.7, 0.5, 0.101, sqrt(sigSq/(1-phi^2))))
    lambda = cbind(mean, var)
    clockTime = proc.time()
    fit = VBIL_PF(yM, lambda, S=1, N=N[j], alpha=0.1, maxIter=5000, threshold=0.01, thresholdIS=0.9)
    runTime = proc.time() - clockTime
      
    yDensity = rep(0, 1000)
    ysupport = seq(min(y)-1.5, max(y)+1.5, length.out=1000)
    xT1draws = rnorm(1000, fit$Mu[4], sqrt(sum(fit$U[,4]^2)))
    for(l in 1:1000){
      yDensity = yDensity + dnorm(ysupport, 0, exp(xT1draws[l]/2))/1000
    }
    time = rbind(time, data.frame(Method = 'BPF', runTime = runTime[3], S = S[i], N = N[j], iter = fit$Iter, logscore = log(yDensity[min(which(y[T+1] < ysupport))])))
    
    lambda = cbind(mean, var)
    clockTime = proc.time()
    fit = VBIL_APF(yM, lambda, S=1, N=N[j], alpha=0.1, maxIter=5000, threshold=0.01, thresholdIS=0.9)
    runTime = proc.time() - clockTime
    
    yDensity = rep(0, 1000)
    ysupport = seq(min(y)-1.5, max(y)+1.5, length.out=1000)
    xT1draws = rnorm(1000, fit$Mu[4], sqrt(sum(fit$U[,4]^2)))
    for(l in 1:1000){
      yDensity = yDensity + dnorm(ysupport, 0, exp(xT1draws[l]/2))/1000
    }
    time = rbind(time, data.frame(Method = 'APF', runTime = runTime[3], S = S[i], N = N[j], iter = fit$Iter, logscore = log(yDensity[min(which(y[T+1] < ysupport))])))
          if(k %% 5 == 0){
      print(paste(i, j, k, sep = " "))
    }
    if(k %% 50 == 0){
      write.csv(time, 'time.csv', row.names=FALSE)
    }
  }
}

time = read.csv('time.csv')

ggplot(time) + geom_point(aes(runTime, logscore, colour=Method))
ggplot(time) + geom_point(aes(runTime, logscore, colour=Method), alpha=0.7) + facet_grid(N ~ S, scales = 'free')
ggplot(time) + geom_boxplot(aes(Method, logscore)) + facet_grid(S ~ N, scales = 'free')
ggplot(filter(time, S == 10)) + geom_point(aes(runTime, logscore, colour=factor(N))) + facet_wrap(~Method, scales = 'free')
ggplot(filter(time, S == 10)) + geom_boxplot(aes(factor(N), logscore)) + facet_wrap(~Method)

time %>% gather(statistic, value, -Method, -N, -S, -iter) -> timelong

ggplot(filter(timelong, S==1)) + geom_boxplot(aes(factor(N), value)) +
  facet_grid(statistic~Method, scales='free_y', labeller = label_parsed) + 
  labs(x = 'Number of Particles', y = NULL)

# Updating Algorithm
S = 50
reps = 500
set.seed(112)
update = data.frame()
update = read.csv('update.csv')
for(i in 1:25){
  x = rep(0, T+S+1)
  y = rep(0, T+S+1)
  for(t in 1:(T+S+1)){
  if(t == 1){
    x[t] = rnorm(1, mu, sqrt(sigSq/(1-phi^2)))
  } else {
    x[t] = mu + phi*(x[t-1]-mu) + sqrt(sigSq) * rnorm(1)
  }
  y[t] = exp(x[t]/2) * rnorm(1)
}

  yStar = log(y[1:(T+S)]^2)
  MCMC = SVM_MCMC(yStar, 15000, 0.05)
  
  yDensity = rep(0, 1000)
  for(j in 5001:15000){
    yDensity = yDensity + dnorm(ysupport, 0, exp(MCMC$x[j,T+S+1]/2))/10000
  }
  logscore = log(yDensity[min(which(y[T+S+1] < ysupport))])
  yCDF = cumsum(yDensity) / sum(yDensity)
  CRPS = -sum((yCDF - (ysupport > y[T+S+1]))^2) * (ysupport[2] - ysupport[1]) 
  update = rbind(update, data.frame(Method = 'MCMC', LogScore = logscore, CRPS = CRPS, MSE = (log(y[T+S+1]^2) - mean(MCMC$x[5001:15000, T+S+1]) + 1.27)^2))
  
  ymT = matrix(y[1:T], T)
  ymS = matrix(y[1:S+T], S)
  mean = c(-4.4, 0, 0.86, mu)
  var = diag(c(0.7, 0.5, 0.101, sqrt(sigSq/(1-phi^2))))
  VBfit = list()
  finalELBO = vector(length=3)
  for(j in 1:3){
    lambda = cbind(mean, var)
    fit = VBIL_PF(ymT, lambda, S=1, alpha=0.25, maxIter=5000, threshold=0.01, thresholdIS=0.9)
    finalELBO[j] = fit$ELBO[fit$Iter]
    VBfit[[j]] = fit
  }
  finalELBO[finalELBO > 0] = -10000
  bestFit = min(which(finalELBO == max(finalELBO, na.rm=TRUE)))
  fitPartial = VBfit[[bestFit]]
  U = fitPartial$U
  Sigma = t(U) %*% U
  SigInv = solve(Sigma[1:3, 1:3])
  muV = fitPartial$Mu
  vComp = Sigma[4, 1:3] %*% SigInv
  conVar = Sigma[4, 4] - vComp %*% Sigma[1:3, 4]

  for(j in 1:3){
    lambda = cbind(mean, var)
    fit = VBIL_PF_Update(ymS, lambda, diag(U)[1:3], muV, SigInv, vComp, conVar, S=1, alpha=0.25, maxIter=5000, threshold=0.01, thresholdIS=0.9)
    finalELBO[j] = fit$ELBO[fit$Iter]
    VBfit[[j]] = fit
  }
  finalELBO[finalELBO > 0] = -10000
  bestFit = min(which(finalELBO == max(finalELBO, na.rm=TRUE)))
  VBfit = VBfit[[bestFit]]
  yDensity = rep(0, 1000)
  xT1draws = rnorm(1000, VBfit$Mu[4], sqrt(sum(VBfit$U[,4]^2)))
  for(j in 1:1000){
    yDensity = yDensity + dnorm(ysupport, 0, exp(xT1draws[j]/2))/1000
  }
  logscore = log(yDensity[min(which(y[T+S+1] < ysupport))])
  yCDF = cumsum(yDensity) / sum(yDensity)
  CRPS = -sum((yCDF - (ysupport > y[T+S+1]))^2) * (ysupport[2] - ysupport[1]) 
  update = rbind(update, data.frame(Method = 'VB-BPF', LogScore = logscore, CRPS = CRPS, MSE = (log(y[T+S+1]^2) - VBfit$Mu[4] + 1.27)^2))
  
  print(i)
  if(i %% 10 == 0){
  #  write.csv(update, 'update.csv', row.names=FALSE)
  }
}

update = read.csv('update.csv')
ggplot(update) + geom_boxplot(aes(Method, logscore)) + labs(y=(expression(paste(Y[T+S+1], ' predictive log-score'))))
