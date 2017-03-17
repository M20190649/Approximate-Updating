library(ggplot2)
library(GGally)
library(coda)
library(truncnorm)
library(VineCopula)
library(fitdistrplus)
library(mixtools)
library(QRM)
library(pscl)

set.seed(2)
T = 50
mu = 2
phi = 0.5
sigmaSqY = 1
sigmaSqX = 1
alphaY = 1
betaY = 1
alphaX = 1
betaX = 1
muBar = 0
muVar = 10

x0 = rnorm(1, 0, sqrt(sigmaSqX))
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

FFBS = function(y, theta){
  T = length(y)
  phi = theta[1]
  mu = theta[2]
  sigmaSqY = theta[3]
  sigmaSqX = theta[4]
  att = rep(0, T)
  ptt = rep(0, T)
  y = y - mu
  at = 0
  pt = phi^2*sigmaSqX + sigmaSqX
  vt = y[1] 
  att[1] = pt*vt/(pt + sigmaSqY)
  ptt[1] = pt - pt^2/(pt + sigmaSqY)
  for(t in 2:T){
    at = phi*att[t-1]
    pt = phi^2*ptt[t-1] + sigmaSqX
    vt = y[t] - at
    att[t] = at + pt*vt/(pt + sigmaSqY)
    ptt[t] = pt - pt^2/(pt + sigmaSqY)
  }
  
  alpha = rep(0, T)
  alpha[T] = rnorm(1, att[T], sqrt(ptt[T]))
  for(t in (T-1):1){
    vstar = alpha[t+1] - phi*att[t]
    fstar = phi^2*ptt[t] + sigmaSqX
    mstar = ptt[t]*phi
    atT = att[t] + mstar*vstar/fstar
    ptT = ptt[t] - mstar^2/fstar
    alpha[t] = rnorm(1, atT, sqrt(ptT))
  }
  a0T = phi*alpha[1]/(phi^2+1)
  p0T = sigmaSqX/(phi^2 + 1)
  alpha0 = rnorm(1, a0T, sqrt(p0T))
  return(c(alpha0, alpha))
}

posterior.statistics = function(x){
  l95 = quantile(x, probs = 0.025)
  mean = mean(x)
  median = median(x)
  u95 = quantile(x, probs = 0.975)
  return(c(l95, mean, median, u95))
}

rep = 100000
xdraw = matrix(0, nrow = rep, ncol = T+1)
theta = matrix(0, ncol = 4, nrow = rep) #phi, mu, sigmaSqY, sigmaSqX

theta[1,] = c(phi, mu, sigmaSqY, sigmaSqX)

for(i in 2:rep){
  #states from the Kalman Filter
  xdraw[i, ] = FFBS(y, theta[i-1, ])
  #phi from trunc normal
  meanPhi = sum(xdraw[i,1:T]*xdraw[i,2:(T+1)]) / sum(xdraw[i,1:T]^2)
  varPhi = theta[i-1, 4] / sum(xdraw[i,1:T]^2)
  theta[i, 1] = rtruncnorm(1, -1, 1, meanPhi, sqrt(varPhi))
  #mu from normal
  meanMu = (theta[i-1, 3]*muBar + muVar*(sum(y - xdraw[i, 2:(T+1)]))) / (muVar*T + theta[i-1, 3])
  varMu = muVar*theta[i-1, 3] / (muVar*T + theta[i-1, 3])
  theta[i, 2] = rnorm(1, meanMu, sqrt(varMu))
  #sigmaSqY from invG
  theta[i, 3] = 1/rgamma(1, shape = T/2 + alphaY, rate = betaY + sum((y - xdraw[i, 2:(T+1)] - theta[i, 2])^2)/2)
  #sigmaSqX from invG
  theta[i, 4] = 1/rgamma(1, shape = (T+1)/2 + alphaX, rate = betaX + (xdraw[i, 1]^2 + sum((xdraw[i, 2:(T+1)] - xdraw[i, 1:T]*theta[i, 1])^2)/2))
}

thetaKeep = theta[seq(0.2*rep + 1, rep, length.out = 1000),]
colnames(thetaKeep) = c("Phi", "Mu", "SigY", "SigX")
effectiveSize(thetaKeep)
ggpairs(thetaKeep)
apply(thetaKeep, 2, posterior.statistics)
xdrawKeep = xdraw[seq(0.2*rep + 1, rep, length.out = 1000),]
colnames(xdrawKeep) = paste0("X", 0:T)
effectiveSize(xdrawKeep)
ggpairs(xdrawKeep[,1:5])
apply(xdrawKeep[,1:5], 2, posterior.statistics)


library(rstan)
DLM_data = list(T = T, y = y)
stanMCMC = stan(file = 'DLM.stan', data = DLM_data, iter = 50000, chains = 1)
print(stanMCMC, digits = 3)

model = stan_model('DLM.stan')
stanMF = vb(model, data = DLM_data, algorithm = "meanfield")
stanFull = vb(model, data = DLM_data, algorithm = "fullrank")

#MVFB

initialMean = apply(cbind(log(thetaKeep[,3:4]), thetaKeep[,1:2], xdrawKeep), 2, mean)
initialSd = apply(cbind(log(thetaKeep[,3:4]), thetaKeep[,1:2], xdrawKeep), 2, sd)
initialL = t(chol(cov(cbind(log(thetaKeep[,3:4]), thetaKeep[,1:2], xdrawKeep))))
options(warn = -1)
set.seed(3)
VB = list()
convergedLB = rep(0, 10)
for(i in 1:10){
  VB[[i]] = FRSGA(y, initialMean, initialL, 1, 0.2, 500, 0.1, 0.1, FALSE, TRUE)
  convergedLB[i] = VB[[i]]$LB[10]
}
BestVB = VB[[which.max(convergedLB)]]
FRVB = FRSGA
MFVB = MFSGA(y, initialMean, initialSd, 1, 0.1, 200, TRUE)
options(warn = 0)


mean = initialMean
sd = initialSd






pobtheta = pobs(thetaKeep)
pobx = pobs(xdrawKeep)

VineMatrix = matrix(0, T+5, T+5)
VineMatrix[T+5, ] = 1
VineMatrix[T+4, 1:(T+4)] = 2
VineMatrix[T+3, 1:(T+3)] = 3
VineMatrix[T+2, 1:(T+2)] = 4
VineMatrix[1:(T+1), 1] = c(T+5, 5:(T+4))
diag(VineMatrix[2:(T+1), 2:(T+1)]) = 5:(T+4)
for(t in 3:(T+1)){
  for(j in 2:(t-1)){
    VineMatrix[t, j] = T + j - t + 5
  }
}

Vine = RVineCopSelect(cbind(pobtheta, pobx), Matrix = VineMatrix, indeptest = TRUE, cores = 4)
Vinebackup = Vine
Vine$family[1:100,] = Vine$par[1:100,] = Vine$par2[1:100, ] = 0

fitdist(thetaKeep[,1], "truncnorm", fix.arg =  list(a = -1, b = 1), start = list(mean = 0, sd = 1), method = "mle")$bic
log(1000)*8 - 2*normalmixEM(thetaKeep[,1], k = 3)$loglik
log(1000)*5 - 2*normalmixEM(thetaKeep[,1], k = 2)$loglik
log(1000)*3 - 2*fit.st(thetaKeep[,2])$ll.max

fitdist(thetaKeep[,2], "norm", method = "mle")$bic
log(1000)*3 - 2*fit.st(thetaKeep[,2])$ll.max
log(1000)*8 - 2*normalmixEM(thetaKeep[,2], k = 3)$loglik
log(1000)*5 - 2*normalmixEM(thetaKeep[,2], k = 2)$loglik


Vine <- readRDS("Vine.rds")
table(Vine$family[101:105,])
table(Vine$family[101,1:100])
