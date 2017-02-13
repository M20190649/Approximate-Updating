library(ggplot2)
library(GGally)
library(coda)
library(truncnorm)

set.seed(31)
T = 100
mu = 3
phi = 0.5
sigmaSqY = 1
sigmaSqX = 1
alphay = 1
betay = 1
alphax = 1
betax = 1

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

FFBS = function(y, T, phi, mu, sigmaSqY, sigmaSqX){
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

rep = 100000
xdraw = matrix(0, nrow = rep, ncol = T+1)
theta = matrix(0, ncol = 4, nrow = rep) #phi, mu, sigmaSqY, sigmaSqX

theta[1,] = c(phi, mu, sigmaSqY, sigmaSqX)

for(i in 2:rep){
  #states from the Kalman Filter
  xdraw[i, ] = FFBS(y, theta[i-1, ])
  #phi from trunc normal
  theta[i, 1] = rtruncnorm(1, -1, 1, sum(xdraw[i,1:T]*xdraw[i,2:(T+1)]) / sum(xdraw[i,1:T]^2), sqrt(theta[i-1, 4]/sum(xdraw[i,1:T]^2)))
  #mu from normal
  theta[i, 2] = rnorm(1, sum(y - xdraw[i, 2:(T+1)])/T, sqrt(theta[i-1, 3]/T))
  #sigmaSqY from invG
  theta[i, 3] = 1/rgamma(1, shape = T/2 + alphay, rate = betay + sum((y - xdraw[i, 2:(T+1)] - theta[i, 2])^2)/2)
  #sigmaSqX from invG
  theta[i, 4] = 1/rgamma(1, shape = (T+1)/2 + alphax, rate = betax + (xdraw[i, 1]^2 + sum((xdraw[i, 2:(T+1)] - xdraw[i, 1:T]*theta[i, 1])^2)/2))
}

thetaKeep = theta[seq(0.2*rep + 1, rep, length.out = 1000),]
colnames(thetaKeep) = c("Phi", "Mu", "SigY", "SigX")
effectiveSize(thetaKeep)
ggpairs(thetaKeep)
xdrawKeep = xdraw[seq(0.2*rep + 1, rep, length.out = 1000),]
colnames(xdrawKeep) = paste0("X", 0:T)
effectiveSize(xdrawKeep)
ggpairs(xdrawKeep[,1:5])

posterior.statistics = function(x){
  l95 = quantile(x, probs = 0.025)
  mean = mean(x)
  median = median(x)
  u95 = quantile(x, probs = 0.975)
  return(c(l95, mean, median, u95))
}

#xd = t(replicate(1000, FFBS(y, T, phi, mu, sigmaSqY, sigmaSqX), simplify = "matrix"))
#l95 = apply(xd[,2:151], 2, quantile, probs = 0.025)
#u95 = apply(xd[,2:151], 2, quantile, probs = 0.975)
#sum(l95 < x & x < u95)/T
#ggplot() + geom_line(aes(1:150, l95), colour = "red") + geom_line(aes(1:150, x)) + geom_line(aes(1:150, u95), colour = "red")

apply(thetaKeep, 2, posterior.statistics)
apply(xdrawKeep[,1:5], 2, posterior.statistics)

library(VineCopula)
pobtheta = pobs(thetaKeep)
VineTheta = RVineStructureSelect(pobtheta, indeptest = TRUE, cores = 4)
VineTheta$Matrix
VineTheta$family

library(fitdistrplus)
library(mixtools)
library(QRM)
fitdist(thetaKeep[,1], "truncnorm", fix.arg =  list(a = -1, b = 1), start = list(mean = 0, sd = 1), method = "mle")$bic
log(1000)*8 - 2*normalmixEM(thetaKeep[,1], k = 3)$loglik
log(1000)*5 - 2*normalmixEM(thetaKeep[,1], k = 2)$loglik
log(1000)*3 - 2*fit.st(thetaKeep[,2])$ll.max

fitdist(thetaKeep[,2], "norm", method = "mle")$bic
log(1000)*3 - 2*fit.st(thetaKeep[,2])$ll.max
log(1000)*8 - 2*normalmixEM(thetaKeep[,2], k = 3)$loglik
log(1000)*5 - 2*normalmixEM(thetaKeep[,2], k = 2)$loglik

#q(phi | truncnorm) * q(mu | t) * q(sigx | IG) * q(sigy | IG) * c(q sigx sigy | gaussian) * c(q phi q sigx Rotated BB8)