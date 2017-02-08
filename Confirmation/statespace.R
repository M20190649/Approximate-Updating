library(FKF)
library(truncnorm)
library(ggplot2)
T = 150
set.seed(51)

mu = 3
rho = 0.8
sigmaSqY = 1
sigmaSqX = 1

x0 = rnorm(1)
x = rep(0, T)
y = rep(0, T)
for(t in 1:T){
  if(t == 1){
    x[t] = rho*x0 + rnorm(1, 0, sigmaSqX) 
  } else {
    x[t] = rho*x[t-1] + rnorm(1, 0, sigmaSqX)
  }
  y[t] = mu + x[t] + rnorm(1, 0, sigmaSqY)
}

#add sigmas
forwardsFilter = function(y, T, rho, mu, sigy, sigx){
  att = rep(0, T)
  ptt = rep(0, T)
  y = y - mu
  at = 0
  pt = rho^2 + 1
  vt = y[1] - at
  att[1] = at + pt*vt/(pt + 1)
  ptt[1] = pt - pt^2/(pt + 1)
  for(t in 2:T){
    at = rho*att[t-1]
    pt = rho^2*ptt[t-1] + 1
    vt = y[t] - at
    att[t] = at + pt*vt/(pt + 1)
    ptt[t] = pt - pt^2/(pt + 1)
  }
  return(list(att = att, ptt = ptt))
}

#add sigmas
backwardsSampler = function(kalman, T, rho, sigx){
  alpha = rep(0, T+1)
  alpha[T+1] = rnorm(1, kalman$att[T], sqrt(kalman$ptt[T]))
  for(t in T:2){
    vstar = alpha[t+1] - rho*kalman$att[t]
    fstar = rho^2*kalman$ptt[t] + 1
    mstar = kalman$ptt[t]*rho
    atT = kalman$att[t] + mstar*vstar/fstar
    ptT = kalman$ptt[t] - mstar^2/fstar
    alpha[t] = rnorm(1, atT, sqrt(ptT))
  }
  a0T = rho*alpha[2]/(rho^2+1)
  p0T = 1 - rho^2/(rho^2 + 1)
  alpha[1] = rnorm(1, a0T, sqrt(p0T))
  return(alpha)
}

rep = 5000
xdraw = matrix(0, nrow = rep, ncol = T+1)
theta = matrix(0, ncol = 4, nrow = rep) #rho, mu, sigy, sigx

theta[1,] = 0.5

for(i in 2:rep){
  kalman = forwardsFilter(y[1:T], T, theta[i-1, 1], theta[i-1, 2], theta[i-1, 3], theta[i-1, 4])
  
  adraw[i, ] = backwardsSampler(kalman, T, theta[i-1], theta[i-1, 4])
  
  #rho from trunc normal
  theta[i, 1] = rtruncnorm(1, -1, 1, sum(xdraw[i,1:T]*xdraw[i,2:(T+1)]) / sum(xdraw[i,1:T]^2), sqrt(1/sum(xdraw[i,1:T]^2)))
  #mu from normal
  theta[i, 2] = rnorm(1)
  #sigy from InvG
  theta[i, 3] = 1/rgamma(1)
  #sigx from InvG
  theta[i, 4] = 1/rgamma(1)
}

