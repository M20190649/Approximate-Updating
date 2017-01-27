library(FKF)
library(truncnorm)
library(ggplot2)
T = 150
J = 150

rho = 0.8
beta = 1

a0 = rnorm(1)
a = rep(0, T+J)
y = rep(0, T+J)
for(t in 1:300){
  if(t == 1){
    a[t] = rho*a0 + rnorm(1) 
  } else {
    a[t] = rho*a[t-1] + rnorm(1)
  }
  y[t] = beta*a[t] + rnorm(1)
}

forwardsFilter = function(y, T, rho, beta){
  att = rep(0, T)
  ptt = rep(0, T)
  
  at = 0
  pt = rho^2 + 1
  vt = y[1] - beta*at
  ft = beta^2*pt + 1
  mt = pt*beta
  att[1] = at + mt*vt/ft
  ptt[1] = pt - mt^2/ft
  for(t in 2:T){
    at = rho*att[t-1]
    pt = rho^2*ptt[t-1] + 1
    vt = y[t] - beta*at
    ft = beta^2*pt + 1
    mt = pt*beta
    att[t] = at + mt*vt/ft
    ptt[t] = pt - mt^2/ft
  }
  return(list(att = att, ptt = ptt))
}

backwardsSampler = function(kalman, T, rho){
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
adraw = matrix(0, nrow = rep, ncol = T+1)
theta = rep(0, rep)

theta[1] = 0.5

for(i in 2:rep){
  kalman = forwardsFilter(y[1:T], T, theta[i-1], beta)
  
  adraw[i, ] = backwardsSampler(kalman, T, theta[i-1])
  
  theta[i] = rtruncnorm(1, -1, 1, sum(adraw[i,1:T]*adraw[i,2:(T+1)]) / sum(adraw[i,1:T]^2), sqrt(1/sum(adraw[i,1:T]^2)))
}

