library(GGally)
library(coda)
library(ggplot2)
library(gridExtra)
library(reshape)
library(mvtnorm)

T = 1000    #Length of y
N = 15000   #MCMC draws

#True parameters
sigma2 <- 1
phi <- 0.97

#Generate data
y = rnorm(1, 0, sqrt(sigma2/(1-phi^2)))
for(i in 2:T){
  y = c(y, phi*y[(i-1)] + rnorm(1, 0, sqrt(sigma2)))
}

#Uninformative prior
p1 = 1
p2 = 1
a = 1
b = 1

#constant
shape = T/2 + a

#lphi kernel
lphiden = function(draw, sig2) {
  t1 = 1/2*log(1-draw^2) -(1-draw^2)*y[1]^2/(2*sig2)  #The normal exponential for y1
  t2 = -sum(((y[2:T]-draw*y[1:(T-1)])^2))/(2*sig2)    #The normal exponentials for the rest of the y's
  t3 = log(1+draw)*(p1-1) + log(1-draw)*(p2-1)   #The stretched beta prior
  return(t1+t2+t3)
}

theta = matrix(0, ncol = 2, nrow = N)
theta[1, ] = c(0.5, 1)
accept = 0

for(i in 2:N){
  #phi MH step
  #phihat = sum(y[2:T]*y[1:(T-1)])/sum(y[1:(T-1)]^2)
  #vphi = theta[(i-1),3]/sum(y[1:(T-1)]^2)
  #phidraw = rnorm(1, phihat, sqrt(vphi))
  phidraw = rnorm(1, theta[(i-1),1], 0.1)     
  
  if(phidraw > -1 & phidraw < 1) {                       #within stationary conditions
    mh.candidate = lphiden(phidraw, theta[(i-1), 2])     #symmetrical candidate so only need to evaluate the true density components
    mh.prev = lphiden(theta[(i-1), 1], theta[(i-1), 2])  #using log density
    alpha = min(1, exp(mh.candidate-mh.prev))
    u = runif(1)
    if(u <= alpha) {             #accept draw
      theta[i, 1] = phidraw
      accept = accept + 1
    } else {                    #reject draw
      theta[i, 1] = theta[(i-1), 1]
    }
  } else {                      #non stationary
    theta[i, 1] = theta[(i-1), 1]
  }
  
  #sigma squared conditional
  inv.scale = b + (y[1]^2*(1-theta[i, 2]^2) + sum((y[2:T] - theta[i, 1]*y[1:(T-1)])^2))/2
  theta[i, 2] = 1/rgamma(1, shape = shape, rate = inv.scale)
}

accept/N
keep = theta[(N/5+1):N,] #Drop first 20% of draws
colnames(keep) = c("Phi", "Sigma2")
effectiveSize(keep)
#ggpairs(keep)

thin = keep[seq(1, 4*N/5, 10),] #Keep every 10th
effectiveSize(thin)
ggpairs(thin)
