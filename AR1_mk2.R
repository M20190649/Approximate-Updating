library(GGally)
library(coda)

T = 100    #Length of y
N = 50000   #MCMC draws

#True parameters
sigma2 <- 2
mu <- 2
phi <- 0.5

#Generate data
y = rnorm(1, mu, sqrt(sigma2/(1-phi^2)))
for(i in 2:T){
  y = c(y, mu + phi*(y[(i-1)]-mu) + rnorm(1, 0, sqrt(sigma2)))
}

#Prior parameters
mubar = 0
lambda = 10
a = 1
b = 1
p1 = 5
p2 = 5

#constant
shape = T/2 + a

#log phi kernel for the MH step
lphiden = function(draw, mu, sig2){
  t1 = 1/2*log(1-draw^2) -(1-draw^2)*(y[1]-mu)^2/(2*sig2)  #The normal exponential for y1
  t2 = -(sum(((y[2:T]-mu)-(draw*(y[1:(T-1)]-mu)))^2))/(2*sig2)    #The normal exponentials for the rest of the y's
  t3 = log(1+draw)*(p1-1) + log(1-draw)*(p2-1)   #The stretched beta prior
  return(t1+t2+t3)
}

#put draws here
theta = matrix(ncol=3, nrow=N)
theta[1,-2] = 1       #Starting parameters
theta[1, 2] = 0.5

#count accepted draws
accept = 0

for(i in 2:N){
  
  #mu conditional
  denom = lambda^2 * ((T-1)*(1-theta[(i-1),2])^2+(1-theta[(i-1),2]^2)) + theta[(i-1), 3]
  sigsqmu = theta[(i-1), 3]*lambda^2 / denom
  muhat = (lambda^2*((1-theta[(i-1),2]^2)*y[1] + (1-theta[(i-1),2])*sum(y[2:T]-theta[(i-1),2]*y[1:(T-1)])) + mubar*theta[(i-1), 3])/denom
  theta[i,1] = rnorm(1, muhat, sqrt(sigsqmu))
  
  #phi conditional
  phihat = sum((y[2:T]-theta[i, 1])*(y[1:(T-1)]-theta[i,1]))/sum((y[1:(T-1)]-theta[i,1])^2)   #from the Catherine's notes for ETC4541, pg 84
  vphi = theta[(i-1),3]/sum((y[1:(T-1)]-theta[i,1])^2)    #acceptance rate for this candidate ~75%, trying random walk instead
  phidraw = rnorm(1, phihat, sqrt(vphi))
  #phidraw = rnorm(1, theta[(i-1),2], 0.4)               #random walk MH candidate draw to tune acceptance ratios to ~25%
  
  #phi MH step
  if(phidraw > -1 & phidraw < 1) {                       #within stationary conditions
    mh.candidate = lphiden(phidraw, theta[i, 1], theta[(i-1), 3])     #symmetrical candidate so only need to evaluate the true density components
    mh.prev = lphiden(theta[(i-1), 2], theta[i, 1], theta[(i-1), 3])  #using log density
    alpha = min(1, exp(mh.candidate-mh.prev))
    u = runif(1)
    if(u <= alpha) {             #accept draw
      theta[i, 2] = phidraw
      accept = accept + 1
    } else {                    #reject draw
      theta[i, 2] = theta[(i-1), 2]
    }
  } else {                      #non stationary
    theta[i, 2] = theta[(i-1), 2]
  }
  
  #sig2 conditional
  inv.scale = b + ((y[1] - theta[i, 1])^2*(1-theta[i, 2]^2) + sum(((y[2:T] - theta[i, 1]) - theta[i, 2]*(y[1:(T-1)]-theta[i, 1]))^2))/2
  theta[i, 3] = 1/rgamma(1, shape = shape, rate = inv.scale)
}
accept/N
keep = theta[(N/5+1):N,] #Drop first 20% of draws
colnames(keep) = c("Mu", "Phi", "Sigma2")
effectiveSize(keep)
#ggpairs(keep)

thin = keep[seq(1, 4*N/5, 10),] #Keep every 10th
effectiveSize(thin)
ggpairs(thin)
