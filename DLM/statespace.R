library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
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
T = 50
J = 50
set.seed(5)

x0 = rnorm(1, 0, sqrt(sigmaSqX))
x = rep(0, T+J)
y = rep(0, T+J)
for(t in 1:(T+J)){
  if(t == 1){
    x[t] = phi*x0 + rnorm(1, 0, sqrt(sigmaSqX)) 
  } else {
    x[t] = phi*x[t-1] + rnorm(1, 0, sqrt(sigmaSqX))
  }
  y[t] = mu + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
}
sourceCpp("DLM_MCMC.cpp")
MCMCdraws = DLM_MCMC(y[1:T], 50000)
thetaKeep = MCMCdraws$theta[25001:50000,]
xdrawKeep = MCMCdraws$x[25001:50000,]
#MVFB

initialMean = apply(cbind(log(thetaKeep[,1:2]), thetaKeep[,3:4], xdrawKeep), 2, mean)
initialSd = apply(cbind(log(thetaKeep[,1:2]), thetaKeep[,3:4], xdrawKeep), 2, sd)

sourceCpp("DLM_SGA.cpp")
MFVB = DLM_SGA(y[1:T], 5, 5000, 0.05, cbind(initialMean, initialSd), TRUE)


sourceCpp("DLM_SGA_FR.cpp")
FRVBadam = DLM_SGA(y[1:T], 1, 5000, 0.01, 0.001)



pobtheta = pobs(thetaKeep[sample(1:25000, 1000),])
pobx = pobs(xdrawKeep[sample(1:25000, 1000),])

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
