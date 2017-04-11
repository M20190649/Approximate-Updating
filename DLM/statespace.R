library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(ggplot2)
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
J = 10
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

# Markov Chain Monte Carlo

sourceCpp("DLM_MCMC.cpp")
MCMCdraws = DLM_MCMC(y[1:T], 50000)
thetaKeep = MCMCdraws$theta[25001:50000,]
xdrawKeep = MCMCdraws$x[25001:50000,]



# J + 1 step ahead forecast with no extra information
ysupport = seq(min(y)-1, max(y)+1, length.out=1000)
ydens = rep(0, 1000)
for(i in 1:10000){
  u = sample(25000, 1)
  gammaDraw = thetaKeep[u, 4]
  phiDraw = thetaKeep[u, 3]
  sigxDraw = thetaKeep[u, 2]
  sigyDraw = thetaKeep[u, 1]
  xTDraw = xdrawKeep[u, T+1]
  phiprod = rep(1, J+1)
  for(i in 2:(J+1)){
    phiprod[i] = phiDraw * phiprod[i-1]
  }
  ydens = ydens + dnorm(ysupport, gammaDraw + phi^J*xTDraw, sqrt(sigyDraw + sum(sigxDraw*phiprod)))/10000
}

# One step ahead forecast using y_{T+1:T+J}, not updating thetas
ydens2 = rep(0, 1000)
for(i in 1:10000){
  u = sample(25000, 1)
  gammaDraw = thetaKeep[u, 4]
  phiDraw = thetaKeep[u, 3]
  sigxDraw = thetaKeep[u, 2]
  sigyDraw = thetaKeep[u, 1]
  xTDraw = xdrawKeep[u, T+1]
  XTS = FFUpdatercpp(y[(T+1):(T+J)], phiDraw, gammaDraw, sigyDraw, sigxDraw, xTDraw)
  ydens2 = ydens2 + dnorm(ysupport, gammaDraw + phiDraw*XTS[1], sqrt(sigyDraw + phi^2 * XTS[2] + sigxDraw))/10000
}

# Rerun MCMC using all y up to time T+J
MCMCupdate = DLM_MCMC(y[1:(T+J)], 50000)
thetaUpdateKeep = MCMCupdate$theta[25001:50000,]
xdrawUpdateKeep = MCMCupdate$x[25001:50000,]

ydens3 = rep(0, 1000)
for(i in 1:10000){
  u = sample(25000, 1)
  gammaDraw = thetaUpdateKeep[u, 4]
  phiDraw = thetaUpdateKeep[u, 3]
  sigxDraw = thetaUpdateKeep[u, 2]
  sigyDraw = thetaUpdateKeep[u, 1]
  xTDraw = xdrawUpdateKeep[u, T+J+1]
  ydens3 = ydens3 + dnorm(ysupport, gammaDraw + phiDraw*xTDraw, sqrt(sigyDraw + sigxDraw))/10000
}

ggplot() + geom_line(aes(ysupport, ydens), colour = "red") +geom_line(aes(ysupport, ydens2), colour = "blue") +
  geom_line(aes(ysupport, ydens3), colour = "darkgreen") + geom_vline(aes(xintercept=y[61]))



# Variational Bayes

sourceCpp("DLM_SGA_FR.cpp")
set.seed(1)
FRVB = list()
FinalElbo = rep(0, 10)
for(i in 1:10){
  FRVB[[i]] = DLM_SGA(y=y[1:T], S=T, M=1, maxIter=5000, initialM = rep(0, 55), initialL = diag(0.1, 55))
  FinalElbo[i] = tail(FRVB[[i]]$ELBO, 1)
}
FRVBInit = FRVB[[which.max(FinalElbo)]]
UpdateM = c(FRVBInit$Mu, rep(0, J))
UpdateL = diag(0.1, T+J+5)
UpdateL[1:(T+5), 1:(T+5)] = FRVBInit$L
FRVBU = list()
UpdateElbo = rep(0, 10)
for(i in 1:10){
  FRVBU[[i]] = DLM_SGA(y=y[1:(T+J)], S=J, M=1, maxIter=5000, initialM=UpdateM, initialL=UpdateL)
  UpdateElbo[i] = tail(FRVBU[[i]]$ELBO, 1)
}
FRVBUpdate = FRVBU[[which.max(UpdateElbo)]]

MFVB = list()
FinalElbo = rep(0, 10)
for(i in 1:10){
  MFVB[[i]] = DLM_SGA(y=y[1:T], S=T, M=1, maxIter=5000, initialM = rep(0, 55), initialL = diag(0.1, 55), meanfield=TRUE)
  FinalElbo[i] = tail(MFVB[[i]]$ELBO, 1)
}
MFVBInit = MFVB[[which.max(FinalElbo)]]
UpdateM = c(MFVBInit$Mu, rep(0, J))
UpdateL = diag(0.1, T+J+5)
UpdateL[1:(T+5), 1:(T+5)] = MFVBInit$L
MFVBU = list()
UpdateElbo = rep(0, 10)
for(i in 1:10){
  MFVBU[[i]] = DLM_SGA(y=y[1:(T+J)], S=J, M=1, maxIter=5000, initialM=UpdateM, initialL=UpdateL, meanfield=TRUE)
  UpdateElbo[i] = tail(MFVBU[[i]]$ELBO, 1)
}
MFVBUpdate = MFVBU[[which.max(UpdateElbo)]]



# Vine Copula

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

VineNoIndep = RVineCopSelect(cbind(pobtheta, pobx), Matrix=VineMatrix, cores=4, trunclevel=5)





# Marginal Selection

fitdist(thetaKeep[,3], "truncnorm", fix.arg =  list(a = -1, b = 1), start = list(mean = 0, sd = 1), method = "mle")$bic
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
