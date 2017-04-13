library(Rcpp)
library(RcppArmadillo)
library(VineCopula)
library(ggplot2)
library(microbenchmark)
library(mvtnorm)
library(tidyr)
sourceCpp("DLM_MCMC.cpp")
sourceCpp("DLM_SGA_FR.cpp")

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
h = 1

set.seed(5)
x0 = rnorm(1, 0, sqrt(sigmaSqX))
x = rep(0, T+J+h)
y = rep(0, T+J+h)
for(t in 1:(T+J+h)){
  if(t == 1){
    x[t] = phi*x0 + rnorm(1, 0, sqrt(sigmaSqX)) 
  } else {
    x[t] = phi*x[t-1] + rnorm(1, 0, sqrt(sigmaSqX))
  }
  y[t] = mu + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
}

# Markov Chain Monte Carlo

MCMCdraws = DLM_MCMC(y[1:T], 50000)
thetaKeep = MCMCdraws$theta[25001:50000,]
xdrawKeep = MCMCdraws$x[25001:50000,]

# J + 1 step ahead forecast with no extra information
ysupport = seq(min(y)-5, max(y)+5, length.out=1000)
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
set.seed(1)
VB = list()
ELBOfit = data.frame()

Mvec = c(1, 5, 10, 25, 50, 100)
for(i in 301:600){
  if(i <= 300){
    VB[[i]] = DLM_SGA(y=y[1:T], S=T, M=Mvec[i%%6+1], maxIter=5000, initialM=c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL=diag(0.1, 55))
    df = data.frame(ELBO=tail(VB[[i]]$ELBO, 1), Iter=VB[[i]]$Iter, M=Mvec[i%%6+1], Covariance="Non-Diagonal", runtime=0)
    rownames(df) = NULL
    ELBOfit = rbind(ELBOfit, df)
  } else {
    VB[[i]] = DLM_SGA(y=y[1:T], S=T, M=Mvec[i%%6+1], maxIter=5000, initialM=c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL=diag(0.1, 55), meanfield=TRUE)
    df = data.frame(ELBO=tail(VB[[i]]$ELBO, 1), Iter=VB[[i]]$Iter, M=Mvec[i%%6+1], Covariance="Diagonal", runtime=0)
    rownames(df) = NULL
    ELBOfit = rbind(ELBOfit, df)
  }
  if(i %% 50 == 0){
    print(i)
  }
}
rownames(ELBOfit) = NULL
ggplot(ELBOfit, aes(x=Iter, y=ELBO, colour=factor(M))) + geom_point() + facet_wrap(~Covariance, scales = "free")

timings = microbenchmark(
               DLM_SGA(y=y[1:T], S=T, M=1, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=5, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=10, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=50, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=100, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=1, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=5, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=10, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=50, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=100, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=1, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=5, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=10, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=50, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=100, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55)),
               DLM_SGA(y=y[1:T], S=T, M=1, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=5, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=10, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=50, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=100, maxIter=1, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=1, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=5, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=10, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=50, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=100, maxIter=10, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=1, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=5, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=10, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=50, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               DLM_SGA(y=y[1:T], S=T, M=100, maxIter=20, initialM = c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE),
               times = 5L, control=list(order='inorder'))
runtime = data.frame(M=rep(c(1, 5, 10, 50, 100), 30),
                     Iter=rep(rep(c(1, 10, 20), rep(5, 3)), 10),
                     Covariance = rep(rep(c("Non-Diagonal", "Diagonal"), rep(15, 2)), 5),
                     time=timings$time) 
timeMod = lm(time ~ M*Iter*Covariance, data=runtime)
ELBOfit$runtime = predict(timeMod, newdata=ELBOfit)/1000000000
ggplot(ELBOfit, aes(x=runtime, y=ELBO, colour=factor(M))) + geom_point() + facet_wrap(~Covariance) +
  labs(x="Total Runtime (seconds)", y="Converged ELBO") + scale_color_discrete(name="Simulations per Iteration")



set.seed(1)
FRVB = list()
FinalElbo = rep(0, 10)
for(i in 1:10){
  FRVB[[i]] = DLM_SGA(y=y[1:T], S=T, M=100, maxIter=5000, initialM=c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL=diag(0.1, 55))
  FinalElbo[i] = tail(FRVB[[i]]$ELBO, 1)
}
FRVBInit = FRVB[[which.max(FinalElbo)]]
UpdateM = c(FRVBInit$Mu, rep(0, J))
UpdateL = diag(0.1, T+J+5)
UpdateL[1:(T+5), 1:(T+5)] = FRVBInit$L
FRVBU = list()
UpdateElbo = rep(0, 10)
for(i in 1:10){
  FRVBU[[i]] = DLM_SGA(y=y[1:(T+J)], S=J, M=100, maxIter=5000, initialM=UpdateM, initialL=UpdateL)
  UpdateElbo[i] = tail(FRVBU[[i]]$ELBO, 1)
}
FRVBUpdate = FRVBU[[which.max(UpdateElbo)]]

MFVB = list()
FinalElbo = rep(0, 10)
for(i in 1:10){
  MFVB[[i]] = DLM_SGA(y=y[1:T], S=T, M=5, maxIter=5000, c(0, 0, 0, mean(y[1:T]), rep(0, 51)), initialL = diag(0.1, 55), meanfield=TRUE)
  FinalElbo[i] = tail(MFVB[[i]]$ELBO, 1)
}
MFVBInit = MFVB[[which.max(FinalElbo)]]
UpdateM = c(MFVBInit$Mu, rep(0, J))
UpdateL = diag(0.1, T+J+5)
diag(UpdateL[1:(T+5), 1:(T+5)]) = MFVBInit$Sd
MFVBU = list()
UpdateElbo = rep(0, 10)
for(i in 1:10){
  MFVBU[[i]] = DLM_SGA(y=y[1:(T+J)], S=J, M=5, maxIter=5000, initialM=UpdateM, initialL=UpdateL, meanfield=TRUE)
  UpdateElbo[i] = tail(MFVBU[[i]]$ELBO, 1)
}
MFVBUpdate = MFVBU[[which.max(UpdateElbo)]]

#J + 1 step ahead FRVB Forecast no extra info

ysupport = seq(min(y)-5, max(y)+5, length.out=1000)
ydensVB = rep(0, 1000)
ydensMF = rep(0, 1000)
frSigmaInit = FRVBInit$L %*% t(FRVBInit$L)

for(i in 1:10000){
  draw = rmvnorm(1, c(FRVBInit$Mu), frSigmaInit)
  phiprod = rep(1, J+1)
  for(i in 2:(J+1)){
    phiprod[i] = draw[3] * phiprod[i-1]
  }
  ydensVB = ydensVB + dnorm(ysupport, draw[4] + draw[3]^J*draw[T+5], sqrt(exp(draw[1]) + sum(exp(draw[2])*phiprod)))/10000
}
for(i in 1:10000){
  sigyd = exp(rnorm(1, MFVBInit$Mu[1], MFVBInit$Sd[1]))
  sigxd = exp(rnorm(1, MFVBInit$Mu[2], MFVBInit$Sd[2]))
  phid = rnorm(1, MFVBInit$Mu[3], MFVBInit$Sd[3])
  gammad = rnorm(1, MFVBInit$Mu[4], MFVBInit$Sd[4])
  XTd = rnorm(1, MFVBInit$Mu[T+5], MFVBInit$Sd[T+5])
  phiprod = rep(1, J+1)
  for(i in 2:(J+1)){
    phiprod[i] = phid * phiprod[i-1]
  }
  ydensMF = ydensMF + dnorm(ysupport, gammad + phid^J*XTd, sqrt(sigyd + sum(sigxd*phiprod)))/10000
}

# One step ahead forecast using y_{T+1:T+J}, not updating thetas
ydensVB2 = rep(0, 1000)
ydensMF2 = rep(0, 1000)
for(i in 1:10000){
  sigyd = exp(rnorm(1, MFVBInit$Mu[1], MFVBInit$Sd[1]))
  sigxd = exp(rnorm(1, MFVBInit$Mu[2], MFVBInit$Sd[2]))
  phid = rnorm(1, MFVBInit$Mu[3], MFVBInit$Sd[3])
  gammad = rnorm(1, MFVBInit$Mu[4], MFVBInit$Sd[4])
  XTd = rnorm(1, MFVBInit$Mu[T+5], MFVBInit$Sd[T+5])
  XTS = FFUpdatercpp(y[(T+1):(T+J)], phid, gammad, sigyd, sigxd, XTd)
  ydensMF2 = ydensMF2 + dnorm(ysupport, gammad + phid*XTS[1], sqrt(sigyd + phid^2 * XTS[2] + sigxd))/10000
}
for(i in 1:10000){
  draw = rmvnorm(1, c(FRVBInit$Mu), frSigmaInit)
  XTS = FFUpdatercpp(y[(T+1):(T+J)], draw[3], draw[4], exp(draw[1]), exp(draw[2]), draw[T+5])
  ydensVB2 = ydensVB2 + dnorm(ysupport, draw[4] + draw[3]*XTS[1], sqrt(exp(draw[1]) + draw[3]^2 * XTS[2] + exp(draw[2])))/10000
}


# Update initial VB using all y_{T+1:T+J}
ydensVB3 = rep(0, 1000)
ydensMF3 = rep(0, 1000)
for(i in 1:10000){
  draw = rmvnorm(1, c(FRVBUpdate$Mu), frSigmaU)
  ydensVB3 = ydensVB3 + dnorm(ysupport, draw[4] + draw[3]*draw[T+J+5], sqrt(exp(draw[1]) + exp(draw[2])))/10000
}
for(i in 1:10000){
  sigyd = exp(rnorm(1, MFVBUpdate$Mu[1], MFVBUpdate$Sd[1]))
  sigxd = exp(rnorm(1, MFVBUpdate$Mu[2], MFVBUpdate$Sd[2]))
  phid = rnorm(1, MFVBUpdate$Mu[3], MFVBUpdate$Sd[3])
  gammad = rnorm(1, MFVBUpdate$Mu[4], MFVBUpdate$Sd[4])
  XTd = rnorm(1, MFVBUpdate$Mu[T+5], MFVBUpdate$Sd[T+5])
  ydensMF3 = ydensMF3 + dnorm(ysupport, gammad + phid*XTd, sqrt(sigyd + sigxd))/10000
}

MCMCforecast = gather(data.frame(ysupport, ydens, ydens2, ydens3), version, density, -ysupport)
MCMCforecast$method = "MCMC"
MFforecast = gather(data.frame(ysupport, ydensMF, ydensMF2, ydensMF3), version, density, -ysupport)
MFforecast$method = "Meanfield"
FRforecast = gather(data.frame(ysupport, ydensVB, ydensVB2, ydensVB3), version, density, -ysupport)
FRforecast$method = "FullRank"
forecasts = rbind(MCMCforecast, MFforecast, FRforecast)
forecasts$version = rep(rep(c("S+1 step", "filtered", "1 step"), rep(1000, 3)), 3)

ggplot(forecasts, aes(x=ysupport, y=density, colour=method)) + facet_wrap(~version) + geom_line()





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
