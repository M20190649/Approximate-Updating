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
phi = 0.8
sigmaSqY = 1
sigmaSqX = 1
alphaY = 1
betaY = 1
alphaX = 1
betaX = 1
muBar = 0
muVar = 10

KLdiv = function(p, q){
  n = length(p)
  out = 0
  for(i in 1:n){
    out = q[i] * log(q[i] / p[i])
  }
  return(out)
}

predictiveDensities = function(T, S, h, MCMCrep){
  x0 = rnorm(1, 0, sqrt(sigmaSqX))
  x = rep(0, T+S+h)
  y = rep(0, T+S+h)
  for(t in 1:(T+S+h)){
    if(t == 1){
      x[t] = phi*x0 + rnorm(1, 0, sqrt(sigmaSqX)) 
    } else {
      x[t] = phi*x[t-1] + rnorm(1, 0, sqrt(sigmaSqX))
    }
    y[t] = mu + x[t] + rnorm(1, 0, sqrt(sigmaSqY))
  }
  
  #ggplot() + geom_line(aes(x=1:T, y=y[1:T])) + geom_line(aes(x=T:(T+S), y=y[T:(T+S)]), colour="red") + 
  #  geom_line(aes(x=(T+S):(T+S+h), y=y[(T+S):(T+S+h)]), colour='blue') + 
  #  theme(legend.position='none') + labs(x='t', y='y')
  
  # Markov Chain Monte Carlo
  
  MCMCdraws = DLM_MCMC(y[1:T], MCMCrep)
  thetaKeep = MCMCdraws$theta[(MCMCrep/2 + 1):MCMCrep,]
  xdrawKeep = MCMCdraws$x[(MCMCrep/2 + 1):MCMCrep,]
  
  #ggplot() + geom_line(aes(ysupport, ydens), colour = "red") +geom_line(aes(ysupport, ydens2), colour = "blue") +
  #  geom_line(aes(ysupport, ydens3), colour = "darkgreen") + geom_vline(aes(xintercept=y[T+S+h]))
  
  initM = c(0, 0, cor(y[2:T],y[2:T-1]), mean(y), rep(0, T+1))
  FRVBInit = DLM_SGA(y=y[1:T], S=T, M=100, maxIter=5000, initialM=initM, initialL=diag(0.1, T+5))
  MFVBInit = DLM_SGA(y=y[1:T], S=T, M=5, maxIter=5000, initialM=initM, initialL = diag(0.1, T+5), meanfield=TRUE)
  
  #S + 1 step ahead FRVB Forecast no extra info
  
  ysupport = seq(min(y)-5, max(y)+5, length.out=1000)
  frSigmaInit = FRVBInit$L %*% t(FRVBInit$L)
  
  ydensVB = rep(0, 1000)
  ydensMF = rep(0, 1000)
  for(i in 1:500){
    draw = rmvnorm(1, c(FRVBInit$Mu), frSigmaInit)
    phiprod = rep(1, S+h)
    for(i in 2:(S+h)){
      phiprod[i] = draw[3]^2 * phiprod[i-1]
    }
    ydensVB = ydensVB + dnorm(ysupport, draw[4] + draw[3]^(S+h)*draw[T+5], sqrt(exp(draw[1]) + sum(exp(draw[2])*phiprod)))/500
  }
  for(i in 1:500){
    sigyd = exp(rnorm(1, MFVBInit$Mu[1], MFVBInit$Sd[1]))
    sigxd = exp(rnorm(1, MFVBInit$Mu[2], MFVBInit$Sd[2]))
    phid = rnorm(1, MFVBInit$Mu[3], MFVBInit$Sd[3])
    gammad = rnorm(1, MFVBInit$Mu[4], MFVBInit$Sd[4])
    XTd = rnorm(1, MFVBInit$Mu[T+5], MFVBInit$Sd[T+5])
    phiprod = rep(1, S+h)
    for(i in 2:(S+h)){
      phiprod[i] = phid^2 * phiprod[i-1]
    }
    ydensMF = ydensMF + dnorm(ysupport, gammad + phid^(S+h)*XTd, sqrt(sigyd + sum(sigxd*phiprod)))/500
  }
  ydens = rep(0, 1000)
  for(i in 1:500){
    u = sample(MCMCrep/2, 1)
    gammaDraw = thetaKeep[u, 4]
    phiDraw = thetaKeep[u, 3]
    sigxDraw = thetaKeep[u, 2]
    sigyDraw = thetaKeep[u, 1]
    xTmean = xdrawKeep[u, T+2]
    xTvar = xdrawKeep[u, T+3]
    phiprod = rep(1, S+h)
    for(i in 2:(S+h)){
      phiprod[i] = phiDraw^2 * phiprod[i-1]
    }
    ydens = ydens + dnorm(ysupport, gammaDraw + phiDraw^(S+h)*xTmean, 
                          sqrt(sigyDraw + sum(sigxDraw*phiprod) + phiDraw^(2*(S+h))*xTvar))/500
  }
  obs = min(which((y[T+S+h] < ysupport) == TRUE))
  logscoreSh =  c(log(ydens[obs]), log(ydensMF[obs]), log(ydensVB[obs]))
  meanMCMC = colMeans(cbind(log(thetaKeep[,1:2]), thetaKeep[,3:4], xdrawKeep[,T+1]))
  meanVB = c(FRVBInit$Mu)[c(1:4, T+5)]
  meanMF = c(MFVBInit$Mu)[c(1:4, T+5)]

  # One step ahead forecast using y_{T+1:T+S}, not updating thetas
  yden = rep(0, 1000)
  ydensMF2 = rep(0, 1000)
  for(i in 1:500){
    sigyd = exp(rnorm(1, MFVBInit$Mu[1], MFVBInit$Sd[1]))
    sigxd = exp(rnorm(1, MFVBInit$Mu[2], MFVBInit$Sd[2]))
    phid = rnorm(1, MFVBInit$Mu[3], MFVBInit$Sd[3])
    gammad = rnorm(1, MFVBInit$Mu[4], MFVBInit$Sd[4])
    XTS = FFUpdatercpp(y[(T+1):(T+S)], phid, gammad, sigxd, sigyd, MFVBInit$Mu[T+5], MFVBInit$Sd[T+5]^2)
    ydensMF = ydensMF + dnorm(ysupport, gammad + phid*XTS[1], sqrt(sigxd + phid^2 * XTS[2] + sigyd))/500
  }
  for(i in 1:500){
    draw = rmvnorm(1, c(FRVBInit$Mu), frSigmaInit)
    XTS = FFUpdatercpp(y[(T+1):(T+S)], draw[3], draw[4], exp(draw[2]), exp(draw[1]), FRVBInit$Mu[T+5], frSigmaInit[T+5, T+5])
    ydensVB = ydensVB + dnorm(ysupport, draw[4] + draw[3]*XTS[1], sqrt(exp(draw[1]) + draw[3]^2 * XTS[2] + exp(draw[2])))/500
  }
  ydens = rep(0, 1000)
  for(i in 1:500){
    u = sample(MCMCrep/2, 1)
    gammaDraw = thetaKeep[u, 4]
    phiDraw = thetaKeep[u, 3]
    sigxDraw = thetaKeep[u, 2]
    sigyDraw = thetaKeep[u, 1]
    xTmean = xdrawKeep[u, T+2]
    xTvar = xdrawKeep[u, T+3]
    XTS = FFUpdatercpp(y[(T+1):(T+S)], phiDraw, gammaDraw, sigyDraw, sigxDraw, xTmean, xTvar)
    ydens = ydens + dnorm(ysupport, gammaDraw + phiDraw*XTS[1], sqrt(sigyDraw + phiDraw^2 * XTS[2] + sigxDraw))/500
  }
  obs = min(which((y[T+S+h] < ysupport) == TRUE))
  logscoreFilter =  c(log(ydens[obs]), log(ydensMF[obs]), log(ydensVB[obs]))
  
  UpdateM = c(FRVBInit$Mu, rep(0, S))
  UpdateL = diag(0.1, T+S+5)
  UpdateL[1:(T+5), 1:(T+5)] = FRVBInit$L
  FRVBUpdate = DLM_SGA(y=y[1:(T+S)], S=S, M=100, maxIter=5000, initialM=UpdateM, initialL=UpdateL)
  
  UpdateM = c(MFVBInit$Mu, rep(0, S))
  UpdateL = diag(0.1, T+S+5)
  diag(UpdateL[1:(T+5), 1:(T+5)]) = MFVBInit$Sd
  MFVBUpdate = DLM_SGA(y=y[1:(T+S)], S=S, M=5, maxIter=5000, initialM=UpdateM, initialL=UpdateL, meanfield=TRUE)
  

  # Update initial VB using all y_{T+1:T+S}
  ydensVB = rep(0, 1000)
  ydensMF = rep(0, 1000)
  frSigmaU = FRVBUpdate$L %*% t(FRVBUpdate$L)
  for(i in 1:500){
    draw = rmvnorm(1, c(FRVBUpdate$Mu), frSigmaU)
    ydensVB = ydensVB + dnorm(ysupport, draw[4] + draw[3]*draw[T+S+5], sqrt(exp(draw[1]) + exp(draw[2])))/500
  }
  for(i in 1:500){
    sigyd = exp(rnorm(1, MFVBUpdate$Mu[1], MFVBUpdate$Sd[1]))
    sigxd = exp(rnorm(1, MFVBUpdate$Mu[2], MFVBUpdate$Sd[2]))
    phid = rnorm(1, MFVBUpdate$Mu[3], MFVBUpdate$Sd[3])
    gammad = rnorm(1, MFVBUpdate$Mu[4], MFVBUpdate$Sd[4])
    XTd = rnorm(1, MFVBUpdate$Mu[T+S+5], MFVBUpdate$Sd[T+S+5])
    ydensMF = ydensMF + dnorm(ysupport, gammad + phid*XTd, sqrt(sigyd + sigxd))/500
  }
  # Rerun MCMC using all y up to time T+S
  MCMCupdate = DLM_MCMC(y[1:(T+S)], MCMCrep)
  thetaUpdateKeep = MCMCupdate$theta[(MCMCrep/2 + 1):MCMCrep,]
  xdrawUpdateKeep = MCMCupdate$x[(MCMCrep/2 + 1):MCMCrep,]
  ydens = rep(0, 1000)
  for(i in 1:500){
    u = sample(MCMCrep/2, 1)
    gammaDraw = thetaUpdateKeep[u, 4]
    phiDraw = thetaUpdateKeep[u, 3]
    sigxDraw = thetaUpdateKeep[u, 2]
    sigyDraw = thetaUpdateKeep[u, 1]
    XTS = xdrawUpdateKeep[u, T+S+2:3]
    ydens = ydens + dnorm(ysupport, gammaDraw + phiDraw*XTS[1], sqrt(sigyDraw + phiDraw^2 * XTS[2] + sigxDraw))/500
  }
  obs = min(which((y[T+S+h] < ysupport) == TRUE))
  logscoreRefit = c(log(ydens[obs]), log(ydensMF[obs]), log(ydensVB[obs]))
  meanMCMC = c(meanMCMC, mean(xdrawUpdateKeep[,T+S+1]))
  meanVB = c(meanVB, c(FRVBUpdate$Mu[T+S+5]))
  meanMF = c(meanMF, c(MFVBUpdate$Mu[T+S+5]))

  #MCMCforecast = gather(data.frame(ysupport, ydens, ydens2, ydens3), version, density, -ysupport)
  #MCMCforecast$method = "MCMC"
  #MFforecast = gather(data.frame(ysupport, ydensMF, ydensMF2, ydensMF3), version, density, -ysupport)
  #MFforecast$method = "Diagonal"
  #FRforecast = gather(data.frame(ysupport, ydensVB, ydensVB2, ydensVB3), version, density, -ysupport)
  #FRforecast$method = "Non-Diagonal"
  #forecasts = rbind(MCMCforecast, MFforecast, FRforecast)
  #forecasts$version = rep(rep(c("S+h step", "h step - Filter", "h step - MCMC"), rep(1000, 3)), 3)
  #forecasts$version = factor(forecasts$version, levels=c("S+h step", "h step - Filter", "h step - Refit"))
  return(c(logscoreSh, logscoreFilter, logscoreRefit, meanMCMC, meanVB, meanMF))
}

set.seed(5)
T100S10 = t(replicate(1000, predictiveDensities(100, 10, 1, 15000), simplify = 'matrix'))
MCMCres = data.frame(T100S10[,c(1, 4, 7, 10:15)])
VBres = data.frame(T100S10[,c(2, 5, 8, 16:21)])
MFres = data.frame(T100S10[,c(3, 6, 9, 22:27)])
colnames(MCMCres) = c('S+h step Logscore', 'h step Filtered Logscore', 'h step refit Logscore', 'lnSigYSq', 'lnSigXSq', 'Phi', 'Gamma', 'X_T', 'X_T+S')
colnames(VBres) = c('S+h step Logscore', 'h step Filtered Logscore', 'h step refit Logscore', 'lnSigYSq', 'lnSigXSq', 'Phi', 'Gamma', 'X_T', 'X_T+S')
colnames(MFres) = c('S+h step Logscore', 'h step Filtered Logscore', 'h step refit Logscore', 'lnSigYSq', 'lnSigXSq', 'Phi', 'Gamma', 'X_T', 'X_T+S')
MCMCres$Method = 'MCMC'
VBres$Method = 'Non-Diagonal'
MFres$Method = 'Diagonal'
T50S10df = rbind(MCMCres, VBres, MFres) %>% gather(variable, value, -Method)
ggplot(T50S10df) + geom_boxplot(aes(y=value, x=Method)) + facet_wrap(~variable, scales='free')

p1 <- ggplot(forecasts, aes(x=ysupport, y=density, colour=method)) + facet_wrap(~version) + geom_line() + 
  labs(x=expression(Y[T+S+h]), y="Predictive Density", colour="Approach") +# theme_bw() +
  scale_color_discrete(labels=c('Non-Diagonal VB', 'MCMC', 'Diagonal VB'))
p2 <- ggplot(forecasts, aes(x=ysupport, y=density, colour=version)) + facet_wrap(~method) + geom_line() + 
  labs(x=expression(Y[T+S+h]), y="Predictive Density", colour="Approach") +# theme_bw() +
  scale_color_discrete(labels=c('S+h', 'h Filter', 'h Refit'))
gridExtra::grid.arrange(p1, p2, ncol=1)

drawsXTS = rmvnorm(10000, FRVBUpdate$Mu, frSigmaU)
qplot(drawsXTS[,T+S+5])
qplot(xdrawUpdateKeep[,T+S])

# Variational Bayes Timings
set.seed(1)
VB = list()
ELBOfit = data.frame()

Mvec = c(1, 5, 10, 25, 50, 100)
for(i in 1:300){
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
ggplot(ELBOfit, aes(x=runtime, y=ELBO, colour=factor(M))) + geom_point() + facet_wrap(~Covariance, scales = 'free_x') +
  labs(x="Total Runtime (seconds)", y="Converged ELBO") + scale_color_discrete(name="N")








# Vine Copula
u = sample(1:25000, 1000)
pobtheta = pobs(thetaKeep[u,])
pobx = pobs(xdrawKeep[u,1:(T+1)])

thetavine = RVineStructureSelect(pobtheta)

VineMatrix = matrix(0, T+5, T+5)
diag(VineMatrix)[1:(T+1)] = 5:(T+5)
for(i in 2:(T+1)){
  for(j in 1:i-1){
    VineMatrix[i, j] = T + j - i + 6
  }
}
VineMatrix[(T+2):(T+5),(T+2):(T+5)] = thetavine$Matrix

tree1tau = TauMatrix(cbind(pobtheta, pobx))


VineMatrix[T+5, 1:(T+1)] = which.max(rowSums(abs(tree1tau[1:4, 5:55]))) # 2
VineMatrix[T+2:4, 1:(T+1)] = c(3, 1, 4)


Vine = RVineCopSelect(cbind(pobtheta, pobx), Matrix=VineMatrix, cores=4, indeptest=TRUE, trunclevel=5)

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
