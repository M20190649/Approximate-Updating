library(ggplot2)
library(mvtnorm)
library(coda)
library(GGally)
library(dplyr)
library(gridExtra)
library(tidyr)
library(pscl)
library(Rcpp)
library(truncnorm)
library(RcppArmadillo)
library(BH)
source("conf_functions.R")
sourceCpp("../AR2.cpp")

set.seed(21)
phi1 = 0.75
phi2 = 0.2
sigma2 = 1
T = 150
J = 150 #for forecasting
y = vector(length = T+J)

gamma0 = sigma2 * (1-phi2) / ((1+phi2)*((1-phi2)^2 - phi1^2))
gamma1 = sigma2 * phi1 / ((1+phi2)*((1-phi2)^2 - phi1^2))
Sigma = matrix(c(gamma0, gamma1, gamma1, gamma0), 2)

y[1:2] = rmvnorm(1, c(0, 0), Sigma)
for(t in 3:(T+J)){
  y[t] = phi1*y[t-1] + phi2*y[t-2] + rnorm(1, sd = sqrt(sigma2))
}

##MCMC
rep = 500000
#T = 30
igShape = T/2

theta = matrix(0, ncol = 3, nrow = rep)
sig2 = 1
phi = c(0.6, 0.1)

tune1 = 0.25
accept1 = 0
tune2 = 0.25
accept2 = 0
set.seed(3)



for(i in 1:rep){
  #Phi 1 Metropolis Hastings
  flag = FALSE
  while(!flag){
    candidate = rnorm(1, phi[1], tune1)
    if(phi[2]  < 1 + candidate & phi[2] < 1 - candidate){
      flag = TRUE
    }
  }
  canpdens = lphidens(y[1:T], candidate, phi[2], sig2)
  oldpdens = lphidens(y[1:T], phi[1], phi[2], sig2)
  canqdens = lqdens(candidate, phi[1], tune1, phi[2], 1)
  oldqdens = lqdens(phi[1], candidate, tune1, phi[2], 1)
  ratio = min(1, exp(canpdens + oldqdens - canqdens - oldpdens))
  if(runif(1) < ratio){
    accept1 = accept1 + 1
    phi[1] = candidate
  }
  
  #Phi 2 Metropolis Hastings
  flag = FALSE
  while(!flag){
    candidate = rnorm(1, phi[2], tune2)
    if(candidate > -1 & candidate  < 1 + phi[1] & candidate < 1 - phi[1]){
      flag = TRUE
    }
  }
  canpdens = lphidens(y[1:T], phi[1], candidate, sig2)
  oldpdens = lphidens(y[1:T], phi[1], phi[2], sig2)
  canqdens = lqdens(candidate, phi[2], tune1, phi[1], 2)
  oldqdens = lqdens(phi[2], candidate, tune1, phi[1], 2)
  ratio = min(1, exp(canpdens + oldqdens - canqdens - oldpdens))
  if(runif(1) < ratio){
    accept2 = accept2 + 1
    phi[2] = candidate
  }
  
  #Sigma-Squared conditional
  rho0 = (1 - phi[2]) / ((1 + phi[2])*((1-phi[2])^2 - phi[1]^2))
  rho1 = phi[1] / ((1 + phi[2])*((1-phi[2])^2 - phi[1]^2))
  igScale = 1/2 * (rho0 * (y[1]^2 + y[2]^2) - 2*rho1*y[1]*y[2]) / (rho0^2 - rho1^2) + 1/2 * sum(((y[3:T] - phi[1]*y[2:(T-1)]- phi[2]*y[1:(T-2)])^2))
  sig2 = 1/rgamma(1, igShape, igScale)
  
  #Store draws
  theta[i, ] = c(phi, sig2)
}
accept1/rep
accept2/rep

stats = apply(theta[5001:500000,], 2, posterior.statistics)
row.names(stats) = c("l95", "mean", "median", "u95")
colnames(stats) = c("Phi1", "Phi2", "SigmaSquared")
stats
colnames(theta) = c("Phi1", "Phi2", "SigmaSquared")
theta %>% as.data.frame() %>% cbind(rep = 1:rep) %>% gather(parameter, draw, -rep) %>% 
  ggplot() + geom_line(aes(rep, draw)) + facet_wrap(~parameter, scales = "free")


thetaKeep = data.frame(theta[seq(50001, 500000, 200),])
effectiveSize(thetaKeep)
ggpairs(thetaKeep)
cov(thetaKeep)

support = seq(min(y)-1, max(y)+1, 0.01)
logscore = vector(length = J)
set.seed(128)
for(t in  (T+1):(T+J)){
  ydens = vector(length = length(support)) 
  for(j in 1:1000){
    index = sample(1:(nrow(thetaKeep)), 1)
    ydens = ydens + dnorm(support, thetaKeep[index, 1]*y[t-1] + thetaKeep[index, 2]*y[t-2], sqrt(thetaKeep[index, 3]))
  }
  ydens = ydens/1000
  if(t == T+1){
    logscore[t-T] = log(ydens[min(which(y[t] < support))])
  } else {
    logscore[t-T] = logscore[t-T-1] + log(ydens[min(which(y[t] < support))])
  }
}
qplot(1:J, logscore, geom = "line")


#VB
library(VineCopula)
pobdraws = pobs(thetaKeep)
Vine = RVineStructureSelect(pobdraws, indeptest = TRUE)
Vine$Matrix
Vine$family
#Gaussian Copula is Optimal
L = t(chol(cov(thetaKeep[, 1:2])))

mu = mean(thetaKeep[,3])
v = var(thetaKeep[,3])
i.a = mu^2/v + 2
i.b = mu*(mu^2/v + 1)
initial = c(mean(thetaKeep[, 1])-0.1, mean(thetaKeep[, 2])-0.1, L[1,1], L[2,1], L[2,2], T/2, i.b)
set.seed(7)

estim = replicate(10, SGA(y[1:150], initial, 1, 0.001, 5000, 0.1, 0.05, 10), simplify = "matrix")
estim
lambda = estim[,which.max(estim[9,])]
lambda

LVB = matrix(c(lambda[3:4], 0, lambda[5]), 2)
Sigma = LVB %*% t(LVB)
Sigma


p1dens = mutate(data.frame(support = seq(0.35, 1.15, length.out = 200)), density = dnorm(support, lambda[1], sqrt(Sigma[1,1])))
p2dens = mutate(data.frame(support = seq(-0.1, 0.55, length.out = 200)), density = dnorm(support, lambda[2], sqrt(Sigma[2,2])))
sig2dens = mutate(data.frame(support = seq(0.501, 2, length.out = 150)), density = densigamma(support, lambda[6], lambda[7]))

plot1 = ggplot() + geom_density(data = thetaKeep, aes(Phi1, y = ..density..)) + theme_bw()
  #geom_line(data = p1dens, aes(support, density), colour = "blue") + theme_bw()
plot2 = ggplot() + geom_density(data = thetaKeep, aes(Phi2, y = ..density..)) + theme_bw()
  #geom_line(data = p2dens, aes(support, density), colour = "blue") + theme_bw()
plot3 = ggplot() + geom_density(data = thetaKeep, aes(SigmaSquared, y = ..density..)) + theme_bw()
 # geom_line(data = sig2dens, aes(support, density), colour = "blue") + theme_bw()

p1p2 = expand.grid(sup1 = seq(0.35, 1.15, length.out = 200), sup2 = seq(-0.1, 0.55, length.out = 200))
p1p2$density = apply(p1p2, 1, mvtnorm::dmvnorm, mean = lambda[1:2], sigma = Sigma) 
p1s2 = expand.grid(sup1 = seq(0.35, 1.15, length.out = 200), sup2 = seq(0.501, 2, length.out = 150))
p1s2$density = apply(p1s2, 1, function(x) dnorm(x[1], lambda[1], sqrt(Sigma[1,1]))*densigamma(x[2], lambda[6], lambda[7]))
p2s2 = expand.grid(sup1 = seq(-0.1, 0.55, length.out = 200), sup2 =seq(0.501, 2, length.out = 150))
p2s2$density = apply(p2s2, 1, function(x) dnorm(x[1], lambda[2], sqrt(Sigma[2,2]))*densigamma(x[2], lambda[6], lambda[7]))

plot4 = ggplot() + geom_density_2d(data=thetaKeep, aes(Phi1, Phi2), colour = "black") +  theme_bw()
  #geom_contour(data=p1p2, aes(sup1, sup2, z = density)) +
plot5 = ggplot() + geom_density_2d(data=thetaKeep, aes(Phi1, SigmaSquared), colour = "black") + theme_bw()
  #geom_contour(data=p1s2, aes(sup1, sup2, z = density)) + theme_bw()
plot6 = ggplot() + geom_density_2d(data=thetaKeep, aes(Phi2, SigmaSquared), colour = "black") + theme_bw()
  #geom_contour(data=p2s2, aes(sup1, sup2, z = density)) 

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)


set.seed(12)
support = seq(min(y)-1, max(y)+1, 0.01)
logscoreVB = vector(length = J)
lambdaUpdate = lambda
for(t in 151:300){
  if(t %% 10 == 0){
    print(paste0("t = ", t))
  }
  ydensVB = vector(length = length(support))
  lambdaUpdate[1:7] = lambdaUpdate[1:7] + c(-0.1, -0.1, 0, 0, 0, 0.5, 0)
  estim = replicate(5, SGA(y[1:t], lambdaUpdate, 1, 0.001, 2000, 0.1, 0.05, 10), simplify = "matrix")
  lambdaUpdate = estim[,which.max(estim[9,])]
  LVBU = matrix(c(lambdaUpdate[3:4], 0, lambdaUpdate[5]), 2)
  SigmaU = LVBU %*% t(LVBU)
  for(j in 1:1000){
    draws = c(rmvnorm(1, lambdaUpdate[1:2], SigmaU), 1/rgamma(1, lambdaUpdate[6], lambdaUpdate[7]))
    ydensVB = ydensVB + dnorm(support, draws[1]*y[t-1] + draws[2]*y[t-2], sqrt(draws[3]))
  }
  ydensVB = ydensVB/1000
  if(t == T+1){
    logscoreVB[t-T] = log(ydensVB[min(which(y[t] < support))])
  } else {
    logscoreVB[t-T] = logscoreVB[t-T-1] + log(ydensVB[min(which(y[t] < support))])
  }
}

df = data.frame(t = 1:J, MCMC = logscore, VB = logscoreVB)
dfl = gather(df, method, score, -t)
ggplot(dfl, aes(t, score, colour = method)) + geom_line() + labs(x = "Extra Data Points", y = "Cumulative Log Score")


denigamma = function(x, shape, scale){ ##integrating constant on pcsl::densigamma was infinite with shape/scale ~ 150
  x ^ (-shape - 1) * exp(-scale/x)
}

p1dens = mutate(p1dens, denUpdate = dnorm(support, lambdaUpdate[1], sqrt(SigmaU[1, 1])))
p2dens = mutate(p2dens, denUpdate = dnorm(support, lambdaUpdate[2], sqrt(SigmaU[2, 2])))
sig2dens = mutate(sig2dens, denUpdate = denigamma(support, lambdaUpdate[6], lambdaUpdate[7]),
                  denUpdate = denUpdate / sum(denUpdate * 0.01))

plot7 = ggplot(data = p1dens) + geom_line(aes(support, denUpdate), colour = "red") + 
  geom_line(aes(support, density), colour = "blue") 
plot8 = ggplot(data = p2dens) + geom_line(aes(support, denUpdate), colour = "red") + 
  geom_line(aes(support, density), colour = "blue") 
plot9 = ggplot(data = sig2dens) + geom_line(aes(support, denUpdate), colour = "red") + 
  geom_line(aes(support, density), colour = "blue") 

p1p2$denUpdate = apply(p1p2, 1, function(x) mvtnorm::dmvnorm(x[1:2], mean = lambdaUpdate[1:2], sigma = SigmaU)) 
p1s2$denUpdate = apply(p1s2, 1, function(x) dnorm(x[1], lambdaUpdate[1], sqrt(SigmaU[1,1]))*denigamma(x[2], lambdaUpdate[6], lambdaUpdate[7]))
p2s2$denUpdate = apply(p2s2, 1, function(x) dnorm(x[1], lambdaUpdate[2], sqrt(SigmaU[2,2]))*denigamma(x[2], lambdaUpdate[6], lambdaUpdate[7]))

plot10 = ggplot(data = p1p2) + geom_contour(aes(sup1, sup2, z = denUpdate), colour = "red") +
  geom_contour(aes(sup1, sup2, z = density), colour = "blue", alpha = 0.5)
plot11 = ggplot(data = p1s2) + geom_contour(aes(sup1, sup2, z = denUpdate), colour = "red") +
  geom_contour(aes(sup1, sup2, z = density), colour = "blue")
plot12 = ggplot(data = p2s2) + geom_contour(aes(sup1, sup2, z = denUpdate), colour = "red") +
  geom_contour(aes(sup1, sup2, z = density), colour = "blue")

grid.arrange(plot7, plot8, plot9, plot10, plot11, plot12, ncol = 3)

draws = array(dim = c(1000, 3, 2))
for(i in 1:1000){
  draws[i,,1] = c(rmvnorm(1, lambda[1:2], Sigma), 1/rgamma(1, lambda[6], lambda[7]))
  draws[i,,2] = c(rmvnorm(1, lambdaUpdate[1:2], SigmaU), 1/rgamma(1, lambdaUpdate[6], lambdaUpdate[7]))
}
colnames(draws) = c("Phi1", "Phi2", "Sig2")
dimnames(draws)[[3]] = c("T=150", "T=300")
draws.l = melt(draws)[,-1]
colnames(draws.l) = c("Param", "Data", "value")
ggplot() + geom_density(data = draws.l, aes(value, colour = Data)) + facet_wrap(~Param, scales = "free")

predictive = data.frame(support = seq(min(y)-1, max(y)+1, length.out = 1000), T150 = 0, T300 = 0)
for(i in 1:1000){
  predictive$T150 = predictive$T150 + dnorm(predictive$support, y[150]*draws[i,1,1] + y[149]*draws[i, 2, 1], sqrt(draws[i, 3, 1]))/1000
  predictive$T300 = predictive$T300 + dnorm(predictive$support, y[150]*draws[i,1,2] + y[149]*draws[i, 2, 2], sqrt(draws[i, 3, 2]))/1000
}
predictive$truth = dnorm(predictive$support, y[150]*0.75 + y[149]*0.2, 1)
predictive.l = gather(predictive, Method, Density, -support)
ggplot(predictive.l) + geom_line(aes(support, Density, colour = Method))

