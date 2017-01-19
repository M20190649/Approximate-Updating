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
phi1 = 0.78
phi2 = 0.2
sigma2 = 2
T = 150
J = 100 #for forecasting
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

igShape = T/2

theta = matrix(0, ncol = 3, nrow = rep)
sig2 = 1
phi = c(0.6, 0.1)

tune1 = 0.15
accept1 = 0
tune2 = 0.15
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

stats = apply(theta[5001:50000,], 2, posterior.statistics)
row.names(stats) = c("l95", "mean", "median", "u95")
colnames(stats) = c("Phi1", "Phi2", "SigmaSquared")
stats
colnames(theta) = c("Phi1", "Phi2", "SigmaSquared")
theta %>% as.data.frame() %>% cbind(rep = 1:rep) %>% gather(parameter, draw, -rep) %>% 
  ggplot() + geom_line(aes(rep, draw)) + facet_wrap(~parameter, scales = "free")


thetaKeep = data.frame(theta[seq(5001, 500000, 600),])
effectiveSize(thetaKeep)
ggpairs(thetaKeep)
cov(thetaKeep)

support = seq(min(y)-1, max(y)+1, 0.01)
logscore = vector(length = J)
set.seed(128)
for(t in  (T+1):(T+J)){
  ydens = vector(length = length(support)) 
  for(j in 1:1000){
    index = sample(1:(nrow(thetaKeep)), 3, replace = TRUE)
    ydens = ydens + dnorm(support, thetaKeep[index[1], 1]*y[t-1] + thetaKeep[index[2], 2]*y[t-2], sqrt(thetaKeep[index[3], 3]))
  }
  ydens = ydens/1000
  if(t == T+1){
    logscore[t-T] = log(ydens[min(which(y[t] < support))])
  } else {
    logscore[t-T] = logscore[t-T-1] + log(ydens[min(which(y[t] < support))])
  }
}
qplot(1:100, logscore, geom = "line")


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
initial = c(mean(thetaKeep[, 1])-0.1, mean(thetaKeep[, 2])-0.1, L[1,1], L[2,1], L[2,2], i.a, i.b)
set.seed(18)

estim = replicate(10, SGA(y[1:T], initial, 5, 0.001, 5000, 0.075, 0.02, 0.1), simplify = "matrix")
estim
lambda = estim[,which.max(estim[9,])]
lambda

LVB = matrix(c(lambda[3:4], 0, lambda[5]), 2)
Sigma = LVB %*% t(LVB)
Sigma

p1dens = mutate(data.frame(support = seq(0.2, 1.2, length.out = 200)), density = dnorm(support, lambda[1], sqrt(Sigma[1,1])))
p2dens = mutate(data.frame(support = seq(-0.3, 0.7, length.out = 200)), density = dnorm(support, lambda[2], sqrt(Sigma[2,2])))
sig2dens = mutate(data.frame(support = seq(0.5, 3.5, length.out = 200)), density = densigamma(support, lambda[6], lambda[7]))


plot1 = ggplot() + geom_density(data = thetaKeep, aes(Phi1, y = ..density..), colour = "red") + 
  geom_line(data = p1dens, aes(support, density), colour = "blue") 
plot2 = ggplot() + geom_density(data = thetaKeep, aes(Phi2, y = ..density..), colour = "red") + 
  geom_line(data = p2dens, aes(support, density), colour = "blue") 
plot3 = ggplot() + geom_density(data = thetaKeep, aes(SigmaSquared, y = ..density..), colour = "red") + 
  geom_line(data = sig2dens, aes(support, density), colour = "blue") 

p1p2 = expand.grid(sup1 = seq(0.2, 1.2, length.out = 200), sup2 = seq(-0.3, 0.7, length.out = 200))
p1p2$density = apply(p1p2, 1, mvtnorm::dmvnorm, mean = lambda[1:2], sigma = Sigma) 
p1s2 = expand.grid(sup1 = seq(0.2, 1.2, length.out = 200), sup2 = seq(0.5, 3.5, length.out = 200))
p1s2$density = apply(p1s2, 1, function(x) dnorm(x[1], lambda[1], sqrt(Sigma[1,1]))*densigamma(x[2], lambda[6], lambda[7]))
p2s2 = expand.grid(sup1 = seq(-0.3, 0.7, length.out = 200), sup2 =seq(0.5, 3, length.out = 200))
p2s2$density = apply(p2s2, 1, function(x) dnorm(x[1], lambda[2], sqrt(Sigma[2,2]))*densigamma(x[2], lambda[6], lambda[7]))

plot4 = ggplot() + geom_density_2d(data=thetaKeep, aes(Phi1, Phi2), colour = "red") + 
  geom_contour(data=p1p2, aes(sup1, sup2, z = density))
plot5 = ggplot() + geom_density_2d(data=thetaKeep, aes(Phi1, SigmaSquared), colour = "red") + 
  geom_contour(data=p1s2, aes(sup1, sup2, z = density))
plot6 = ggplot() + geom_density_2d(data=thetaKeep, aes(Phi2, SigmaSquared), colour = "red") + 
  geom_contour(data=p2s2, aes(sup1, sup2, z = density))

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)

set.seed(15)
support = seq(min(y)-1, max(y)+1, 0.01)
logscoreVB = vector(length = J)
lambdaUpdate = initial
for(t in (T+1):(T+J)){
  ydensVB = vector(length = length(support))
  if(t > T+1){
    lambdaUpdate[1:2] = lambdaUpdate[1:2] - 0.1
  }
  estim = replicate(3, SGA(y[1:t], lambdaUpdate[1:7], 5, 0.01, 5000, 0.075, 0.02, 0.1), simplify = "matrix")
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


sup = seq(0.3, 1.2, 0.01)
den1 = dnorm(sup, lambda[1], sqrt(Sigma[1, 1]))
den2 = dnorm(sup, lambdaUpdate[1], sqrt(SigmaU[1, 1]))
ggplot() + geom_line(aes(sup, den1), colour = "red") + geom_line(aes(sup, den2), colour = "blue")
