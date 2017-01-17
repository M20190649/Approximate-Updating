library(ggplot2)
library(mvtnorm)
library(GGally)
library(dplyr)
library(gridExtra)
library(tidyr)
library(pscl)
library(Rcpp)
library(RcppArmadillo)
library(BH)
sourceCpp("../AR2.cpp")

set.seed(82)
phi1 = 0.8
phi2 = -0.3
sigma2 = 2
T = 150
y = vector(length = T)

gamma0 = sigma2 * (1-phi2) / ((1+phi2)*((1-phi2)^2 - phi1^2))
gamma1 = sigma2 * phi1 / ((1+phi2)*((1-phi2)^2 - phi1^2))
Sigma = matrix(c(gamma0, gamma1, gamma1, gamma0), 2)

y[1:2] = rmvnorm(1, c(0, 0), Sigma)
for(t in 3:T){
  y[t] = phi1*y[t-1] + phi2*y[t-2] + rnorm(1, sd = sqrt(sigma2))
}

posterior.statistics = function(x){
  l95 = quantile(x, probs = 0.025)
  mean = mean(x)
  median = median(x)
  u95 = quantile(x, probs = 0.975)
  return(c(l95, mean, median, u95))
}

set.seed(123)
yForecast = c(y, rep(0, 100))
for(t in (T+1):(T+100)){
  #if(t == 180){
  #  phi1 = - 0.3
  #  phi2 = 0.6
  #  sigma2 = 1
  #}
  yForecast[t] = phi1*yForecast[t-1] + phi2*yForecast[t-2] + rnorm(1, sd = sqrt(sigma2))
}

##MCMC
rep = 50000

phi1b = 0
phi2b = 0
phi1l = 3
phi2l = 3
a = 1
b = 1
igShape = T/2 - 1 + a

ym22 = sum(y[1:(T-2)]^2)
y2 = sum(y[3:T]^2)
ym12 = sum(y[2:(T-1)]^2)
yym1 = sum(y[3:T]*y[2:(T-1)])
yym2 = sum(y[3:T]*y[1:(T-2)])
ym1ym2 = sum(y[2:(T-1)]*y[1:(T-2)])

theta = matrix(0, ncol = 3, nrow = rep)
theta[1, 3] = 1
theta[1, 1:2] = 0.3

set.seed(1230)
for(i in 2:rep){
  #phi 1 conditional
  mean1 = (phi1l*(yym1 - theta[i-1, 2]*ym1ym2) + theta[i-1, 3]*phi1b)/(phi1l*ym12 + theta[i-1, 3])
  var1 = theta[i-1, 3]*phi1l/(phi1l*ym12+theta[i-1, 3]) 
  theta[i, 1] = rnorm(1, mean1, sqrt(var1))
  
  #phi2 conditional
  mean2 = (phi2l*(yym2 - theta[i, 1]*ym1ym2) + theta[i-1, 3]*phi2b)/(phi2l*ym22 + theta[i-1, 3])
  var2 = theta[i-1, 3]*phi2l/(phi2l*ym22+theta[i-1, 3]) 
  theta[i, 2] = rnorm(1, mean2, sqrt(var2))
  
  #sigma2 conditional
  igScale = b + sum(((y[3:T] - theta[i, 1]*y[2:(T-1)]- theta[i, 2]*y[1:(T-2)])^2))/2
  theta[i, 3] = 1/rgamma(1, igShape, igScale)
}

stats = apply(theta[(0.5*rep + 1) : rep,], 2, posterior.statistics)
row.names(stats) = c("l95", "mean", "median", "u95")
colnames(stats) = c("Phi1", "Phi2", "SigmaSquared")
stats
theta %>% as.data.frame() %>% cbind(rep = 1:rep) %>% gather(parameter, draw, -rep) %>% 
  ggplot() + geom_line(aes(rep, draw)) + facet_wrap(~parameter, scales = "free")

thetaKeep = data.frame(theta[(0.5*rep+1):rep,])
colnames(theta) = c("Phi1", "Phi2", "SigmaSquared")
ggpairs(thetaKeep)
cov(thetaKeep)

support = seq(min(y)-1, max(y)+1, 0.01)
logscore = vector(length = 100)
set.seed(128)
for(t in  (T+1):(T+100)){
  ydens = vector(length = length(support)) 
  for(j in 1:1000){
    index = sample(1:(0.5*rep), 3, replace = TRUE)
    ydens = ydens + dnorm(support, thetaKeep[index[1], 1]*yForecast[t-1] + thetaKeep[index[2], 2]*yForecast[t-2], sqrt(thetaKeep[index[3], 3]))
  }
  ydens = ydens/1000
  if(t == T+1){
    logscore[t-T] = log(ydens[min(which(yForecast[t] < support))])
  } else {
    logscore[t-T] = logscore[t-T-1] + log(ydens[min(which(yForecast[t] < support))])
  }
}
qplot(1:100, logscore, geom = "line")


#VB
library(VineCopula)
pobdraws = pobs(thetaKeep)
Vine = RVineStructureSelect(pobdraws, indeptest = TRUE)
Vine$Matrix
Vine$family
#BVN + IG is Optimal
ELBOr = function(lambda, y, n = 1000){
  epsilon = simc(n)
  theta = rep(0, 3)
  out = rep(0, n)
  for(i in 1:n){
    theta[1] = lambda[1] + lambda[3]*epsilon[i, 1]
    theta[2] = lambda[2] + lambda[4]*epsilon[i, 1] + lambda[5]*epsilon[i, 2]
    theta[3] = qigammac(epsilon[i, 3], lambda[6], lambda[7])
    out[i] = logjointc(y, theta[1], theta[2], theta[3]) - logqc(theta, lambda)
  }
  return(mean(out))
}

L = t(chol(cov(thetaKeep[, 1:2])))

mu = mean(thetaKeep[,3])
v = var(thetaKeep[,3])
i.a = mu^2/v + 2
i.b = mu*(mu^2/v + 1)
initial = c(mean(thetaKeep[, 1]), mean(thetaKeep[, 2]), L[1,1], L[2,1], L[2,2], i.a, i.b)
set.seed(11)


estim = replicate(10, SGA(y, initial, 1, 0.001, 5000, 0.25, ELBOr), simplify = "matrix")
lambda = estim[,which.max(estim[9,])]

LVB = matrix(c(lambda[3:4], 0, lambda[5]), 2)
Sigma = LVB %*% t(LVB)

p1dens = mutate(data.frame(support = seq(0.2, 1, length.out = 200)), density = dnorm(support, lambda[1], sqrt(Sigma[1,1])))
p2dens = mutate(data.frame(support = seq(-0.6, 0.2, length.out = 200)), density = dnorm(support, lambda[2], sqrt(Sigma[2,2])))
sig2dens = mutate(data.frame(support = seq(0.5, 3, length.out = 200)), density = densigamma(support, lambda[6], lambda[7]))


plot1 = ggplot() + geom_density(data = thetaKeep, aes(Phi1, y = ..density..), colour = "red") + 
  geom_line(data = p1dens, aes(support, density), colour = "blue") 
plot2 = ggplot() + geom_density(data = thetaKeep, aes(Phi2, y = ..density..), colour = "red") + 
  geom_line(data = p2dens, aes(support, density), colour = "blue") 
plot3 = ggplot() + geom_density(data = thetaKeep, aes(SigmaSquared, y = ..density..), colour = "red") + 
  geom_line(data = sig2dens, aes(support, density), colour = "blue") 

p1p2 = expand.grid(sup1 = seq(0.2, 1, length.out = 200), sup2 = seq(-0.6, 0.2, length.out = 200))
p1p2$density = apply(p1p2, 1, dmvnorm, mean = lambda[1:2], sigma = Sigma) 
p1s2 = expand.grid(sup1 = seq(0.2, 1, length.out = 200), sup2 = seq(0.5, 3, length.out = 200))
p1s2$density = apply(p1s2, 1, function(x) dnorm(x[1], lambda[1], sqrt(Sigma[1,1]))*densigamma(x[2], lambda[6], lambda[7]))
p2s2 = expand.grid(sup1 = seq(-0.6, 0.2, length.out = 200), sup2 =seq(0.5, 3, length.out = 200))
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
logscoreVB = vector(length = 100)
lambda = initial
for(t in (T+1):(T+100)){
  ydens = vector(length = length(support))
  estim = replicate(5, SGA(yForecast[1:(t-1)], lambda[1:7], 1, 0.001, 5000, 0.1), simplify = "matrix")
  lambda = estim[,which.max(estim[9,])]
  LVB = matrix(c(lambda[3:4], 0, lambda[5]), 2)
  Sigma = LVB %*% t(LVB)
  for(j in 1:1000){
    draws = c(rmvnorm(1, lambda[1:2], Sigma), 1/rgamma(1, lambda[6], lambda[7]))
    ydens = ydens + dnorm(support, draws[1]*yForecast[t-1] + draws[2]*yForecast[t-2], sqrt(draws[3]))
  }
  ydens = ydens/1000
  if(t == T+1){
    logscoreVB[t-T] = log(ydens[min(which(yForecast[t] < support))])
  } else {
    logscoreVB[t-T] = logscoreVB[t-T-1] + log(ydens[min(which(yForecast[t] < support))])
  }
}

