library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(dplyr)
library(ggplot2)
library(GGally)
sourceCpp("kalman.cpp")

#true vales
t.mu = 1.5
t.gamma = 0.8
t.eta2 = 1.5  
t.sigma2 = 1
T = 500

#x0 mean/var
a = 1
b = 1

#generate data
t.x0 = rnorm(1, 1, 1)
t.x = vector(length=T)
y = vector(length=T)
for(t in 1:T) {
  if(t == 1) {
    t.x[t] = 0.5*t.x0 + rnorm(1, 0, sqrt(t.sigma2))
    y[t] = t.mu + t.gamma*t.x[t] + rnorm(1, 0, sqrt(t.eta2))
  } else {
    t.x[t] = 0.5*t.x[(t-1)] + rnorm(1, 0, sqrt(t.sigma2))
    y[t] = t.mu + t.gamma*t.x[t] + rnorm(1, 0, sqrt(t.eta2))
  }
}
sum.y = sum(y)

### MCMC ###

#put draws in here
x.draws = matrix(0, ncol=(T+1), nrow = 500000)  #x0 in col 1, t = 1 to 100 in col 2 to 101
theta.draws = matrix(0, ncol = 4, nrow = 500000) #mu, gamma, eta2, sigma2

#initial values
i.mu = 1
i.gamma = 1
i.eta2 = 1
i.sigma2 = 1

#prior values
mubar = 0
muvar = 10
gammabar = 0
gammavar = 10
etaa = 1              #change these until narrow enough to be stable
etab = 1
siga = 1
sigb = 1

#shape parameter is constant
shape.eta = T/2 + etaa
shape.sigma = T/2 + siga

#Gibbs initial draws
kalman = KalmanLoop(a, b, 0, 0.5, i.sigma2, i.mu, i.gamma, i.eta2, y)
kalman = rbind(c(a, b), kalman)
x.draws[1, ] = backwardsloop(kalman, 0, 0.5, i.sigma2)

#Calculate some stats
sum.x = sum(x.draws[1, 2:(T+1)])
sum.xy = sum(y*x.draws[1, 2:(T+1)])
sum.x2 = sum(x.draws[1, 2:(T+1)]^2)

#Draw Theta conditionals

mean.mu = (muvar*(sum.y-i.gamma*sum.x)+i.eta2*mubar)/(T*muvar+i.eta2)
var.mu = (i.eta2*muvar)/(T*muvar+i.eta2)
theta.draws[1, 1] = rnorm(1, mean.mu, sqrt(var.mu))

mean.gamma = (gammavar*(sum.xy-theta.draws[1, 1]*sum.x)+i.eta2*gammabar)/(sum.x2*gammavar+i.eta2)
var.gamma = (i.eta2*gammavar)/(sum.x2*gammavar+i.eta2)
theta.draws[1, 2] = rnorm(1, mean.gamma, sqrt(var.gamma))

scale.eta = sum((y - (theta.draws[1, 1] + theta.draws[1, 2]*x.draws[1, 2:(T+1)]))^2)/2 + etab
theta.draws[1, 3] = 1/rgamma(1, shape = shape.eta, rate = scale.eta)

scale.sigma = sum((x.draws[1, 2:(T+1)] - 0.5*x.draws[1, 1:T])^2)/2 + sigb
theta.draws[1, 4] = 1/rgamma(1, shape = shape.sigma, rate = scale.sigma)

#Gibbs loop
for(i in 2:500000) {
  kalman = KalmanLoop(a, b, 0, 0.5, theta.draws[(i-1), 4], theta.draws[(i-1), 1], theta.draws[(i-1), 2], theta.draws[(i-1), 3], y)
  kalman = rbind(c(a, b), kalman)
  x.draws[i, ] = backwardsloop(kalman, 0, 0.5, theta.draws[(i-1), 4])
  sum.x = sum(x.draws[i, 2:(T+1)])
  sum.xy = sum(y*x.draws[i, 2:(T+1)])
  sum.x2 = sum(x.draws[i, 2:(T+1)]^2)
  
  mean.mu = (muvar*(sum.y- theta.draws[(i-1), 2]*sum.x)+ theta.draws[(i-1), 3]*mubar)/(T*muvar+ theta.draws[(i-1), 3])
  var.mu = (theta.draws[(i-1), 3]*muvar)/(T*muvar+ theta.draws[(i-1), 3])
  theta.draws[i, 1] = rnorm(1, mean.mu, sqrt(var.mu))
  
  mean.gamma = (gammavar*(sum.xy-theta.draws[i, 1]*sum.x)+theta.draws[(i-1), 3]*gammabar)/(sum.x2*gammavar+theta.draws[(i-1), 3])
  var.gamma = (theta.draws[(i-1), 3]*gammavar)/(sum.x2*gammavar+theta.draws[(i-1), 3])
  theta.draws[i, 2] = rnorm(1, mean.gamma, sqrt(var.gamma))
  
  scale.eta = sum((y - (theta.draws[i, 1] + theta.draws[i, 2]*x.draws[i, 2:(T+1)]))^2)/2 + etab
  theta.draws[i, 3] = 1/rgamma(1, shape = shape.eta, rate = scale.eta)
  
  scale.sigma = sum((x.draws[i, 2:(T+1)] - 0.5*x.draws[i, 1:T])^2)/2 + sigb
  theta.draws[i, 4] = 1/rgamma(1, shape = shape.sigma, rate = scale.sigma)
  
}

colnames(x.draws) = paste0("Xt", 0:T)
colnames(theta.draws) = c("mu", "gamma", "eta_squared", "sigma_squared")

theta.keep = theta.draws[100001:500000,]
x.keep = x.draws[100001:500000,]
#qplot(1:5000, theta.draws[45001:50000,4])
summary(theta.keep)
effectiveSize(theta.keep)

theta.thin = theta.keep[seq(10, 400000, 10),]

#ggpairs(theta.keep)
#qplot(x.keep[, 124], x.keep[, 125])
#qplot(1:25000, theta.keep[,4])


### Mean Field ###

#iterations
N = 75

#hyperparam containers, initialised at 1
mubarbar = c(1, vector(length = (N-1)))
lambda.mu = c(1, vector(length = (N-1)))
gammabarbar = c(1, vector(length = (N-1)))
lambda.gamma = c(1, vector(length = (N-1)))
beta.eta = c(1, vector(length = (N-1)))
beta.sigma = c(1, vector(length = (N-1)))
x0bar = c(1, vector(length = (N-1)))
lambda.x0 = c(1, vector(length = (N-1)))
xtbar = rbind(rep(1, T), matrix(0, ncol = T, nrow = (N-1)))
lambda.xt = rbind(rep(1, T), matrix(0, ncol = T, nrow = (N-1)))

#do not depend on other params
alpha.eta = T/2 + etaa
alpha.sigma = T/2 + siga

#main loop
for(i in 2:N){
  #mubarbar
  numer = muvar*alpha.eta/beta.eta[(i-1)]*(sum.y - gammabarbar[(i-1)]*sum(xtbar[(i-1),])) + mubar
  denom = T*muvar*alpha.eta/beta.eta[(i-1)] + 1
  mubarbar[i] = numer/denom
  
  #lambda.mu
  lambda.mu[i] = muvar/denom
  
  #gammabarbar
  numer = gammavar*alpha.eta/beta.eta[(i-1)]*(sum(xtbar[(i-1),]*y) - sum(xtbar[(i-1),])*mubarbar[i]) + gammabar
  denom = sum(lambda.xt[(i-1),] + xtbar[(i-1),]^2)*gammavar*alpha.eta/beta.eta[(i-1)] + 1
  gammabarbar[i] = numer/denom
  
  #lambda.gamma
  lambda.gamma[i] = gammavar/denom
  
  #beta.eta
  first = sum(y^2) + T*(mubarbar[i]^2 + lambda.mu[i]) + (gammabarbar[i]^2 + lambda.gamma[i])*sum(lambda.xt[(i-1),] + xtbar[(i-1),]^2)
  second = mubarbar[i]*sum.y + gammabarbar[i]*sum(xtbar[(i-1),]*y) + mubarbar[i]*gammabarbar[i]*sum(xtbar[(i-1),])
  beta.eta[i] = first/2 - second + etab
  
  #beta.sigma
  beta.sigma[i] = (sum(lambda.xt[(i-1),] + xtbar[(i-1),]^2) + 1/4*(sum(lambda.xt[(i-1),1:99] + xtbar[(i-1),1:99]^2) + 
                  lambda.x0[(i-1)] + x0bar[(i-1)]) - sum(xtbar[(i-1)]*c(x0bar[(i-1)], xtbar[(i-1), 1:99])))/2 + sigb
  
  #x0bar
  x0bar[i] = 4/5 * (a + 1/2*xtbar[(i-1),1])
  
  #lambda.x0
  lambda.x0[i] = 4/5 * beta.sigma[i]/alpha.sigma
  
  #The other states
  for(j in 1:T) {
    if(j == 1){
      numer = gammabarbar[i]*alpha.eta/beta.eta[i]*(-mubarbar[i]+y[j]) + 1/2 * alpha.sigma/beta.sigma[i] * (x0bar[i] + xtbar[(i-1), (j+1)])
      denom = (gammabarbar[i]^2 + lambda.gamma[i])*alpha.eta/beta.eta[i] + 5/4 * alpha.sigma/beta.sigma[i]
      xtbar[i,j] = numer/denom
      lambda.xt[i,j] = 1/denom
    } else if (j > 1 &  j < T) {
      numer = gammabarbar[i]*alpha.eta/beta.eta[i]*(-mubarbar[i]+y[j]) + 1/2 * alpha.sigma/beta.sigma[i] * (xtbar[i, (j-1)] + xtbar[(i-1), (j+1)])
      denom = (gammabarbar[i]^2 + lambda.gamma[i])*alpha.eta/beta.eta[i] + 5/4 * alpha.sigma/beta.sigma[i]
      xtbar[i,j] = numer/denom
      lambda.xt[i,j] = 1/denom
    } else {
      numer = gammabarbar[i]*alpha.eta/beta.eta[i]*(-mubarbar[i]+y[j]) + 1/2 * alpha.sigma/beta.sigma[i] * xtbar[i, (j-1)]
      denom = (gammabarbar[i]^2 + lambda.gamma[i])*alpha.eta/beta.eta[i] + alpha.sigma/beta.sigma[i]
      xtbar[i,j] = numer/denom
      lambda.xt[i,j] = 1/denom
    }
  }
}
