library(tidyr)
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(dplyr)
library(ggplot2)
sourceCpp("heirarchical.cpp")

#Drawing data
t.mubar = 3.5
t.tau2 = 1
J = 100
N = 1000
t.muj = rnorm(J, t.mubar, sqrt(t.tau2))
y = sapply(t.muj, rnorm, n = N, sd = 1)
ybar = apply(y, 2, mean)

#Hyperpriors
mubarbar = 2
lambda = 10
a = 0.0001
b = 0.0001

#Gibbs setup
M = 5000 #retained
B = 5000 #burn in
draws = matrix(0, ncol=J+2, nrow=M+B) #1-100: muj, 101 = mubar, 102 = tau2
#Initials
draws[1,] = 1

#Gibbs Chain
draws = gibbschain(y, M, B, mubarbar, lambda, a, b)
results = mcmc(draws[(B+1):(B+M),])
qplot(draws[(B+1):(M+B), 15], draws[(B+1):(M+B), 102]) + labs(x="Mu bar", y="Tau squared")

#Variational setup
n = 10
cv = varfit(y, n, mubarbar, lambda, a, b)
colnames(cv) = c(1:100, "lstar", "mubar", "lambdabar", "alpha", "beta")
cv = as.data.frame(cv)

#Plotting results
support = seq(2.9, 4.1, 0.01)
vdens = dnorm(support, cv$mubar[n], sqrt(cv$lambdabar[n]))
mcdens = dnorm(support, mean(draws[(B+1):(M+B), 101]), sd(draws[(B+1):(M+B), 101]))
p.mubar = data.frame(variational = vdens, mcmc = mcdens, support = support)
p.mubar = gather(p.mubar, key=method, value=density, -support)
ggplot(data=p.mubar, aes(y=density, x=support, colour=method))+geom_line()


#Draw new data
muj.new = rnorm(10000, t.mubar, sqrt(t.tau2))
ynew = sapply(muj.new, rnorm, n = 1, sd = 1)

#Calculate log score as an expectation over the posterior
samp = sample((B+1):(M+B), 1000)
theta.s = draws[samp, 101:102]

#To use apply with a vector argument
vec.rnorm = function(theta, n = 1){
  rnorm(n, theta[1], sqrt(theta[2]))
}

#Draw muj|theta.s
mu.draw = apply(theta.s, 1, vec.rnorm)

score = vector(length=10000)

for(i in 1:10000){
  score[i] = log(sum(sapply(mu.draw, dnorm, x=ynew[i], sd = 1))/1000)
}

#Log scores
theta.v = cbind(rnorm(1000, vmubar.star[n], sqrt(vlbar.star[n])),
                1/rgamma(1000, shape = valpha, rate=vbeta[n]))

mu.draw.v = apply(theta.v, 1, vec.rnorm)
score.v = vector(length=10000)

for(i in 1:10000){
  score.v[i] = log(sum(sapply(mu.draw.v, dnorm, x=ynew[i], sd = 1))/1000)
}

scores = data.frame(mcmc = score, variational=score.v, data = ynew)
scoresl = gather(scores, key=method, value=log.score, -data)
ggplot(data=scoresl, aes(x=data, y=log.score, colour=method)) + geom_line()
