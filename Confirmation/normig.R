library(ggplot2)
library(pscl)

set.seed(15)
mu = 2
sig2 = 1
T = 100
y = rnorm(T, mu, sig2)

gamma = 0
tau = 10
alpha = 1
beta = 1

mubar = (gamma*tau + T * mean(y)) / (tau + T)
vbarbar = T + 2*alpha
sbarsq = (2*beta + sum(y^2) + tau*gamma^2 - (tau + T)*mubarbar^2)/vbarbar

sigMC = 1/rgamma(100000, alpha + T/2, vbarbar*sbarsq/2)
muMC = rep(0, 10000)
for(i in 1:100000){
  muMC[i] = rnorm(1, (tau * gamma + T * mean(y))/(tau+T), sqrt(sigMC[i] / (tau + T)))
}

sigsup = seq(0.8, 2.5, 0.01)
musup = seq(1.4, 2.4, 0.01)

sigtrue = densigamma(sigsup, alpha + T/2, vbarbar*sbarsq/2)
mutrue = rep(0, length(musup))
for(i in 1:length(musup)){
  mutrue[i] = dent(musup[i], mubarbar, sbarsq/(T + tau) * vbarbar/(vbarbar-2), vbarbar)
}
