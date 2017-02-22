library(ggplot2)
library(pscl)
library(gridExtra)
dent = function(x, mean, var, df, log = FALSE){
  front = gamma((df + 1)/ 2) / (gamma(1/2) * gamma(df/2))
  middle = (df * var)^(-1/2)
  end = (1 + 1/df * (x - mean)^2/var)^(-(df+1)/2)
  if(log){
    return(log(front) + log(middle) + log(end))
  } else {
    return(front*middle*end)
  }
}

set.seed(15)
mu = 2
sig2 = 1
T = 100
y = rnorm(T, mu, sig2)

gamma = 0
tau = 1
alpha = 1
beta = 1

mubar = (gamma*tau + T * mean(y)) / (tau + T)
vbarbar = T + 2*alpha
sbarsq = (2*beta + sum(y^2) + tau*gamma^2 - (tau + T)*mubar^2)/vbarbar

sigMC = 1/rgamma(100000, alpha + T/2, vbarbar*sbarsq/2)
muMC = rep(0, 10000)
for(i in 1:100000){
  muMC[i] = rnorm(1, (tau * gamma + T * mean(y))/(tau+T), sqrt(sigMC[i] / (tau + T)))
}

sigsup = seq(0.5, 2.2, 0.01)
musup = seq(1.6, 2.6, 0.01)

sigtrue = densigamma(sigsup, alpha + T/2, vbarbar*sbarsq/2)
mutrue = rep(0, length(musup))
for(i in 1:length(musup)){
  mutrue[i] = dent(musup[i], mubar, sbarsq/(T + tau) * vbarbar/(vbarbar-2), vbarbar)
}

tildegamma = (sum(y) + tau*gamma)/(tau + T)
tildealpha = alpha + (T+1)/2

tildes = matrix(0, ncol = 2, nrow = 10)
tildes[1,] = 1
for(i in 2:10){
  tildes[i, 1] = tildes[i-1, 2] / (tildealpha *(T + tau))
  tildes[i, 2] = beta + 1/2 * ((T + tau)*(tildegamma^2 + tildes[i, 1]) - 
                                 2 * tildegamma * (sum(y) + tau*gamma) + sum(y^2) + tau * gamma^2)
}

muvb = dnorm(musup, tildegamma, sqrt(tildes[10, 1]))
sigvb = densigamma(sigsup, tildealpha, tildes[10, 2])

p1 = ggplot() + geom_line(aes(musup, mutrue)) + geom_density(aes(muMC), colour = "red") + 
  geom_line(aes(musup, muvb), colour = "blue") + labs(x = expression(mu), y = "Density") + theme_bw()
p2 = ggplot() + geom_line(aes(sigsup, sigtrue)) + geom_density(aes(sigMC), colour = "red") + 
  geom_line(aes(sigsup, sigvb), colour = "blue") + labs(x = expression(sigma^2), y = "Density") + theme_bw()
grid.arrange(p1, p2, ncol = 2)


