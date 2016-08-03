##Generic regression
library(tidyr)
library(ggplot2)
library(dplyr)
library(MASS)

T <- 250
sigma <- 0.75
b1 <- 1.3
b2 <- 2.5
eps <- rnorm(T)
x1 <- rnorm(T)
x2 <- rnorm(T)
y <- b1*x1+b2*x2+sigma*eps

X <- cbind(x1, x2)
XX <- t(X)%*%X
b <- solve(XX)%*%t(X)%*%y

#Variational

k <- 2
sigsqhat <- t(y-X%*%b)%*%(y-X%*%b)/(T-k)

scale <- rep(0, 25)
scale[1] <- (T-k)*sigsqhat/2
for(i in 2:25){
  scale[i] <- (T-k)*sigsqhat/2 + scale[(i-1)]/T
}
betavar2 <- mvrnorm(2000, b, solve(t(X)%*%X)*scale[25]*2/T)
sigmavar2 <- 1/(rgamma(2000, T/2, scale[25]))
variational <- data.frame(b1 = betavar2[,1], b2 = betavar2[,2], sigma.squared = sigmavar2, method="Variational")

#Exact

sigmaexact <- 1/rgamma(2000, (T-k)/2, sigsqhat*(T-k)/2)
sigtransfex <- array(sapply(sigmaexact, function(x) x*solve(XX)), dim=c(2,2,2000))
betaexact <- t(apply(sigtransfex, 3, mvrnorm, n=1, mu=b))

exact <- data.frame(b1 = betaexact[,1], b2 = betaexact[,2], sigma.squared = sigmaexact, method="Exact")

combined <- rbind(exact, variational)
combl <- gather(combined, theta, density, -method)
ggplot(data=combl, aes(x=density, colour=method)) + geom_density() + facet_wrap(~theta, scales="free")


