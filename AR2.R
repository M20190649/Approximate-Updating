##AR(2) Model
library(tidyr)
library(ggplot2)
library(dplyr)
library(MASS)

T <- 50
J <- 10
sigma <- 0.5
mu <- 2
psi1 <- 0.4
psi2 <- 0.4
eps <- rnorm(T+J, 0, sigma)

y <- mu + psi1*(0-mu) + psi2*(0-mu) + eps[1]
y[2] <- mu + psi1*(y[1]-mu) + psi2*(0-mu) + eps[2]
for(t in 3:60){
  y[t] <- mu + psi1*(y[(t-1)]-mu) + psi2*(y[(t-2)]-mu) + eps[t]
}
X <- as.matrix(cbind(rep(1,58), y[2:59], y[1:58]))
XX <- t(X)%*%X
b <- solve(XX)%*%t(X)%*%y[3:60]                   
v <- T+J-5

reps <- 2000
dimbeta <- 3

#exact posterior
sigsqexact <- t(y[3:60]-X%*%b)%*%(y[3:60]-X%*%b)/v
sigmaexact <- 1/sqrt(rgamma(reps, v/2, sigsqexact*v/2))
sigtransfex <- array(sapply(sigmaexact, function(x) x*solve(XX)), dim=c(dimbeta,dimbeta,reps))
betaexact <- t(apply(sigtransfex, 3, mvrnorm, n=1, mu=b))

exact <- data.frame(mu = betaexact[,1], psi1 = betaexact[,2], psi2 = betaexact[,3], sigma = sigmaexact, method="Exact")

#setting up the variational algorithm
sigsqhat <- var(y[3:60])
covar <- array(0, dim=c(dimbeta,dimbeta,10))
sigsq <- rep(0, 10)

transf <- function(x){
  out <- t(x-b)%*%XX%*%(x-b)
  return(out)
}

#the iteration
for(i in 1:10){
  if(i ==1){
    taudraw <- rgamma(10000, (T+J-2)/2, (T+J-2)*sigsqhat/2)
    precision <- XX*mean(taudraw)
    covar[,,1] <- solve(precision)
    betadraw <- mvrnorm(10000, b, covar[,,1])
    transform <- mean(apply(betadraw, 1, transf))
    sigsq[1] <- (v*sigsqhat+transform)/(T+J-2)
  } else {
    taudraw <- rgamma(10000, (T+J-2)/2, (T+J-2)*sigsq[(i-1)]/2)
    precision <- XX*mean(taudraw)
    covar[,,i] <- solve(precision)
    betadraw <- mvrnorm(10000, b, covar[,,i])
    transform <- mean(apply(betadraw, 1, transf))
    sigsq[i] <- (v*sigsqhat+transform)/(T+J-2)
  }
}
betavar <- mvrnorm(2000, b, covar[,,10])
sigmavar <- 1/sqrt(rgamma(2000, (T+J-2)/2, (T+J-2)*sigsq[10]/2))
variational <- data.frame(mu = betavar[,1], psi1 = betavar[,2], psi2 = betavar[,3], sigma = sigmavar, method="Variational")


#Other variational bayes algorithm
k <- 3
v <- T+J-2
sigsqhat <- t(y[3:60]-X%*%b)%*%(y[3:60]-X%*%b)/(v-k)

scale <- rep(0, 25)
scale[1] <- (v-k)*sigsqhat/2
for(i in 2:25){
  scale[i] <- (T-k)*sigsqhat/2 + scale[(i-1)]/T
}
betavar2 <- mvrnorm(2000, b, solve(t(X)%*%X)*scale[25]*2/v)
sigmavar2 <- 1/sqrt(rgamma(2000, v/2, scale[25]))
var2 <- data.frame(mu = betavar2[,1], psi1 = betavar2[,2], psi2 = betavar2[,3], sigma = sigmavar2, method="Variational 2")


#The ABC Algorithm

sample.stats <- function(x) {
  out <- c(mean(x), cor(x[1:9], x[2:10]), cor(x[1:8], x[3:10]), var(x))
  return(out)
}

index.bottom.N = function(x, N){
  o = order(x, na.last=FALSE, decreasing=TRUE)
  o.length = length(o)
  o[((o.length-N+1):o.length)]
}

distance <- function(x) {
  out <- sqrt(sum((yss-x)^2))
  return(out)
}

yss <- sample.stats(y[51:60])
XXprior <- t(X[1:48,])%*%X[1:48,]
bprior <- solve(XXprior)%*%t(X[1:48,])%*%y[3:50]
sigsqprior <- as.numeric(t(y[3:50]-X[1:48,]%*%bprior)%*%(y[3:50]-X[1:48,]%*%bprior)/(T-5))
sigmaprior <- 1/sqrt(rgamma(200000, (T-5)/2, sigsqprior*(T-5)/2))
sigtransf <- array(sapply(sigmaprior, function(x) x*solve(XXprior)), dim=c(3,3,200000))
betaprior <- t(apply(sigtransf, 3, mvrnorm, n=1, mu=bprior))
theta <- cbind(betaprior, sigmaprior)
z <- matrix(0, nrow=10, ncol=200000)
z[1,] <- y[51]
z[2,] <- y[52]
eps <- sapply(theta[,4], rnorm, n=8, mean=0)

for(t in 3:10){
  z[t,] <- theta[,1] + theta[,2]*(z[(t-1)]-theta[,1]) + theta[,3]*(z[(t-2)]-theta[,1]) + eps[(t-2),]
}
zss <- apply(z, 2, sample.stats)
dist <- apply(zss, 2, distance)
accept <- theta[index.bottom.N(dist, 2000),]
abc <- data.frame(mu = accept[,1], psi1 = accept[,2], psi2 = accept[,3], sigma = accept[,4], method = "ABC")

##the prior
prior <- data.frame(mu = theta[1:2000,1], psi1 = theta[1:2000,2], psi2 = theta[1:2000,3], sigma = theta[1:2000,4], method = "Prior")

AR2 <- rbind(prior, exact, variational, abc)
AR2 <- filter(AR2, psi1 + psi2 < 1 & psi2 - psi1 < 1 & abs(psi2) < 1)
#AR2 <- mutate(AR2, mu = mu/(1-psi1-psi2))
AR2l <- gather(AR2, theta, density, -method)

ggplot(data=AR22l, aes(x=density, colour=method)) + geom_density() + facet_wrap(~theta, scales="free")

AR22 <- rbind(exact, variational, var2)
AR22l <- gather(AR22, theta, density, -method)


