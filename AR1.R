##AR(1) Model
library(Rcpp)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RcppArmadillo)
sourceCpp("AR1c.cpp")

T <- 50
J <- 10
sigma <- 0.5
mu <- 2
rho <- 0.7
eps <- rnorm(T+J, 0, sigma)

#y <- genAR1(eps, rho, mu)

y <- mu + rho*(0-mu) + eps[1]
for(t in 2:60){
  y[t] <- mu + rho*(y[(t-1)]-mu) + eps[t]
}
X <- as.matrix(cbind(rep(1,59), y[1:59]))
#XX <- mult(t(X), X)
XX <- t(X)%*%X
b <- solve(XX)%*%t(X)%*%y[2:60]                   
v <- T+J-3

#exact posterior
sigsqexact <- t(y[2:60]-X%*%b)%*%(y[2:60]-X%*%b)/(T+J-3)
sigmaexact <- 1/sqrt(rgamma(2000, v/2, sigsqexact*v/2))
sigtransfex <- array(sapply(sigmaexact, function(x) x*solve(XX)), dim=c(2,2,2000))
betaexact <- t(apply(sigtransfex, 3, mvrnorm, n=1, mu=b))
exact <- data.frame(mu = betaexact[,1], psi = betaexact[,2], sigma = sigmaexact, method="Exact")

#setting up the variational algorithm
sigsqhat <- var(y[2:60])
covar <- array(0, dim=c(2,2,10))
sigsq <- rep(0, 10)

transf <- function(x){
  out <- t(x-b)%*%XX%*%(x-b)
  return(out)
}

#the iteration
for(i in 1:10){
  if(i ==1){
    taudraw <- rgamma(10000, (T+J-1)/2, (T+J-1)*sigsqhat/2)
    precision <- XX*mean(taudraw)
    covar[,,1] <- solve(precision)
    betadraw <- mvrnorm(10000, b, covar[,,1])
    transform <- mean(apply(betadraw, 1, transf))
    sigsq[1] <- ((T+J-3)*sigsqhat+transform)/(T+J-1)
  } else {
    taudraw <- rgamma(10000, (T+J-1)/2, (T+J-1)*sigsq[(i-1)]/2)
    precision <- XX*mean(taudraw)
    covar[,,i] <- solve(precision)
    betadraw <- mvrnorm(10000, b, covar[,,i])
    transform <- mean(apply(betadraw, 1, transf))
    sigsq[i] <- ((T+J-3)*sigsqhat+transform)/(T+J-1)
  }
}
betavar <- mvrnorm(2000, b, covar[,,10])
sigmavar <- 1/sqrt(rgamma(2000, (T+J-1)/2, (T+J-1)*sigsq[10]/2))
variational <- data.frame(mu = betavar[,1], psi = betavar[,2], sigma = sigmavar, method="Variational")

#The ABC Algorithm

sample.stats <- function(x) {
  out <- c(mean(x), cor(x[1:9], x[2:10]), var(x))
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

ABCR <- function(){
yss <- sample.stats(y[51:60])
XXprior <- t(X[1:49,])%*%X[1:49,]
bprior <- solve(XXprior)%*%t(X[1:49,])%*%y[2:50]
sigsqprior <- as.numeric(t(y[2:50]-X[1:49,]%*%bprior)%*%(y[2:50]-X[1:49,]%*%bprior)/(T-3))
sigmaprior <- 1/sqrt(rgamma(200000, (T-3)/2, sigsqprior*(T-3)/2))
sigtransf <- array(sapply(sigmaprior, function(x) x*solve(XXprior)), dim=c(2,2,200000))
betaprior <- t(apply(sigtransf, 3, mvrnorm, n=1, mu=bprior))
theta <- cbind(betaprior, sigmaprior)
z <- matrix(0, nrow=10, ncol=200000)
z[1,] <- y[51]
eps <- sapply(theta[,3], rnorm, n=9, mean=0)

for(t in 2:10){
  z[t,] <- theta[,1] + theta[,2]*(z[(t-1)]-theta[,1]) + eps[(t-1),]
}
zss <- apply(z, 2, sample.stats)
dist <- apply(zss, 2, distance)
accept <- theta[index.bottom.N(dist, 2000),]
abc <- data.frame(mu = accept[,1], psi = accept[,2], sigma = accept[,3], method = "ABC")
}
##the prior
prior <- data.frame(mu = theta[1:2000,1], psi = theta[1:2000,2], sigma = theta[1:2000,3], method = "Prior")

AR1 <- rbind(prior, exact, variational, abc)
AR1 <- mutate(AR1, mu.actual = mu/(1-psi))
AR1l <- gather(AR1, theta, density, -method)

ggplot(data=AR1l, aes(x=density, colour=method)) + geom_density() + facet_wrap(~theta, scales="free")

##One step variational forecast
T <- 100
tmu <- 3
tpsi <- 0.5
tsig <- 1
y <- genAR1(T, c(tmu, tpsi, tsig))

#Empty vectors
N <- 1000
barmu = vector(length=N)
lambda = vector(length=N)
hatpsi = vector(length=N)
gamma = vector(length=N)
a = vector(length=N)
b = vector(length=N)
yhat = vector(length=N)
d = vector(length=N)

b0 = 1
a0 = 1
yhat0 = 1
psi0 = 0.8
gamma0 = 1
barmu0 = 1
lambda0 = 1
d0 = 1
p1 = 2
p2 = 3

#psi kernel function
psi.k <- function(psi, y, a, b, yhat, mu, p1, p2) {
  psibar = (sum(y[2:T]*y[1:T-1])+yhat*y[T]-mu*sum(y))/sum(y^2)
  out = sum(y^2)*(-b/a)*(psi-psibar)^2+(p1-1)*log((1+psi)/2)+(p2-1)*log((1-psi)/2)
  return(out)
}

psi.hessian <- function(psi, y, a, b, p1, p2) {
  const = sum(y^2)*(-b/a)*2
  out = const-(p1-1)/(psi+1)^2-(p2-1)/(1-psi)^2
  return(out)
}

#Iteration 1
barmu[1] = 1/T*(sum(y[2:T]) + yhat0 - psi0*sum(y))
lambda[1] = b0/(a0*T)
hatpsi[1] = optimise(psi.k, c(-0.99, 0.99), y, a0, b0, yhat0, barmu[1], p1, p2, maximum=TRUE)$maximum
gamma[1] = -1/(psi.hessian(hatpsi[1], y, a0, b0, p1, p2))
a = (T-2)/2
b[1] = 1/2*(sum(y[2:T]^2-2*y[2:T]*(barmu[1]+hatpsi[1]*y[1:(T-1)])) + T*(barmu[1]^2+lambda[1]) +
            sum((hatpsi[1]^2+gamma[1])*y[1:T]^2+2*barmu[1]*hatpsi[1]*y[1:T]) + yhat0^2+d0-2*yhat0*(barmu[1]+hatpsi[1]*y[T]))
yhat[1] = barmu[1] + hatpsi[1]*y[T]
d[1] = b[1]/a

#Loop
for(i in 2:N){
  barmu[i] = 1/T*(sum(y[2:T]) + yhat[(i-1)] - hatpsi[(i-1)]*sum(y))
  lambda[i] = b[(i-1)]/(a*T)
  hatpsi[i] = optimise(psi.k, c(0, 1), y, a, b[(i-1)], yhat[(i-1)], barmu[i], p1, p2, maximum=TRUE)$maximum
  gamma[i] = -1/(psi.hessian(hatpsi[i], y, a, b[(i-1)], p1, p2))
  b[i] = 1/2*(sum(y[2:T]^2-2*y[2:T]*(barmu[i]+hatpsi[i]*y[1:(T-1)])) + T*(barmu[i]^2+lambda[i]) +
                sum((hatpsi[i]^2+gamma[i])*y[1:T]^2+2*barmu[i]*hatpsi[i]*y[1:T]) +
                yhat[(i-1)]^2+d[(i-1)]-2*yhat[(i-1)]*(barmu[i]+hatpsi[i]*y[T]))
  yhat[i] = barmu[i] + hatpsi[i]*y[T]
  d[i] = b[i]/a
}

support <- seq(0, 6, 0.005)
dens <- dnorm(support, yhat[1000], sqrt(d[1000]))
qplot(support, dens, geom="line")

##MCMC
B <- 1000
K <- 50000
chain <- matrix(0, ncol=3, nrow=(B+K))
mu0 = 1
sig0 = 1
psi0 = 0.8

psi.kmcmc <- function(psi, mu, sigma2, y, p1, p2) {
  out = -1/(2*sigma2)*sum((y[2:T]-(mu+psi*y[1:(T-1)]))^2) + (p1-1)*log(psi) + (p2-1)*log(1-psi)
  return(exp(out))
}
#Starting iteration
mubar = sum(y[2:T]-chain[1,2]*(y[1:(T-1)]))/T
chain[1,1] = rnorm(1, mubar, sqrt(sig0/T))
old = psi0
candidate = rnorm(1, old, 0.15)
if(candidate < 0 | candidate > 1) {
  chain[1, 2] = old
} else {
  alpha = psi.kmcmc(candidate, chain[1,1], sig0, y, p1, p2)/psi.kmcmc(old, chain[1,1], sig0, y, p1, p2)
  alpha = min(alpha, 1)
  test = runif(1)
  if(alpha > test) {
    chain[1,2] = candidate
  } else {
    chain[1,2] = old
  }
}
chain[1,3] = 1/rgamma(1, shape=(T-1)/2, rate=sum((y[2:T]-(chain[1,1]+chain[1,2]*y[1:(T-1)]))^2)/2)

#Loop
for(i in 2:(B+K)) {
  mubar = sum(y[2:T]-chain[(i-1),2]*(y[1:(T-1)]))/T
  chain[i,1] = rnorm(1, mubar, sqrt(chain[(i-1),3]/T))
  old = chain[(i-1),2]
  candidate = rnorm(1, old, 0.15)
  if(candidate < 0 | candidate > 1) {
    chain[i, 2] = old
  } else {
    alpha = psi.kmcmc(candidate, chain[i,1], chain[(i-1),3], y, p1, p2)/psi.kmcmc(old, chain[i,1], chain[(i-1),3], y, p1, p2)
    alpha = min(1, alpha)
    test = runif(1)
    if(alpha > test) {
      chain[i,2] = candidate
    } else {
      chain[i,2] = old
    }
  }
  chain[i,3] = 1/rgamma(1, shape=(T-1)/2, rate=sum((y[2:T]-(chain[i,1]+chain[i,2]*y[1:(T-1)]))^2)/2)
}

ypred = data.frame(support = seq(-1, 8, 0.005))
ypred$density = rao(chain[1001:51000,], ypred$support, 50000, y[T])

p1 <- qplot(chain[,1], geom="density", xlim=c(1, 4)) + labs(x="p(Mu) Support", y="Density")
p2 <- qplot(chain[,2], geom="density", xlim=c(-0.25, 0.25))+ labs(x="p(Psi) Support", y="Density")
p3 <- qplot(chain[,3], geom="density", xlim=c(0, 2))+ labs(x="p(Sigma Squared) Support", y="Density")
q1 <- qplot(seq(1, 4, 0.005), dnorm(seq(1, 4, 0.005), barmu[1000], sqrt(lambda[1000])), geom="line")+ labs(x="q(Mu) Support", y="Density")
q2 <- qplot(seq(0, 1, 0.005), dnorm(seq(0, 1, 0.005), hatpsi[1000], sqrt(gamma[1000])), geom="line")+ labs(x="q(Psi) Support", y="Density")
q3 <- qplot(seq(0.005, 2, 0.005), dgamma(1/seq(0.005, 2, 0.005), shape=a, rate=b[1000]), geom="line")+ labs(x="q(Sigma Squared) Support", y="Density")
grid.arrange(p1, q1, p2, q2, p3, q3, ncol=2)
densities <- cbind(q = dens, p = ypred$density)
p4 <- ggplot() + geom_line(data=data.frame(sup=support, dens=dens), aes(x=sup, y=dens), colour="blue") + theme_bw() + 
  geom_line(data=data.frame(s = support, d = ypred$density), aes(x=s, y=d), colour="red") + labs(x="Support", y="Density")



  
                 
