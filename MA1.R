##MA(1) Model
library(tidyr)
library(ggplot2)
library(dplyr)
library(MASS)
library(coda)


psi <- 0.5
sigma <- 0.8
T <- 50
J <- 10
eps <- rnorm(T+J, 0, sigma)


y <- eps[1]
for(t in 2:60){
  y[t] <- eps[t] + psi*eps[(t-1)]
}


#kalman filter, from t = 2 to n
kalman <- function(psi, sigma, n){
  eps <- matrix(0, ncol=n-1, nrow=2) #eps(t|t)
  epsl <- matrix(0, ncol=n-1, nrow=2) #eps(t|t-1)
  P <- array(0, dim=c(2,2,n-1)) #P(t|t)
  Pl <- array(0, dim=c(2,2,n-1)) #P(t|t-1)
  R <- matrix(c(0,1,0,0),2) #formerly F
  H <- c(1, psi)
  Q <- matrix(c(sigma,0,0,0), 2)
  v <- rep(0,n-1) #y - mean(y)
  F <- rep(0,n-1) #var(y)
  M <- matrix(0, ncol=n-1, nrow=2)
  
  #time 1 - #eps0 (error in time1) starts at y1, P0 is a zero matrix (everything is known about y1)
  epsl[,1] <- R%*%c(y[1],0)
  Pl[,,1] <- Q
  v[1] <- y[2] - t(H)%*%epsl[,1]
  F[1] <- t(H)%*%Pl[,,1]%*%H
  M[,1] <- Pl[,,1]%*%H
  eps[,1] <- epsl[,1] + M[,1]%*%solve(F[1])%*%v[1]
  P[,,1] <- Pl[,,1] - M[,1]%*%solve(F[1])%*%t(M[,1])
  
  #update loop
  for(t in 2:(n-1)){
    epsl[,t] <-  R%*%eps[,(t-1)]
    Pl[,,t] <- R%*%P[,,(t-1)]%*%t(R)+Q
    v[t] <- y[(t+1)] - t(H)%*%epsl[,t]
    F[t] <- t(H)%*%Pl[,,t]%*%H
    M[,t] <- Pl[,,t]%*%H
    eps[,t] <- epsl[,t] + M[,t]%*%solve(F[t])%*%v[t]
    P[,,t] <- Pl[,,t] - M[,t]%*%solve(F[t])%*%t(M[,t])
  }
  epsdraw <- c(eps[2,1], eps[1,])
  epsdraw
}


##MCMC Scheme
init.psi <- 0.8
init.sig <- 1
init.eps <- kalman(init.psi, init.sig, 60)
init.X <- init.eps[1:59]
init.XX <- t(init.X)%*%init.X
init.b <- solve(init.XX)%*%t(init.X)%*%y[2:60]

sigdraws <- rep(0, 21000)
psidraws <- rep(0, 21000)

for(i in 1:21000){
  if(i == 1){
    sigsqhat <- sum((y[2:60]-init.X*init.b)^2)/58
    sigdraws[1] <- 1/sqrt(rgamma(1, 58/2, sigsqhat*58/2))
    psidraws[1] <- rnorm(1, init.b, sqrt(sigdraws[1]^2/init.XX))
  } else {
    eps.star <- kalman(psidraws[(i-1)], sigdraws[(i-1)], 60)
    X <- eps.star[1:59]
    XX <- t(X)%*%X
    b <- solve(XX)%*%t(X)%*%y[2:60]
    sigsqhat <- sum((y[2:60]-X*b)^2)/58
    sigdraws[i] <- 1/sqrt(rgamma(1, 58/2, sigsqhat*58/2))
    psidraws[i] <- rnorm(1, b, sqrt(sigdraws[i]^2/XX))
  }
}

Gibbs <- data.frame(psi = mcmc(psidraws[1001:21000]), sig = mcmc(sigdraws[1001:21000]))
colnames(Gibbs) <- c("psi", "sigma")
effectiveSize(Gibbs$psi)
effectiveSize(Gibbs$sigma)

#Prior
init.psi <- 0.8
init.sig <- 1
init.eps <- kalman(init.psi, init.sig, 50)
init.X <- init.eps[1:49]
init.XX <- t(init.X)%*%init.X
init.b <- solve(init.XX)%*%t(init.X)%*%y[2:50]

sigdraws <- rep(0, 11000)
psidraws <- rep(0, 11000)

for(i in 1:11000){
  if(i == 1){
    sigsqhat <- sum((y[2:50]-init.X*init.b)^2)/48
    sigdraws[1] <- 1/sqrt(rgamma(1, 48/2, sigsqhat*48/2))
    psidraws[1] <- rnorm(1, init.b, sqrt(sigdraws[1]^2/init.XX))
  } else {
    eps.star <- kalman(psidraws[(i-1)], sigdraws[(i-1)], 50)
    X <- eps.star[1:49]
    XX <- t(X)%*%X
    b <- solve(XX)%*%t(X)%*%y[2:50]
    sigsqhat <- sum((y[2:50]-X*b)^2)/48
    sigdraws[i] <- 1/sqrt(rgamma(1, 48/2, sigsqhat*48/2))
    psidraws[i] <- rnorm(1, b, sqrt(sigdraws[i]^2/XX))
  }
}

Prior <- data.frame(psi = mcmc(psidraws[1001:11000]), sig = mcmc(sigdraws[1001:11000]))
colnames(Prior) <- c("psi", "sigma")
effectiveSize(Prior$psi)
effectiveSize(Prior$sigma)

##ABC Scheme
sum.stats <- function(x) {
  out <- c(cor(x[1:9], x[2:10]), var(x))
  out
}

ma1 <- function(theta){
  out <- rep(0, 10)
  out[1] <- theta[1]*y[50] + rnorm(1, sd=theta[2])
  for(i in 2:10){
    out[i] <- theta[1]*out[(i-1)] + rnorm(1, sd=theta[2])
  }
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

yss <- sum.stats(y[51:60])
samp <- sample(1:10000, 100000, replace=TRUE)
theta <- Prior[samp,]
z <- apply(theta, 1, ma1)
zss <- apply(z, 2, sum.stats)
dist <- apply(zss, 2, distance)
accept <- theta[index.bottom.N(dist, 1000),]
abc <- data.frame(psi = accept[,1], sigma = accept[,2], method = "ABC")
prior <- data.frame(Prior[sample(1:10000, 1000, replace=TRUE),], method = rep("Prior", 1000))
exact <- data.frame(Gibbs[sample(1:10000, 1000, replace=TRUE),], method = rep("Exact", 1000))
MA1 <- rbind(prior, exact, abc)
MA1l <- gather(MA1, theta, density, -method)
ggplot(data=MA1l, aes(x=density, colour=method)) + geom_density() + facet_wrap(~theta, scales="free")
