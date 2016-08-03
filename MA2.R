##MA(2) Model
library(tidyr)
library(ggplot2)
library(dplyr)
library(MASS)
library(coda)


psi1 <- 0.5
psi2 <- 0.3
sigma <- 0.8
T <- 50
J <- 10
eps <- rnorm(T+J, 0, sigma)


y <- eps[1]
y[2] <- psi1*eps[1]+eps[2]
for(t in 3:60){
  y[t] <- eps[t] + psi1*eps[(t-1)] + psi2*(eps[(t-2)])
}


#kalman filter, from t = 2 to n
kalman <- function(theta, n){
  psi1 <- theta[1]
  psi2 <- theta[2]
  sigma <- theta[3]
  eps <- matrix(0, ncol=n-1, nrow=3) #eps(t|t)
  epsl <- matrix(0, ncol=n-1, nrow=3) #eps(t|t-1)
  P <- array(0, dim=c(3,3,n-1)) #P(t|t)
  Pl <- array(0, dim=c(3,3,n-1)) #P(t|t-1)
  R <- matrix(0,3,3) #formerly F
  R[2,1] <- 1
  R[3,2] <- 1
  H <- c(1, psi1, psi2)
  Q <- matrix(c(sigma,rep(0,8)), 3)
  v <- rep(0,n-1) #y - mean(y)
  F <- rep(0,n-1) #var(y)
  M <- matrix(0, ncol=n-1, nrow=3)
  
  #time 1 - #eps0 (error in time1) starts at y1, P0 is a zero matrix (everything is known about y1)
  epsl[,1] <- R%*%c(y[1],0,0)
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
init.psi1 <- 0.3
init.psi2 <- 0.3
init.sig <- 1
init.eps <- kalman(c(init.psi1, init.psi2, init.sig), 60)
init.X <- cbind(init.eps[2:59], init.eps[1:58])
init.XX <- t(init.X)%*%init.X
init.b <- solve(init.XX)%*%t(init.X)%*%y[3:60]

sigdraws <- rep(0, 21000)
psidraws <- matrix(0, ncol=2, nrow=21000)

for(i in 1:21000){
  if(i == 1){
    sigsqhat <- sum((y[3:60]-init.X%*%init.b)^2)/56
    sigdraws[1] <- 1/sqrt(rgamma(1, 56/2, sigsqhat*56/2))
    psidraws[1,] <- mvrnorm(1, init.b, solve(init.XX)*sigdraws[1]^2)
  } else {
    eps.star <- kalman(c(psidraws[(i-1),], sigdraws[(i-1)]), 60)
    X <- cbind(eps.star[2:59], eps.star[1:58])
    XX <- t(X)%*%X
    b <- solve(XX)%*%t(X)%*%y[3:60]
    sigsqhat <- sum((y[3:60]-X%*%b)^2)/56
    sigdraws[i] <- 1/sqrt(rgamma(1, 56/2, sigsqhat*56/2))
    psidraws[i,] <- mvrnorm(1, b, solve(XX)*sigdraws[i]^2)
  }
}

Exact <- data.frame(psi = mcmc(psidraws[1001:21000,]), sig = mcmc(sigdraws[1001:21000]))
colnames(Exact) <- c("psi1", "psi2", "sigma")
effectiveSize(Exact$psi1)
effectiveSize(Exact$psi2)
effectiveSize(Exact$sigma)

#Prior
init.psi1 <- 0.4
init.psi2 <- 0.4
init.sig <- 1
init.eps <- kalman(c(init.psi1, init.psi2, init.sig), 50)
init.X <- cbind(init.eps[2:49], init.eps[1:48])
init.XX <- t(init.X)%*%init.X
init.b <- solve(init.XX)%*%t(init.X)%*%y[3:50]

sigdraws <- rep(0, 11000)
psidraws <- matrix(0, ncol=2, nrow=11000)

for(i in 1:11000){
  if(i == 1){
    sigsqhat <- sum((y[3:50]-init.X%*%init.b)^2)/46
    sigdraws[1] <- 1/sqrt(rgamma(1, 46/2, sigsqhat*46/2))
    psidraws[1,] <- mvrnorm(1, init.b, solve(init.XX)*sigdraws[1]^2)
  } else {
    eps.star <- kalman(c(psidraws[(i-1),], sigdraws[(i-1)]), 50)
    X <- cbind(eps.star[2:49], eps.star[1:48])
    XX <- t(X)%*%X
    b <- solve(XX)%*%t(X)%*%y[3:50]
    sigsqhat <- sum((y[3:50]-X%*%b)^2)/46
    sigdraws[i] <- 1/sqrt(rgamma(1, 46/2, sigsqhat*46/2))
    psidraws[i,] <- mvrnorm(1, b, solve(XX)*sigdraws[i]^2)
  }
}

Prior <- data.frame(psi = mcmc(psidraws[1001:11000,]), sig = mcmc(sigdraws[1001:11000]))
colnames(Prior) <- c("psi1", "psi2", "sigma")
effectiveSize(Prior$psi1)
effectiveSize(Prior$psi2)
effectiveSize(Prior$sigma)

##ABC Scheme
sum.stats2 <- function(x) {
  out <- c(cor(x[1:9], x[2:10]), cor(x[1:8], x[3:10]), cor(x[1:7], x[4:10]), var(x))
  out
}

ma2 <- function(theta){
  out <- rep(0, 10)
  out[1] <- theta[1]*y[50] + theta[2]*y[49] + rnorm(1, sd=theta[3])
  out[2] <- theta[1]*out[1] + theta[2]*y[50] + rnorm(1, sd=theta[3])
  for(i in 3:10){
    out[i] <- theta[1]*out[(i-1)] + theta[2]*out[(i-2)] + rnorm(1, sd=theta[3])
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

yss <- sum.stats2(y[51:60])
samp <- sample(1:10000, 100000, replace=TRUE)
theta <- Prior[samp,]
z <- apply(theta, 1, ma2)
zss <- apply(z, 2, sum.stats2)
dist <- apply(zss, 2, distance)
accept <- theta[index.bottom.N(dist, 1000),]
abc <- data.frame(psi1 = accept[,1], psi2 = accept[,2], sigma = accept[,3], method = "ABC")
prior <- data.frame(Prior[sample(1:10000, 1000, replace=TRUE),], method = rep("Prior", 1000))
exact <- data.frame(Exact[sample(1:10000, 1000, replace=TRUE),], method = rep("Exact", 1000))
MA2 <- rbind(prior, exact, abc)
MA2l <- gather(MA2, theta, density, -method)
ggplot(data=MA2l, aes(x=density, colour=method)) + geom_density() + facet_wrap(~theta, scales="free")
