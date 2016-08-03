#ARMA(1,1)
library(tidyr)
library(ggplot2)
library(dplyr)
library(coda)
library(MASS)

T <- 50
J <- 10
sigma <- 1.2
psi <- 0.3
rho <- 0.5
eps <- rnorm(T+J, 0, sigma)

y <- eps[1]
for(t in 2:60){
  y[t] <- rho*y[(t-1)] + psi*eps[(t-1)] + eps[t]
}

kalman <- function(psi, rho, sigma, n){
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
  v[1] <- y[2] -  rho*y[1] - t(H)%*%epsl[,1]
  F[1] <- t(H)%*%Pl[,,1]%*%H
  M[,1] <- Pl[,,1]%*%H
  eps[,1] <- epsl[,1] + M[,1]%*%solve(F[1])%*%v[1]
  P[,,1] <- Pl[,,1] - M[,1]%*%solve(F[1])%*%t(M[,1])
  
  #update loop
  for(t in 2:(n-1)){
    epsl[,t] <-  R%*%eps[,(t-1)]
    Pl[,,t] <- R%*%P[,,(t-1)]%*%t(R)+Q
    v[t] <- y[(t+1)] -  rho*y[(t)] - t(H)%*%epsl[,t]
    F[t] <- t(H)%*%Pl[,,t]%*%H
    M[,t] <- Pl[,,t]%*%H
    eps[,t] <- epsl[,t] + M[,t]%*%solve(F[t])%*%v[t]
    P[,,t] <- Pl[,,t] - M[,t]%*%solve(F[t])%*%t(M[,t])
  }
  epsdraw <- c(eps[2,1], eps[1,])
  epsdraw
}

psiinit <- 0.5
rhoinit <- 0.5
sigmainit <- 1
v <- T+J-2
betadraws <- matrix(0, nrow=51000, ncol=2)
sigdraws <- rep(0, 51000)

for(i in 1:51000){
  if(i ==1){
    eps <- kalman(psiinit, rhoinit, sigmainit, 60)
    X <- as.matrix(cbind(y[1:59], eps[1:59]))
    XX <- t(X)%*%X
    b <- solve(XX)%*%t(X)%*%y[2:60]
    sigsqhat <- sum((y[2:60]-X%*%b)^2)/v
    sigdraws[i] <- 1/sqrt(rgamma(1, v/2, sigsqhat*v/2))
    betadraws[i,] <- mvrnorm(1, b, solve(XX)*sigdraws[i]^2)
  } else {
    eps <- kalman(betadraws[(i-1), 2], betadraws[(i-1), 1], sigdraws[(i-1)], 60)
    X <- as.matrix(cbind(y[1:59], eps[1:59]))
    XX <- t(X)%*%X
    b <- solve(XX)%*%t(X)%*%y[2:60]
    sigsqhat <- sum((y[2:60]-X%*%b)^2)/v
    sigdraws[i] <- 1/sqrt(rgamma(1, v/2, sigsqhat*v/2))
    betadraws[i,] <- mvrnorm(1, b, solve(XX)*sigdraws[i]^2)
  }
}

Exact <- data.frame(psi = mcmc(betadraws[1001:51000,2]), rho= mcmc(betadraws[1001:51000,1]),sig = mcmc(sigdraws[1001:51000]))
colnames(Exact) <- c("psi", "rho", "sigma")
effectiveSize(Exact$psi)
effectiveSize(Exact$rho)
effectiveSize(Exact$sigma)

arma11 <- Exact[sample(50000, 2000, replace=TRUE),]

ar1l <- gather(AR1, theta, density)
ar2l <- gather(AR2, theta, density)
ma1l <- gather(MA1, theta, density)
ma2l <- gather(MA2, theta, density)
arma11l <- gather(arma11, theta, density)
ar1l <- mutate(ar1l, theta = ifelse(theta=="psi", "rho", theta))
ar2l <- mutate(ar2l, theta = ifelse(theta=="psi1", "rho", theta))
ar2l <- mutate(ar2l, theta = ifelse(theta=="psi2", "rho2", theta))
ma1l <- mutate(ma1l, theta = ifelse(theta=="psi1", "psi", theta))
ma2l <- mutate(ma2l, theta = ifelse(theta=="psi1", "psi", theta))
ar1l$model <- "AR1"
ar2l$model <- "AR2"
ma1l$model <- "MA1"
ma2l$model <- "MA2"
arma11l$model <- "ARMA"
all <- rbind(ar1l, ar2l, ma1l, ma2l, arma11l)
