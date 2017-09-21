library(tidyverse)
set.seed(465786)
sigSqV <- 0.001
sigSqD <- 0.001
phi <- 0.9
gamma <- 0.9
T <- 500
pos <- data.frame()

v = 5
vLag1 = 4.9
d = pi/2
x = 5
y = 0
st = 0


for(t in 1:T){
  d = pi/2 + gamma * (d - pi/2) + rnorm(1, 0, sqrt(sigSqD))
  vLag2 = vLag1
  vLag1 = v
  v = vLag1 + phi * (vLag1 - vLag2) + rnorm(1, 0, sqrt(sigSqV))
  x = x + v * cos(d)
  y = y + v * sin(d)
  pos <- rbind(pos, data.frame(x=x, y=y, v=v, d=d, t=t))
}

ggplot(pos) + geom_path(aes(x, y))

mu = rep(0, 4)
u = diag(0.1, 4)
lambda = as.matrix(c(mu, u))
hyper = c(3, 0.003, 3, 0.003, 20, 1.5, 20, 1.5)
data = data.frame(dx = pos$x[2:T] - pos$x[1:(T-1)], dy =  pos$y[2:T] - pos$y[1:(T-1)])
data = mutate(data, vt = sqrt(dx^2 + dy^2), dt = atan2(dy, dx)) %>% select(vt, dt) %>% as.matrix() 
alpha = 0.01
beta1 = 0.9
beta2 = 0.99
threshold = 0.01
S = 100
maxIter = 100

carsBasic <- function(data, lambda, hyper, S, maxIter, alpha, beta1, beta2, threshold = 0.01, thresholdIS = 0){
  sobol <- sobol_points(100+S, 4)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- numeric(20)
  V <- numeric(20)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    if(any(is.na(lambda))){
      break
    }
    grad <- matrix(0, 20, S)
    eval <- numeric(S)
    q <- numeric(S)
    unif <- shuffle(sobol)
    epsilon <- qnorm(unif[101:(100+S), ])
    for(s in 1:S){
      logpj <- PJDeriv(data, lambda, epsilon[s,], hyper)
      eval[s] <- logpj$val
      grad[,s] <- logpj$grad
      q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
    }
    gradient <- rowMeans(grad, na.rm = TRUE)
    gradientSq <- rowMeans(grad^2, na.rm=TRUE)
    LB[iter] <- mean(eval - q, na.rm=TRUE) 
  
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 25 == 0){
      print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  return(lambda)
}

fit <- carsBasic(data, lambda, hyper, 25, 5000, 0.15, 0.25, 0.75, 0.005)

muF <- fit[1:4]
U <- matrix(fit[5:20], 4)
Sig <- t(U) %*% U
samples <- mvtnorm::rmvnorm(2000, muF, Sig)
transform <- function(x){
  result <- numeric(4)
  result[1:2] = exp(x[1:2])
  result[3:4] = 1/(1 + exp(-x[3:4]))
  result
}
samples <- t(apply(samples, 1, transform))
colnames(samples) <- c('sigma^2[epsilon]', 'sigma^2[eta]', 'phi', 'gamma')
samplesL <- gather(as.data.frame(samples), var, draw)
ggplot(samplesL) + geom_density(aes(draw)) + 
  geom_vline(data = data.frame(true = c(sigSqV, sigSqD, phi, gamma), var = colnames(samples)),
             aes(xintercept=true), colour = 'red') + 
  facet_wrap(~var, scales='free') + theme_bw()
samplesL %>% group_by(var) %>% summarise(mean(draw))



## Model Two

sigSqV <- 0.001
sigSqD <- 0.00001
phi <- 0
zeta <- 0.05
a0 <- -5.72
b0 <- 0.88
a1 <- 2
b1 <- 12
T <- 500
pos <- data.frame()

v <- 5
vLag1 <- 4.9
d <- pi/2
x <- 5
y <- 0
st <- 0

for(t in 1:T){
  if(st == 0){
    if(runif(1) < 1 / (1 + exp(-(a0 + b0 * (x - 5)^2)))){
      st <- 1
    }
  } else {
    if(runif(1) > 1 / (1 + exp(-(a1 + b1 * (x - 5)^2)))){
      st <- 0
    }
  }
  d <- (1 - st) * d + st * (pi/2 + zeta * (x - 5)) + sqrt(sigSqD) * rt(1, 15)#rnorm(1, 0, sqrt(sigSqD))
  vLag2 <- vLag1
  vLag1 <- v
  v <- vLag1 + phi * (vLag1 - vLag2) + rnorm(1, 0, sqrt(sigSqV))
  x <- x + v * cos(d)
  y <- y + v * sin(d)
  pos <- rbind(pos, data.frame(x=x, y=y, v=v, d=d, st=st, t=t))
}

ggplot(pos) + geom_path(aes(x, y, group = 1, colour = factor(st))) + xlim(0, 10)

df <- data.frame(xLag=pos$x[1:(T-1)], v=sqrt((pos$x[2:T] - pos$x[1:(T-1)])^2 + (pos$y[2:T] - pos$y[1:(T-1)])^2), d=atan2(pos$y[2:T] - pos$y[1:(T-1)], pos$x[2:T] - pos$x[1:(T-1)])) 
data <- as.matrix(df)
lambda <- matrix(c(rep(0, 7), diag(0.01, 7)), 56)
hyper <- c(3, 0.03, 3, 0.00003, 0, 0.01, -5, 10, 1, 10, 0, 10, 10, 10)

carsHF <- function(data, lambda, hyper, S, maxIter, alpha, beta1, beta2, threshold = 0.01){
  sobol <- sobol_points(100+S, 7)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- numeric(56)
  V <- numeric(56)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    if(any(is.na(lambda))){
      break
    }
    grad <- matrix(0, 56, S)
    eval <- numeric(S)
    q <- numeric(S)
    unif <- shuffle(sobol)
    epsilon <- qnorm(unif[101:(100+S), ])
    for(s in 1:S){
      logpj <- HFDeriv(data, lambda, epsilon[s,], hyper)
      eval[s] <- logpj$val
      grad[,s] <- logpj$grad
      q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
    }
    gradient <- rowMeans(grad, na.rm = TRUE)
    gradientSq <- rowMeans(grad^2, na.rm=TRUE)
    LB[iter] <- mean(eval - q, na.rm=TRUE) 
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 25 == 0){
      print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  return(lambda)
}

fit <- carsHF(data, lambda, hyper, 5, 5000, 1, 0.01, 0.999, 0.01)
muF <- fit[1:7]
U <- matrix(fit[8:56], 7)
Sig <- t(U) %*% U
samples <- mvtnorm::rmvnorm(2000, muF, Sig)
transform <- function(x){
  result <- numeric(7)
  result[1:2] = exp(x[1:2])
  result[3:7] = x[3:7]
  result
}
samples <- t(apply(samples, 1, transform))
colnames(samples) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
samplesL <- gather(as.data.frame(samples), var, draw)
ggplot(samplesL) + geom_density(aes(draw)) + 
  geom_vline(data = data.frame(true = c(sigSqV, sigSqD, zeta, a0, b0, a1, b1), var = colnames(samples)),
             aes(xintercept=true), colour = 'red') + 
  facet_wrap(~var, scales='free') + theme_bw()
samplesL %>% group_by(var) %>% summarise(mean(draw))

