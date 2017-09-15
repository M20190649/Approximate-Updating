library(tidyverse)
Sigma = 0.01 * matrix(c(1, 0.5, 0.5, 1), 2)
sigE = 0.01
sigZ = 0.0005
phi = 0.9
lambda = 0.95
gamma = 0.005
alpha0 = -4
beta0 = 10
alpha1 = -4
beta1 = 2

a0 = rnorm(1, 0, sigE / sqrt(1 - phi^2))
v0 = 5
x0 = 5
y0 = 0
t0 = pi/2

T = 500
path = data.frame()
for(i in 1:1){
  position = data.frame()
  x0 = 5*i
  notice = 0 
  delta = pi/2 + lambda * (t0 - pi/2)  +  sigZ * rt(1, 5)
  acc = phi * a0  +  sigE * rt(1, 5)
  vel = v0  +  acc
  error = mvtnorm::rmvnorm(1, c(0, 0), Sigma)
  x = x0  +  vel * cos(delta)
  y = y0  +  vel * sin(delta)
  measurex = x + error[1]
  measurey = y + error[2]
  position = rbind(position,
                   data.frame(x=x, y=y, measurex=measurex, measurey=measurey, vel=vel, delta=delta, acc=acc, notice=notice, car = i))
  while(y < 2000) {
    notice = 1 / (1 + exp(-(1-notice)*(alpha0 + beta0*abs(x-x0)) - notice*(alpha1 + beta1*abs(x-x0)) + rnorm(1, 0, 0.1)))
    delta = pi/2  + (1-notice) * lambda * (delta- pi/2)  + notice * gamma * (x-x0) +  sigZ * rt(1, 5)
    acc = phi * acc  +  sigE * rt(1, 5)
    vel = vel  +  acc
    error = mvtnorm::rmvnorm(1, c(0, 0), Sigma)
    x = x  +  vel * cos(delta)
    y = y  +  vel * sin(delta)
    measurex = x + error[1]
    measurey = y + error[2]
    position = rbind(position,
                     data.frame(x=x, y=y, measurex=measurex, measurey=measurey, vel=vel, delta=delta, acc=acc, notice=notice, car = i))
  }
  path = rbind(path, position)
}
p1 <- ggplot(path) + geom_path(aes(x, y, group=car, colour=notice)) + 
  labs(x="Simulated", y=NULL)
print(p1)
#cars %>%
#  filter(changed == FALSE) %>%
#  filter(ID %in% head(sort(unique(.$ID)), 10)) %>%
#  ggplot() + geom_path(aes(x, y, colour = factor(ID))) + theme(legend.position = 'none') +
#  labs(x="Actual", y=NULL) -> p2
#gridExtra::grid.arrange(p1, p2)




hyper = c(10, .1, 200, 0.1, 10, .1, 10, .1, 20, 1.5, 20, 1.5, 0, 0.1, -3, 5, 0, 5, -3, 5, 0, 5)
lambda <- as.matrix(c(-4, -8, -4, -4, 2, 2, 0, -3, 1, -3, 1, rep(0.1, 11)), 22)
#VB_Cars2(as.matrix(position[1:20, 3:4]), lM, hyper, maxIter=1, S=1, N=1)
#posi = as.matrix(position[, 3:4]) 
#MCMC = PMMH(pos= posi, reps=10, N=1, hyperParams=hyper, stepSize=rep(0.005, 8))

data <- as.matrix(position[1:100, 3:4])
N = 100
S = 10
maxIter = 100
alpha = 0.1
threshold = 0.01

carsVB <- function(data, lambda, hyper, N, S, maxIter, alpha, threshold){
  T <- nrow(data)
  sobol <- sobol_points(100+S, 11)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- numeric(22)
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    if(any(is.na(lambda))){
      break
    }
    unif <- shuffle(sobol)
    epsilon <- qnorm(unif[101:(100+S), ])
    gradP <- matrix(0, 22, S)
    valP <- numeric(S)
    gradJ <- matrix(0, 22, S)
    valJ <- numeric(S)
    q <- numeric(S)
    for(s in 1:S){
      p <- pDeriv(data[1:100,], lambda, epsilon[s,], hyper, N)
      valP[s] <- p$val
      gradP[,s] <- p$grad
      j <- jDeriv(lambda, epsilon[s,])
      valJ[s] <- j$val
      gradJ[,s] <- j$grad
      q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
    }
    gradient <- rowMeans(gradP + gradJ, na.rm = TRUE)
    valP[valP == -Inf] = NA
    LB[iter] <- mean(valP + valJ - q, na.rm=TRUE) 
    M <- M + gradient^2
    if(iter > 1){
      lambda = lambda + alpha * M^(-0.5) * gradient
    }
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
      print(meanLB)
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  return(lambda)
}

carsVB(data, lambda, hyper, N=100, S=10, maxIter=1500, alpha=0.25, threshold=0.01)
