library(tidyverse)
Sigma = 0.01 * matrix(c(1, 0, 0, 1), 2)
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
d0 = pi/2

T = 400
path = data.frame()
for(i in 1:10){
  position = data.frame()
  x0 = 5
  notice = 0 
  delta = pi/2 + lambda * (d0 - pi/2)  +  sigZ * rt(1, 5)
  acc = phi * a0  +  sigE * rt(1, 5)
  vel = v0  +  acc
  error = mvtnorm::rmvnorm(1, c(0, 0), Sigma)
  x = x0  +  vel * cos(delta)
  y = y0  +  vel * sin(delta)
  measurex = x + error[1]
  measurey = y + error[2]
  position = rbind(position,
                   data.frame(x=x, y=y, measurex=measurex, measurey=measurey, vel=vel, delta=delta, acc=acc, notice=notice, car = i))
  for(t in 2:T){
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
p1 <- ggplot(path) + geom_path(aes(x, y, group=car, colour=factor(car))) + 
  labs(x="Simulated", y=NULL)
print(p1)
#cars %>%
#  filter(changed == FALSE) %>%
#  filter(ID %in% head(sort(unique(.$ID)), 10)) %>%
#  ggplot() + geom_path(aes(x, y, colour = factor(ID))) + theme(legend.position = 'none') +
#  labs(x="Actual", y=NULL) -> p2
#gridExtra::grid.arrange(p1, p2)




hyper = c(10, .1, 200, 0.1, 10, .1, 10, .1, 20, 1.5, 20, 1.5, 0, 0.1, -3, 5, 10, 5, -3, 5, 0, 5)
lambda <- as.matrix(c(-4, -8, -4, -4, 2, 2, 0, -3, 5, -3, 1, rep(0.1, 11)), 22)
#VB_Cars2(as.matrix(position[1:20, 3:4]), lM, hyper, maxIter=1, S=1, N=1)
#posi = as.matrix(position[, 3:4]) 
#MCMC = PMMH(pos= posi, reps=10, N=1, hyperParams=hyper, stepSize=rep(0.005, 8))

data = as.matrix(select(car10Path, x, y))#matrix(position[,3:4])
N = 100
S = 10
maxIter = 1000
alpha = 0.1
beta1 = 0.9
beta2 = 0.99
threshold = 0.01

carsVB <- function(data, lambda, hyper, N, S, maxIter, alpha, beta1, beta2, threshold){
  T <- nrow(data)
  sobol <- sobol_points(100+S, 11)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- numeric(22)
  V <- numeric(22)
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
    unif <- shuffle(sobol)
    epsilon <- qnorm(unif[101:(100+S), ])
    grad <- matrix(0, 22, S)
    eval <- numeric(S)
    q <- numeric(S)
    #car <- filter(data, car==1)#sample(1:max(data$car), 1))
    for(s in 1:S){
      logpj <- PFDeriv(data[1:100,], lambda, epsilon[s,], hyper, 1000)
      eval[s] <- logpj$val
      grad[,s] <- logpj$grad
      q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
    }
    gradient <- rowMeans(grad, na.rm = TRUE)
    gradientSq <- rowMeans(grad^2, na.rm=TRUE)
    valP[valP == -Inf] = NA
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
      print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  return(lambda)
}

fit <- carsVB(data, lambda, hyper, N=100, S=10, maxIter=100, alpha=0.25, beta1=0.9, beta2=0.99, threshold=0.01)


muFit <- fit[1:11]# c(-5.1963258, -7.5416080, -4.8451232, -4.1386897, 3.1992270, 2.3012937, 1.1962711, -4.1964376, -0.1966839, -4.1963091, -0.1966602)
sigmaFit <- fit[12:22] #c(0.5590673, -0.3585571, 0.1492671, 0.1903151, -1.0920861, -0.2012837, 0.5588602, 0.5592747, 1.2967022, -0.3601970,  0.5594020)
Sigma <- diag(sigmaFit^2)
transform <- function(x){
  result <- numeric(11)
  result[1:4] = exp(x[1:4])
  result[5:6] = 1/(1 + exp(-x[5:6]))
  result[7:11] = x[7:11]
  result
}
samples <- mvtnorm::rmvnorm(1000, muFit, Sigma)
samples <- t(apply(samples, 1, transform))
colnames(samples) <- c('SigA', 'SigD', 'SigX', 'SigY', 'lambda', 'phi', 'gamma', 'a0', 'b0', 'a1', 'b1')
samples <- gather(as.data.frame(samples), var, draw)
ggplot(samples) + geom_density(aes(draw)) + facet_wrap(~var, scales='free')
