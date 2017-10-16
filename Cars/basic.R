library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
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
sourceCpp('basic.cpp')
sigSqV <- 0.001
sigSqD <- 0.00001
phi <- 0
zeta <- 0.025
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
  d <- (1 - st) * d + st * (pi/2 + zeta * (x - 5)) + rnorm(1, 0, sqrt(sigSqD))
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
lambda <- matrix(c(-6, -11, 0.025, -5, 1, 2, 12, diag(0.001, 7)), 56)
hyper <- c(3, 0.03, 3, 0.00003, 0, 0.01, -5, 10, 1, 10, 0, 10, 10, 10)

carsHF <- function(data, lambda, hyper, S, maxIter, alpha, beta1, beta2, threshold = 0.01, Midpoint=5){
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
      logpj <- HFDeriv(data, lambda, epsilon[s,], hyper, Midpoint)
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
  return(list(lambda = lambda,
              ELBO = LB[1:min(iter, maxIter)],
              iter = iter))
}

transform <- function(x){
  x[1:2] = exp(x[1:2])
  x
}

fit <- carsHF(data, lambda, hyper, 100, 5000, 0.01, 0.9, 0.999, 0.01)
muF <- fit$lambda[1:7]
U <- matrix(fit$lambda[8:56], 7)
Sig <- t(U) %*% U
samples <- mvtnorm::rmvnorm(2000, muF, Sig)
samples <- t(apply(samples, 1, transform))
colnames(samples) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')

support <- matrix(0, ncol=7, nrow=5000)
dens <- matrix(0, ncol=7, nrow=5000)
for(i in 1:7){
  support[,i] <- seq(min(samples[,i]), max(samples[,i]), length.out=5000)
  if(i < 3){
    dens[,i] <- dlnorm(support[,i], muF[i], sqrt(Sig[i, i]))
  } else {
    dens[,i] <- dnorm(support[,i], muF[i], sqrt(Sig[i, i]))
  }
}
support <- data.frame(support)
dens <- data.frame(dens)
colnames(support) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
colnames(dens) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
supL <- gather(support, var, support)
densL <- gather(dens, var, density)
VB <- cbind(supL, density= densL$density)

ggplot(VB) + geom_line(aes(support, density)) +
  geom_vline(data = data.frame(true = c(sigSqV, sigSqD, zeta, a0, b0, a1, b1), var = colnames(samples)),
             aes(xintercept=true), colour = 'red') + 
  facet_wrap(~var, scales='free') + theme_bw() + 
  labs(title = 'VB - Fit is similar regardless of starting values') -> plotVB



reps <- 1000000

MCMC <- HFMCMC(data, reps, c(0.00005, 0.0000005, 0.001, 0.25, 0.25, 0.25, 0.25), hyper)

coda::effectiveSize(MCMC[(reps/2 +1):reps,])
MCMCdf <- as.data.frame(MCMC)
colnames(MCMCdf) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
MCMCdf$iter <- 1:reps
MCMCl <- gather(MCMCdf, var, draw, -iter)
MCMCl %>% 
  filter(iter > reps/2) %>%
  ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, ncol=1, scales='free')
MCMCl %>% 
  filter(iter > reps/2 & iter %% 10 == 0) %>%
  ggplot() + geom_density(aes(draw)) +
  geom_line(data=VB, aes(support, density), colour = 'blue') +
  geom_vline(data = data.frame(true = c(sigSqV, sigSqD, zeta, a0, b0, a1, b1), var = colnames(MCMCdf)[1:7]), aes(xintercept = true), colour='red') +
  facet_wrap(~var, scales='free') + 
  theme_bw() + 
  labs(title = 'PMMH - 100,000 draws, discard 50,000 - Acceptance rate ~ 15%') -> plotMH

gridExtra::grid.arrange(plotVB, plotMH, ncol=1)


N <- nrow(data)
xi <- numeric(N)
rho01 <- numeric(N)

rho11 <- numeric(N)
pv <- numeric(N)
pd0 <- numeric(N)
pd1 <- numeric(N)
likelihood <- numeric(N)
loglik = 0
xi[1] <- 0.95
for(t in 2:N){
  rho01[t] = 1.0 / (1 + exp(-(a0 + b0 * (data[t, 1] - 5)^2)))
  rho11[t] = 1.0 / (1 + exp(-(a1 + b1 * (data[t, 1] - 5)^2)))
  pv[t] = - 0.5 * log(2 * pi * sigSqV)  -  (data[t, 2] - data[t-1, 2])^2 / (2*sigSqV)
  pd0[t] = -0.5 * log(2 * pi * sigSqD)  -  (data[t, 3] - data[t-1, 3])^2 / (2*sigSqD)
  pd1[t] = -0.5 * log(2 * pi * sigSqD)  -  (data[t, 3] - pi/2 - zeta * (data[t, 1] - 5))^2 / (2*sigSqD)
  likelihood[t] =  xi[t-1] * (1 - rho01[t]) * exp(pv[t] + pd0[t])  +
    xi[t-1] * rho01[t] * exp(pv[t] + pd1[t])  +                        
    (1 - xi[t-1]) * (1 - rho11[t]) * exp(pv[t] + pd0[t])  +           
    (1 - xi[t-1]) * rho11[t] * exp(pv[t] + pd1[t])                   
  xi[t] = ((1 - rho01[t]) * xi[t-1] * exp(pv[t] + pd0[t])  +  (1 - rho11[t]) * (1 - xi[t-1]) * exp(pv[t] + pd0[t])) / likelihood[t]
  loglik = loglik + log(likelihood[t]);
}

df <- mutate(df, state1 = pi/2 + zeta * (xLag - 5),
             dens1 = 1/sqrt(2*pi*sigSqD) * exp(-(d - state1)^2 / (2*sigSqD)))
df$dlag <- c(0, df$d[1:(nrow(df)-1)])
df <- mutate(df, dens0 = 1/sqrt(2*pi*sigSqD) * exp(-(d - dlag)^2 / (2*sigSqD)))
ggplot(df) + geom_line(aes(1:nrow(df), dens1), colour='red') + geom_line(aes(1:nrow(df), dens2))

data <- as.matrix(df)
lambda <- matrix(c(-6, -11.7, 0.025, -5, 1, 3, 10, diag(0.001, 7)), 56)
hyper <- c(3, 0.03, 3, 0.00003, 0, 0.01, -5, 10, 1, 10, 0, 10, 10, 10)
HFDeriv(data, lambda, rnorm(7), hyper)




car10 <- filter(cars, ID == 10) %>% 
ggplot(car10) + geom_path(aes(x, y))
T <- nrow(car10)
df10 <- data.frame(xLag=car10$x[1:(T-1)], v=sqrt((car10$x[2:T] - car10$x[1:(T-1)])^2 + (car10$y[2:T] - car10$y[1:(T-1)])^2), d=atan2(car10$y[2:T] - car10$y[1:(T-1)], car10$x[2:T] - car10$x[1:(T-1)])) 
data <- as.matrix(df10)
ggplot(df10) + geom_line(aes(2:T, d))

lambda <- matrix(c(-5, -5, 0, -2, 1, -2, 1, diag(0.01, 7)), 56)
midpoints <- cars %>% 
  group_by(lane) %>%
  summarise(mid = mean(x))
midpoints

reps <- 500000
VB10 <- carsHF(data, lambda, hyper, 100, 5000, 0.01, 0.9, 0.99, 0.01, 6.97)
MCMC10 <- HFMCMC(data, reps, c(0.0005, 0.000005, 0.001, 0.35, 0.35, 0.35, 0.35), hyper, 6.97)

mu10 <- VB10$lambda[1:7]
U10 <- matrix(VB10$lambda[8:56], 7)
Sig10 <- t(U10) %*% U10
samples <- mvtnorm::rmvnorm(2000, mu10, Sig10)
samples <- t(apply(samples, 1, transform))
colnames(samples) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')

support <- matrix(0, ncol=7, nrow=5000)
dens <- matrix(0, ncol=7, nrow=5000)
for(i in 1:7){
  support[,i] <- seq(min(samples[,i]), max(samples[,i]), length.out=5000)
  if(i < 3){
    dens[,i] <- dlnorm(support[,i], mu10[i], sqrt(Sig10[i, i]))
  } else {
    dens[,i] <- dnorm(support[,i], mu10[i], sqrt(Sig10[i, i]))
  }
}
support <- data.frame(support)
dens <- data.frame(dens)
colnames(support) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
colnames(dens) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
supL <- gather(support, var, support)
densL <- gather(dens, var, density)
VB <- cbind(supL, density= densL$density)


coda::effectiveSize(MCMC10[(reps/2 +1):reps,])
MCMCdf <- as.data.frame(MCMC10)
colnames(MCMCdf) <- c('sigmaSq_epsilon', 'sigmaSq_eta', 'zeta', 'a0', 'b0', 'a1', 'b1')
MCMCdf$iter <- 1:reps
MCMCl <- gather(MCMCdf, var, draw, -iter)
MCMCl %>% 
  filter(iter > reps/2) %>%
  ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, ncol=1, scales='free')
MCMCl %>% 
  filter(iter > reps/2 & iter %% 10 == 0) %>%
  ggplot() + geom_density(aes(draw)) +
  geom_line(data=VB, aes(support, density), colour = 'blue') +
  #geom_vline(data = data.frame(true = c(sigSqV, sigSqD, zeta, a0, b0, a1, b1), var = colnames(MCMCdf)[1:7]), aes(xintercept = true), colour='red') +
  facet_wrap(~var, scales='free') + 
  theme_bw()


### Check out car lane change steering angles
cars %>% 
  filter(enterExit == FALSE) %>% 
  group_by(ID) %>%
  mutate(n = seq_along(x),
         xlag = ifelse(n == 1, 0, lag(x)),
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x - xlag)^2 + (y - ylag)^2),
         delta = atan2(y-ylag, x-xlag)) %>%
  filter(n > 1) -> carsAug

ggplot(carsAug) + geom_histogram(aes(delta)) + xlim(1.25, 1.9)

cars %>%
  filter(changed == TRUE &
         enterExit == FALSE) %>%
  group_by(ID) %>%
  mutate(n = seq_along(x),
         changeTime = ifelse(n == 1, 0, lane != lag(lane))) %>%
  ungroup() -> changing

changing %>% 
  filter(changeTime == 1) %>%
  select(ID, n) -> changeTimes

changing$nextChange = 9999
changingAug <- data.frame()
IDlist <- unique(changing$ID)
for(i in 1:length(IDlist)){
  temp1 <- filter(changeTimes, ID == IDlist[i])
  temp2 <- filter(changing, ID == IDlist[i])
  if(nrow(temp1) == 1){
    temp2$nextChange[temp2$n <= temp1$n] = temp1$n
  } else {
    temp2$nextChange[temp2$n <= temp1$n[1]] = temp1$n[1]
    for(j in 2:nrow(temp1)){
      temp2$nextChange[temp2$n <= temp1$n[j] & temp2$n > temp1$n[j-1]] = temp1$n[j]
    }
  }
  changingAug <- rbind(changingAug, temp2)
}

changingAug %>% 
  select(ID, nextChange, n, x, y, lane) %>%
  mutate(xlag = ifelse(n == 1, 0, lag(x)),
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x - xlag)^2 + (y - ylag)^2),
         delta = atan2(y-ylag, x-xlag), 
         timeToChange = nextChange - n) %>%
  filter(n > 1) -> deltaChange

deltaChange %>%
  filter(timeToChange <= 50) %>%
  mutate(frame = factor(51 - timeToChange)) %>%
  ggplot() + geom_density(aes(x=delta, frame=frame)) + xlim(pi/2 - 0.4, pi/2 + 0.4) -> p
animation::ani.options(interval = 0.1)
gganimate::gganimate(p)


filter(carsAug, ID <= 10) %>%
  select(ID, x, delta, frame) %>%
  gather(var, value, -ID, -frame) %>%
  ggplot() + geom_path(aes(value, frame)) + facet_grid(ID ~ var, scales = 'free') + theme_bw()
  
filter(carsAug, ID <= 10) %>%
  select(ID, x, v, xlag, delta, frame) %>%
  mutate(v = round(v, 1)) %>%
  ggplot() + geom_path(aes(x - xlag, delta, colour=factor(v))) + facet_wrap(~ID, scales='free')


#### General Exploration
cars %>% 
  filter(changed == FALSE) %>%
  .$ID %>%
  unique() -> noChange

compareData <- function(sigSqV = 0.01, sigSqD = 0.00001, zeta = 0.025, 
                        a0 = -5.72, b0 = 0.88, a1 = 2, b = 12, dfD = 5, dfV = 5, error = 'Normal'){
  pos <- data.frame()
  v <- 5
  vLag1 <- 4.9
  d <- pi/2
  x0 <- 0
  x <- x0
  y <- 0
  st <- 0
  while(y < 500){
    if(t > 1000){
      break
    }
    if(st == 0){
      if(runif(1) < 1 / (1 + exp(-(a0 + b0 * (x - x0)^2)))){
        st <- 1
      }
    } else {
      if(runif(1) > 1 / (1 + exp(-(a1 + b1 * (x - x0)^2)))){
        st <- 0
      }
    }
    d <- (1 - st) * d + st * (pi/2 + zeta * (x - x0))
    vLag2 <- vLag1
    vLag1 <- v
    v <- vLag1 + phi * (vLag1 - vLag2)
    if(error == 'Normal'){
      d <- d + rnorm(1, 0, sqrt(sigSqD))
      v <- v + rnorm(1, 0, sqrt(sigSqV))
    } else {
      d <- d + sqrt(sigSqD) * rt(1, dfD)
      v <- v + sqrt(sigSqV) * rt(1, dfV)
    }
    x <- x + v * cos(d)
    y <- y + v * sin(d)
    pos <- rbind(pos, data.frame(x=x, y=y, frame=t, data = 'sim'))
  }
  cars %>% 
    filter(ID == sample(noChange, 1) &
           y < 500) %>%
    mutate(x = x - recode(lane, 6.97, 18.65, 29.932, 41.03, 52.93, 63.66, 65.88, 67.5)) %>%
    select(x, y, frame) %>%
    mutate(frame = frame - min(.$frame) + 1,
           data = 'real') %>%
    rbind(pos) %>%
    group_by(data) %>%
    mutate(x = x - mean(x),
           xlag = ifelse(frame == 1, 0, lag(x)),
           ylag = ifelse(frame == 1, 0, lag(y)),
           v = sqrt((x - xlag)^2 + (y - ylag)^2),
           delta = atan2(y-ylag, x-xlag)) %>%
    ungroup() %>%
    filter(frame > 1) %>%
    select(x, y, v, delta, frame, data) %>%
    ggplot() + geom_point(aes(delta, lag(x)^2)) + facet_wrap(~data) %>%
    #gather(var, value, -frame, -data, -y) %>%
    #ggplot() + geom_path(aes(value, y, colour = data)) + facet_wrap(~var, ncol=1, scales='free') %>%
    print()
}
compareData(dfV = 5, dfD = 5, sigSqD = 1e-05, error = 't')





