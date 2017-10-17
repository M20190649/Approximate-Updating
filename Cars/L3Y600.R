library(tidyverse)
library(forecast)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
source('carsVBfuns.R')
sourceCpp('basic.cpp')
sourceCpp('heirarchical.cpp')
L3Y600 <- readr::read_csv('L3Y600.csv')
L3Y600 %>% 
  group_by(ID) %>%
  filter(min(v) > 0.1 & changed == FALSE) %>%
  .$ID %>%
  unique() -> noChange

plots{
L3Y600 %>%
  filter(changed == FALSE) %>%
  group_by(ID) %>%
  mutate(n = seq_along(time)) %>%
  filter(min(delta) > 0) %>%
  filter(ID %in% head(unique(.$ID), 24)) %>%
  ggplot() + geom_line(aes(n, delta)) +
    facet_wrap(~ID, scales='free_x', ncol = 8) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.ticks.y = element_blank(),
          strip.background = element_blank()) +
    labs(x = NULL) -> p1

L3Y600 %>%
  group_by(ID) %>%
  mutate(n = seq_along(time), lv = lag(v)) %>%
  filter(n > 1 & changed == TRUE & min(delta) > 0) %>%
  filter(ID %in% tail(unique(.$ID), 24)) %>%
  ggplot() + geom_line(aes(n, v - lv, colour = (ttChange < 50 | tsChange < 50) , group = 1)) +
    facet_wrap(~ID, scales='free', ncol = 8) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank()) + 
    labs(x = NULL,y = 'Change in V') -> p2

L3Y600 %>%
  group_by(ID) %>%
  mutate(n = seq_along(time)) %>%
  filter(n > 1 & changed == TRUE & min(delta) > 0) %>%
  filter(ID %in% tail(unique(.$ID), 24)) %>%
  ggplot() + geom_line(aes(n, -relX, colour = (ttChange < 50 | tsChange < 50), group = 1)) +
  facet_wrap(~ID, scales='free', ncol = 8) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank()) + 
  labs(x = NULL,y = 'X position (flipped)') -> p3
gridExtra::grid.arrange(p1, p2, p3, ncol = 1)

L3Y600 %>%
  group_by(ID) %>%
  filter(changed == FALSE) %>%
  mutate(n = seq_along(time),
         deltaT = log((delta/pi)/(1-delta/pi))) %>% 
  filter(ID %in% head(unique(.$ID), 24)) %>%
  ggplot() + geom_line(aes(n, deltaT)) +
  facet_wrap(~ID, scales='free_x', ncol = 8) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank(),
        strip.background = element_blank()) + 
  labs(x = NULL,y = 'delta (unconstrained)') -> p4

gridExtra::grid.arrange(p1, p4, ncol=1)
}

fitARIMAfc{
tsFits <- data.frame()
for(i in 0:5){
  #for(k in 0:1){
    for(l in seq_along(noChange)){
      car <- filter(L3Y600, ID == noChange[l])
      train <- car[1:(nrow(car)-50),]
      test <- car[(nrow(car)-49):nrow(car),]
      train$delta %>% 
        Arima(order = c(i, 0, 0), method = 'ML') %>%
        forecast(h = 50) -> fc
      mse <- mean((fc$mean - test$delta)^2)
      ls <- sum(dnorm(test$delta, fc$mean, sqrt(fc$model$sigma2), log=TRUE))
      tsFits <- rbind(tsFits, 
                      data.frame(var = 'delta', ID = noChange[l], ar = i, mse = mse,
                                 ll = fc$model$loglik, aic = fc$model$aic, bic = fc$model$bic, logscore = ls))
      train$v %>%
        Arima(order = c(i, 1, 0), method = 'ML') %>%
        forecast(h = 50) -> fc
      mse <- mean((fc$mean - test$v)^2)
      ls <- sum(dnorm(test$v, fc$mean, sqrt(fc$model$sigma2), log=TRUE))
      tsFits <- rbind(tsFits, 
                      data.frame(var = 'v', ID = noChange[l], ar = i, mse = mse,
                                 ll = fc$model$loglik, aic = fc$model$aic, bic = fc$model$bic, logscore = ls))
      
        if(l %% 200 == 0){
        print(l)
      }
    }
    print(i)
  #}
}
arFits <- NULL
for(l in seq_along(noChange)){
  car <- filter(L3Y600, ID == noChange[l])
  fitD <- auto.arima(car$delta, max.q = 0)
  fitV <- auto.arima(car$v, max.q = 0)
  arFits <- rbind(arFits,
                  data.frame(v = length(fitV$model$phi), d = length(fitD$model$phi)))
}


for(i in 0:2){
  for(k in 0:2){
    for(l in seq_along(noChange)){
      car <- filter(L3Y600, ID == noChange[l])
      train <- car$v[1:(nrow(car)-50)]
      test <- car$v[(nrow(car)-49):nrow(car)]
      train %>% 
        Arima(order = c(i, 1, k)) %>%
        forecast(h = 50) -> fc
      mse <- mean((fc$mean - test)^2)
      ls <- sum(dnorm(test, fc$mean, sqrt(fc$model$sigma2), log=TRUE))
      tsFits <- rbind(tsFits, 
                      data.frame(var = 'v', ID = noChange[l], ar = i, i = 1, ma = k, mse = mse,
                                 ll = fc$model$loglik, aic = fc$model$aic, bic = fc$model$bic, logscore = ls))
      if(l %% 500 == 0){
        print(l)
      }
    }
    print(paste(i, k))
  }
}

tsFits %>%
  #mutate(model = paste(ar, i, ma)) %>%
  group_by(var, ar) %>%
  summarise(mean(mse), mean(aic), mean(bic), mean(ll), mean(logscore))
# Delta - ARMA(1, 1)
# Velocity - ARIMA(1, 1, 1)

vParam = data.frame()
for(l in seq_along(noChange)){
  car <- filter(L3Y600, ID == noChange[l])
  train <- car[1:(nrow(car)-50),]
  train$v %>%
    Arima(order = c(2, 1, 1)) -> fc
  vParam = rbind(vParam, data.frame(ar1 = fc$coef[1], ar2 = fc$coef[2], ma = fc$coef[3], sig2 = fc$sigma2))
}
vParam %>%
  gather(var, est) %>%
  group_by(var) %>% 
  summarise(med = median(est), mean = mean(est))
}

simulateData{
dar = 0.67
dma = 0.29
var = 0.68
vma = 0.53
sigSqD = 5.93e-05
sigSqV = 6.62e-04

T = 500

v <- 3 + arima.sim(list(order = c(1, 1, 1), ar= var, ma = vma), n = T-1, n.start = 100, sd = sqrt(sigSqV))
d <- pi/2 + arima.sim(list(ar = dar, ma = dma), n = T, n.start = 100, sd = sqrt(sigSqD))
xdiff <- v * cos(d)
ydiff <- v * sin(d)
x <- cumsum(xdiff)
y <- 600 + cumsum(ydiff)
pos <- data.frame(x = x, y = y, v = v, delta = d, method = 'simulated')
d <- pos$delta[2:T]
vdiff <- pos$v[2:T] - pos$v[1:(T-1)]
data <- cbind(vdiff, d)

L3Y600 %>% 
  filter(ID == noChange[sample(1:length(noChange), 1)]) %>%
  ggplot() + geom_path(aes(xRel, y)) + 
  geom_path(data = pos, aes(x, y), colour = 'red')
}

fitARIMAvb{

lMean <- c(0, 0, 0, 0, 0, 0)
lVar <- diag(0.5, 6)
hyper <- c(2, 1e-3, 2, 1e-4, 1, 1, 1, 1, 1, 1, 1, 1)
lambda <- as.matrix(c(lMean, lVar), ncol=1)

trueV <- data.frame(true = c(sigSqV, sigSqD, dar, var, dma, vma),
                    var = c('sigSqV', 'sigSqD', 'arD', 'arV', 'maD', 'maV'))



fit <- carsArima(data, lambda, hyper, 10, 2000, 0.15, 0.9, 0.99, 0.01)

vbDensity(fit, 
          c(rep('exp', 2), rep('sigmoid', 4)),
          c('sigSqV', 'sigSqD', 'arV', 'maV', 'arD', 'maD')) -> densities

densities %>%
  group_by(var) %>%
  summarise(mean = sum(support * density)*(support[2] - support[1]),
            map = support[which.max(density)]) %>%
  merge(trueV, by = 'var')

densities %>%
  ggplot() + geom_line(aes(support, density)) +  
    geom_vline(data = trueV, aes(xintercept=true), colour = 'red') + 
    facet_wrap(~var, scales='free')

VBSims = data.frame()
for(i in 1:100){
  v <- 3 + arima.sim(list(order = c(1, 1, 1), ar= var, ma = vma), n = T-1, n.start = 100, sd = sqrt(sigSqV))
  d <- pi/2 + arima.sim(list(ar = dar, ma = dma), n = T, n.start = 100, sd = sqrt(sigSqD))
  xdiff <- v * cos(d)
  ydiff <- v * sin(d)
  x <- cumsum(xdiff)
  y <- 600 + cumsum(ydiff)
  pos <- data.frame(x = x, y = y, v = v, delta = d, method = 'simulated')
  d <- pos$delta[2:T]
  vdiff <- pos$v[2:T] - pos$v[1:(T-1)]
  data <- cbind(vdiff, d)
  fit <- carsVB(data, lambda, 10, 2000, 0.15, 0.9, 0.99, 0.01, hyper=hyper, model = arimaDeriv)
  vbDensity(fit, 
            c(rep('exp', 2), rep('sigmoid', 4)),
            c('sigSqV', 'sigSqD', 'arV', 'maV', 'arD', 'maD')) -> densities
  densities %>%
    group_by(var) %>%
    summarise(mean = sum(support * density)*(support[2] - support[1]),
              map = support[which.max(density)],
              iter = i) -> df
  VBSims <- rbind(VBSims, df)
  print(i)
}

VBSims %>%
  gather(metric, value, -var, -iter) %>%
  ggplot() + geom_boxplot(aes(metric, value)) + 
  geom_hline(data=trueV, aes(yintercept = true), colour = 'red') + facet_wrap(~var, scales = 'free')
}

laneSwitching{


sigSqV = 0.00059
sigSqD = 0.0001
zeta = 0.05
M3 = 0
M2 = -10
M4 = 10
v = 3
d = pi/2
x = 0
y = 0
T = 500
pos = NULL
var1 = 1.17
var2 = -0.52
vma = 0.14
v <- 3 + arima.sim(list(order = c(2, 1, 1), ar= c(var1, var2), ma = vma), n = T-1, n.start = 100, sd = sqrt(sigSqV))
for(t in 1:T){
  if(t < 60){
    M = M3
  } else if(t < 160) {
    M = M2
  } else if(t < 300) {
    M = M3
  } else {
    M = M4
  }
  if(is.na(acos((x-M)/(50*v[t])))) break
  d = acos((M - x)/(30*v[t])) + rnorm(1, 0, sqrt(sigSqD))
  x = x + v[t] * cos(d)
  y = y + v[t] * sin(d)
  pos = rbind(pos, data.frame(x=x, y=y, d=d, v=v[t], t=t, goal=M))
}
ggplot(pos) + geom_path(aes(x, y, group = 1, colour = factor(goal)))


L3Y600 %>% 
  filter(changed == TRUE) %>%
  .$ID %>%
  unique() -> changedIDs

i = sample(changedIDs, 1)

L3Y600 %>% 
  filter(ID == i) -> carChange
carChange %>%
  ggplot() + geom_path(aes(relX, y)) -> p1
ll <- matrix(0, nrow(carChange)-1, 3)
for(i in 2:nrow(carChange)){
  ll[i-1, 1] = dnorm(carChange$delta[i], acos(-carChange$relX[i-1]/(30*carChange$v[i])), sqrt(0.0005))
  ll[i-1, 2] = dnorm(carChange$delta[i], acos((-10-carChange$relX[i-1])/(30*carChange$v[i])), sqrt(0.0005))
  ll[i-1, 3] = dnorm(carChange$delta[i], acos((10-carChange$relX[i-1])/(30*carChange$v[i])), sqrt(0.0005))
}
ll <- as.data.frame(ll)
colnames(ll) <- c('M3', 'M2', 'M4')
ll$t <- 2:nrow(carChange)
ll %>% 
  mutate(l2 = M2 / (M2 + M3 + M4), l3 = M3 / (M2 + M3 + M4), l4 = M4 / (M2 + M3 + M4)) %>%
  select(l2, l3, l4, t) %>%
  gather(lane, ll, -t) %>%
  ggplot() + geom_path(aes(ll, t, colour = lane)) -> p2

gridExtra::grid.arrange(p1, p2, ncol=2)

}

AR1{
i <- sample(noChange, 1)
L3Y600 %>%
  filter(ID == i[1]) %>%
  select(v, delta) -> car

data <- cbind(car$v[2:nrow(car)] - car$v[1:(nrow(car)-1)], car$delta[2:nrow(car)])
mu <- c(0, 0, 0, 0)
sd <- c(1, 1, 1, 1)
lambda <- matrix(c(mu, diag(sd)), nrow=20)
hyper <- c(2, 0.0002, 2, 0.00002, 1, 1, 1, 1)

fit <- carsVB(data, lambda, hyper=hyper, S=5, maxIter=5000, alpha=0.01, beta1=0.9, beta2=0.99,
              dimTheta=4, model = ar1Deriv, threshold=0.01)

density <- vbDensity(fit, c(rep('exp', 2), rep('stretchedSigmoid', 2)), c('sigma^2[V]', 'sigma^2[D]', 'phi[V]', 'phi[D]'))
density %>%
  ggplot() + geom_line(aes(support, density)) +
  facet_wrap(~var, scales='free') + 
  theme_bw() + 
  theme(strip.background = element_blank())

density %>% 
  group_by(var) %>%
  summarise(map = support[which.max(density)])
}

heterogeneity{

N <- 20
idSubset <- sort(sample(noChange, N))
  
L3Y600 %>%
  filter(ID %in% idSubset) %>%
  select(v, delta, ID, class) -> carhetero

carhetero %>%
  group_by(ID) %>%
  mutate(n = seq_along(v)) %>%
  filter(n == 1) %>%
  .$class %>%
  table()

carhetero %>%
  group_by(ID) %>%
  mutate(n = seq_along(v)) %>%
  filter(n > 1) %>%
  summarise(n =n()) %>%
  .$n %>%
  cumsum() -> obsSum

carhetero %>%
  group_by(ID) %>%
  mutate(n = seq_along(v),
         vlag = ifelse(n == 1, 0, lag(v)),
         vdiff = v - vlag) %>%
  filter(n > 1) %>%
  ungroup() %>%
  select(vdiff, delta) %>%
  as.matrix() -> data

    
dim = 6
dim2 = 0.5 * dim * (dim + 1)
varSeq <- 0.1
for(i in 1:(dim - 1)){
  varSeq <- c(varSeq, rep(0, i), 0.1)
}
lambda <- matrix(c(rep(0, dim*(N+1)), rep(varSeq, N+1)), ncol=1)
hyperMean <- c(-5, 5, rep(0, dim-2))
hyperVar <- diag(5, dim)
hyperLinv <- solve(chol(hyperVar))
Linv <- hyperLinv

heirFit <- carsVB(data = data,
                  lambda = lambda,
                  S = 15,
                  maxIter = 5000,
                  model = heirAr2Deriv,  ## Remember this line
                  dimTheta = dim*(N+1),
                  dimLambda = length(lambda),
                  hyperMean = hyperMean,
                  hyperLinv = hyperLinv,
                  Linv = Linv,
                  obsSum = obsSum,
                  threshold = 0.25)


heirList <- list(list(mean = heirFit[1:dim], U = heirFit[(N+1)*dim+1:dim2]))
for(i in 2:(N+1)){
  heirList[[i]] <- list(mean = heirFit[(i-1)*dim + 1:dim], U = heirFit[(N+1)*dim+dim2*(i-1) + 1:dim2]) 
}


compareCar <- sample(1:N, 1)
compareModels(heirList, idSubset, compareCar)

maps <- data.frame()
for(i in 1:N){
  dens <- vbDensity(heirList[[i+1]],
                    c(rep('exp', 2), rep('stretchedSigmoid', 4)),
                    c('sigma^2[V]', 'sigma^2[D]', 'phi1[V]', 'phi1[D]', 'phi2[V]', 'phi2[D]'))
  map <- dens %>% group_by(var) %>% summarise(map = support[which.max(density)])
  df <- data.frame(ID = idSubset[i], map, 
                   class = carhetero %>%
                     filter(ID == idSubset[i]) %>%
                     head(1) %>%
                     .$class)
  maps <- rbind(maps, df)
}

maps %>%
  ggplot() + geom_boxplot(aes(x = class, y = map, fill = class)) + facet_wrap(~var, scales = 'free')

maps %>%
  spread(var, map) %>%
  select(-ID) %>%
  rename(sigV = `sigma^2[V]`, sigD = `sigma^2[D]`, phiV = `phi1[V]`, phiD = `phi1[D]`, phi2V = `phi2[V]`, phi2D = `phi2[D]`) %>%
  GGally::ggpairs(aes(colour = class, alpha = 0.8))



heirDensity <- vbDensity(heirList[[2]], 
                         c(rep('exp', 2), rep('stretchedSigmoid', 4)),
                         c('sigma^2[V]', 'sigma^2[D]', 'phi1[V]', 'phi1[D]', 'phi2[V]', 'phi2[D]'))

heirDensity$method <- 'heirarchical model - one car'
globalDensity <- vbDensity(heirList[[1]],
                           c(rep('exp', 2), rep('stretchedSigmoid', 4)),
                           c('sigma^2[V]', 'sigma^2[D]', 'phi1[V]', 'phi1[D]', 'phi2[V]', 'phi2[D]'))
globalDensity$method <- 'heirarchical model - global'
heirDensity %>%
  rbind(globalDensity) %>%
  ggplot() + geom_line(aes(support, density)) + 
  facet_wrap(method~var, scales = 'free', ncol =dim) + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  labs(x = NULL, y = NULL)
}
