library(tidyverse)
library(GGally)
library(forecast)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
source('mixtureMCMC.R') # Hierarchical Mixture Model with MCMC- no noise
#source('carsVBfuns.R') # Functions to assist VB 
#source(greta.R) # Contains Hamiltonian MCMC for an AR(1) with the full states / theta 
#sourceCpp('heirarchical.cpp') # AR(p) and Heir-AR(p) VB models
#sourceCpp('arvdMCMC.cpp') #AR(p) MCMC
#sourceCpp('hamiltonianPF.cpp') # Hamiltonian MCMC with Particle Filter
#sourceCpp('AR1PMMH.cpp) # Random Walk PMMH for AR1
#sourceCpp('Ar1ParticleFilter.cpp) # AR1 VB with a Particle Filter, and with full states

## Gibbs One at a time MCMC is easier to deal with things even if inefficient
## Paper from Tas
## Blei VB Consistency

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
    fit <- carsVB(data, lambda, 10, 2000, 0.15, 0.9, 0.99, 0.01, hyper=hyper, model = arimaDeriv)$lambda
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

AR{
lags <- 1
ID <- sample(noChange, 1)
data <- sampleCars(L3Y600, ID)$data
mu <- rep(0, 2 + 2 * lags)
sd <- rep(0.2, 2 + 2 * lags)
lambda <- matrix(c(mu, diag(sd)), ncol=1)
hyper <- c(2, 0.0002, 2, 0.00002, rep(c(0, 1), 2 * lags))

fit <- carsVB(data = data, 
              lambda = lambda,
              hyper = hyper,
              dimTheta = 2 + 2*lags,
              model = arDeriv,
              lags = lags)$lambda

transform <- c(rep('exp', 2), rep('identity', 2 * lags))
names <- c('sigma2V', 'sigma2D')
for(i in 1:lags){
  names <- c(names, paste0('phi', i, 'V'), paste0('phi', i, 'D'))
}

fitList <- list(mean = fit[1:(2+2*lags)], U = fit[(3+2*lags):length(fit)])

density <- vbDensity(fitList, transform, names)
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
  
N <- 200
lags <- 2
diag <- TRUE

idSubset <- sample(noChange, N)
data <- sampleCars(L3Y600, idSubset)
obsSum <- data$obsSum
class <- data$class
data <- data$data

dim <- 2 + 2 * lags
dim2 <- 0.5 * dim * (dim + 1)
varSeq <- 0.1
for(i in 1:(dim - 1)){
  varSeq <- c(varSeq, rep(0, i), 0.1)
}
lambda <- matrix(c(rep(0, dim*(N+1)), rep(varSeq, N+1)), ncol=1)
hyperMean <- c(-7, 7, rep(0, dim-2))
hyperVar <- diag(5, dim)
hyperLinv <- solve(chol(hyperVar))
Var <- diag(0.1^2, dim)
Linv <- solve(chol(Var))

heirFit <- carsVB(data = data,
                  lambda = lambda,
                  S = 10,
                  model = heirArDeriv,
                  dimTheta = dim*(N+1),
                  dimLambda = length(lambda),
                  hyperMean = hyperMean,
                  hyperLinv = hyperLinv,
                  Linv = Linv,
                  lags = lags,
                  obsSum = obsSum,
                  threshold = 1,
                  diag = diag)$lambda
  
  
heirList <- list(list(mean = heirFit[1:dim], U = heirFit[(N+1)*dim+1:dim2]))
for(i in 2:(N+1)){
  heirList[[i]] <- list(mean = heirFit[(i-1)*dim + 1:dim], U = heirFit[(N+1)*dim+dim2*(i-1) + 1:dim2]) 
}
  
transform <- c(rep('exp', 2), rep('identity', 2 * lags))
names <- c('sigma2V', 'sigma2D')
for(i in 1:lags){
  names <- c(names, paste0('phi', i, 'V'), paste0('phi', i, 'D'))
}

compareCar <- sample(1:N, 1)
compareModels(heirList, idSubset, compareCar, lags, transform, names, TRUE)
  
maps <- data.frame()
for(i in 1:N){
  dens <- vbDensity(heirList[[i+1]],
                    transform, 
                    names)
  map <- dens %>% group_by(var) %>% summarise(map = support[which.max(density)])
  df <- data.frame(ID = idSubset[i], map, 
                   class = rep(class[i], dim))
  maps <- rbind(maps, df)
}
  
maps %>%
  ggplot() + geom_boxplot(aes(x = class, y = map, fill = class)) + facet_wrap(~var, scales = 'free')
  
maps %>%
  spread(var, map) %>%
  select(class, sigma2V, phi1V, phi2V, sigma2D, phi1D, phi2D) %>%
  ggparcoord(columns = 2:7,
             mapping = aes(colour = factor(class)))

maps %>%
  spread(var, map) %>%
  select(-ID) %>%
  ggpairs(mapping = aes(colour = class),
          diag = list(continuous = wrap('densityDiag', alpha = 0.75))) + 
  theme_bw()
  
  
maps %>%
  spread(var, map) %>%
  select(sigma2V, phi1V, phi2V) %>%
  kmeans(centers = 2) %>%
  .$cluster -> Vmeans

maps %>%
  .$ID %>%
  unique() %>%
  cbind(Vmeans) %>%
  as.data.frame() -> Vmeans
colnames(Vmeans) <- c('ID', 'cluster')

maps %>%
  spread(var, map) %>%
  cbind(cluster = factor(Vmeans$cluster)) %>%
  ggplot() + geom_path(aes(sigma2V, phi1V, group = 1, colour = factor(cluster)))
  ggpairs(columns = 3:8,
          mapping = aes(colour = cluster),
          diag = list(continuous = wrap('densityDiag', alpha = 0.75))) + 
  theme_bw()
      
L3Y600 %>%
  filter(ID %in% idSubset) %>%
  left_join(Vmeans, by = 'ID', copy = TRUE) %>%
  ggplot() + geom_line(aes(y, v, group = ID, colour = factor(cluster)))

  
heirDensity <- vbDensity(heirList[[sample(N, 1)+1]], transform, names)
heirDensity$method <- 'heirarchical model - one car'

globalDensity <- vbDensity(heirList[[1]], transform, names)
globalDensity$method <- 'heirarchical model - global'

heirDensity %>%
  rbind(globalDensity) %>%
  ggplot() + geom_line(aes(support, density)) + 
  facet_wrap(method~var, scales = 'free', ncol =dim) + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  labs(x = NULL, y = NULL)
}

arUpdate{
  
lags <- 1
N <- 10
  
transform <- c(rep('exp', 2), rep('identity', 2 * lags))
names <- c('sigma2V', 'sigma2D')
for(i in 1:lags){
  names <- c(names, paste0('phi', i, '[V]'), paste0('phi', i, '[D]'))
}
  
idSubset <- sample(noChange, N)
data <- sampleCars(L3Y600, idSubset)
obsSum <- data$obsSum
data <- data$data
  
dim <- 2 + 2 * lags
dim2 <- 0.5 * dim * (dim + 1)
varSeq <- 0.1
for(i in 1:(dim - 1)){
  varSeq <- c(varSeq, rep(0, i), 0.1)
}
  
lambda <- matrix(c(rep(0, dim*(N+1)), rep(varSeq, N+1)), ncol=1)
hyperMean <- c(-7, 7, rep(0, dim-2))
hyperVar <- diag(5, dim)
hyperLinv <- solve(chol(hyperVar))
Linv <- hyperLinv

heirFit <- carsVB(data = data,
                  lambda = lambda,
                  S = 15,
                  model = heirArDeriv,
                  dimTheta = dim*(N+1),
                  dimLambda = length(lambda),
                  hyperMean = hyperMean,
                  hyperLinv = hyperLinv,
                  Linv = Linv,
                  obsSum = obsSum,
                  lags = lags,
                  threshold = 0.25)$lambda

heirGlobal <- list(mean = heirFit[1:dim], U = heirFit[(N+1)*dim+1:dim2])

id <- sample(noChange, 1)
data <- sampleCars(L3Y600, id)$data
  
lambdaUpdate <- heirGlobal$mean
for(i in 1:dim){
  lambdaUpdate <- c(lambdaUpdate, heirGlobal$U[sum(0:(i-1)) + 1:i], rep(0, dim-i))
}
lambdaUpdate <- matrix(lambdaUpdate, ncol=1)
fit <- updateVB(data, lambdaUpdate, hyper, stepsize = 10, lags = lags)
fitUpdate <- list(mean = fit[1:dim], U = fit[(dim+1):length(fit)])

compareModels(fitUpdate, id, 1, lags, transform, names, heir = FALSE, S = 25)
  
upDensity <- vbDensity(fitUpdate, transform, names)
upDensity$method <- 'update'

globalDensity <- vbDensity(heirGlobal, transform, names)
globalDensity$method <- 'global'

upDensity %>%
  rbind(globalDensity) %>%
  ggplot() + geom_line(aes(support, density)) + 
  facet_wrap(method~var, scales = 'free', ncol =dim) + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  labs(x = NULL, y = NULL)

mu <- rep(0, 2 + 2 * lags)
sd <- rep(0.2, 2 + 2 * lags)
lambda <- matrix(c(mu, diag(sd)), ncol=1)
hyper <- c(2, 0.0002, 2, 0.00002, rep(c(0, 1), 2 * lags))

fitSingle <- carsVB(data, lambda, hyper=hyper, S=5, maxIter=5000, alpha=0.01, beta1=0.9, beta2=0.99,
              dimTheta= 2 + 2*lags, model = arDeriv, lags = lags, threshold=0.01)$lambda

fitSL <- list(mean = fitSingle[1:dim], U = fitSingle[(dim+1):length(fit)])
sDensity <- vbDensity(fitSL, transform, names)
sDensity$method <- 'single'

upDensity %>%
  rbind(sDensity) %>%
  rbind(globalDensity) %>%
  ggplot() + geom_line(aes(support, density)) + 
  facet_wrap(method~var, scales = 'free', ncol =dim) + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  labs(x = NULL, y = NULL)

}

ar1NoiseVB{
  lags <- 1
  T <- 200
  id <- sample(noChange, 1)
  L3Y600 %>% 
    filter(ID == id) %>%
    .$v %>%
    head(1) -> initV
  L3Y600 %>%
    filter(ID == id) %>%
    select(relX, y) %>%
    mutate(n = seq_along(y)) %>%
    filter(n <= T) %>%
    select(relX, y) %>%
    as.matrix() -> data 
  mu <- c(-8, -8, -5, 0, 0, 0)
  sd <- rep(0.2, 4 + 2 * lags)
  #states <- rep(rep(c(0, 0.1), rep(nrow(data), 2)), 2)
  lambda <- matrix(c(mu, diag(sd)), ncol=1)# states), ncol=1)
  hyper <- c(rep(c(2, 0.0002), 3), 1, 1)#, rep(c(0, 1), 2 * lags))
  
  fit <- carsVB(data = data, 
                lambda = lambda,
                hyper = hyper,
                dimTheta = 6,#4 + 2*lags + 2 * nrow(data),
                dimLambda = 42,# + 4 * nrow(data),
                model = arPFDeriv,
                maxIter = 5000,
                S = 100,
                P = 100,
                initV = initV,
                threshold = 0.25,
                lags = 1)
  qplot(1:fit$iter, fit$LB, geom = 'line')
  
  #fitTheta <- fit[1:42]
  #fitA <- fit[43:442]
  #fitD <- fit[443:842]
  L3Y600 %>%
    filter(ID == id) %>%
    mutate(n = seq_along(v)) %>%
    filter(n <= 200) %>%
    cbind(fit = fitD[seq(1, 399, 2)],
          sd = abs(fitD[seq(2, 400, 2)])) %>%
    mutate(upper = fit + 1.96 * sd,
          lower = fit - 1.96 * sd) %>%
    ggplot() + geom_ribbon(aes(x=n, ymin = lower, ymax = upper), fill = 'grey70') +
    geom_path(aes(n, delta-pi/2)) + geom_path(aes(n, fit), colour = 'red') 
  
  L3Y600 %>%
    filter(ID == id) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n==1, 0, lag(v)),
           a = v - vl) %>%
    filter(n > 1 & n <= 201) %>%
    cbind(fit = fitA[seq(1, 399, 2)],
          sd = abs(fitA[seq(2, 400, 2)])) %>%
    mutate(upper = fit + 1.96 * sd,
           lower = fit - 1.96 * sd) %>%
    ggplot() + geom_ribbon(aes(x=n, ymin = lower, ymax = upper), fill = 'grey70') +
    geom_path(aes(n, a)) + geom_path(aes(n, fit), colour = 'red')
    
  
  
  fitList <- list(mean = fit$lambda[1:(4+2*lags)], U = fit$lambda[(5+2*lags):length(fit$lambda)])
  
  transform <- c(rep('exp', 4), rep('stretchedSigmoid', 2 * lags))
  names <- c('sigma2V', 'sigma2D', 'sigma2X', 'sigma2Y')
  for(i in 1:lags){
    names <- c(names, paste0('phi', i, 'V'), paste0('phi', i, 'D'))
  }

  
  density <- vbDensity(fitList, transform, names)
  density %>%
    ggplot() + geom_line(aes(support, density)) +
    facet_wrap(~var, scales='free') + 
    theme_bw() + 
    theme(strip.background = element_blank())
  
  density %>% 
    group_by(var) %>%
    summarise(map = support[which.max(density)])
}

individualNIG{
  
N <- 200
lags <- 2
idSubset <- sample(noChange, N)

names <- c('sigma2V', 'sigma2D')
for(i in 1:lags){
  names <- c(names, paste0('phi', i, 'V'))
}
for(i in 1:lags){
  names <- c(names, paste0('phi', i, 'D'))
}
  
reps <- 10000
MCMCres <- NULL
hyper <- c(2, 2e-04, 2, 2e-05, rep(c(0, 1), 2*lags))
for(i in 1:N){
  data %>%
    filter(ID == idSubset[i]) %>%
    select(vd, d) %>%
    as.matrix() -> tempDat
  
  draws <- ARVD_MCMC(tempDat, hyper, reps, 2)
  maps <- NULL
  for(k in 1:(2 + 2*lags)){
    dens <- density(draws[(reps/2+1):reps,k])
    maps <- c(maps, dens$x[which.max(dens$y)])
  }
  
  MCMCres <- rbind(MCMCres, 
                   data.frame(ID = idSubset[i], var = names, mean = maps))
}

MCMCres %>%
  spread(var, mean) %>%
  select(-ID) %>%
  kmeans(centers = 3) %>%
  .$cluster -> Vmeans

MCMCres %>%
  .$ID %>%
  unique() %>%
  cbind(Vmeans) %>%
  as.data.frame() -> Vmeans
colnames(Vmeans) <- c('ID', 'cluster')

MCMCres %>%
  spread(var, mean) %>%
  cbind(cluster = factor(Vmeans$cluster)) %>%
  ggpairs(columns = c(2, 4, 3, 5:8),
          mapping = aes(colour = cluster),
          diag = list(continuous = wrap('densityDiag', alpha = 0.75))) + 
  theme_bw()



  
  
  
  
  
}

ar1NoiseMCMC{
T <- 200
idSubset <- sample(noChange, 1)

L3Y600 %>%
  filter(ID == idSubset) %>%
  .$v %>%
  head(1) -> initV
L3Y600 %>%
  filter(ID == idSubset) %>%
  select(relX, y) -> tempDat

xM <- tempDat$relX[1:T]
yM <- tempDat$y[1:T]

data <- list(T = T,
             xM = xM,
             yM = yM,
             initV = initV)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

initf <- function(){
  list(arD = 0.8, arV = 0.7, sigSqXInv = 100, sigSqYInv = 1, sigSqVInv = 1000, sigSqDInv = 1000)
}
initL <- lapply(1:4, function(x) initf())

fit <- stan('singleCar.stan', data = data, init = initL)

hyper <- c(1, 0.01, 1, 0.01, 1, 0.01, 1, 1, 1, 1, 1, 1)
stepSize <- c(1e-6, 1e-7, 1e-6, 1e-3, 1e-3, 1e-3)
initial <- c(8.5e-04, 7.6e-05, 1e-3, 1, 0.75, 0.8)

fitPMMH <- PMMH(cbind(xM, yM), 100, 2000, hyper, stepSize, initial, initV)

fitHamPMMH <- hamiltonianPF(data = cbind(xM, yM),
                            hyper = hyper,
                            initV = initV,
                            P = 50,
                            reps = 1, 
                            initTheta = matrix(initial, ncol=1),
                            M = rep(1, 6),
                            L = 20,
                            epsilon = 0.01)


}

Mixture{
  N <- 125
  T <- 250
  idSubset <- sample(noChange)
  data <- sampleCars(L3Y600, idSubset)
  i <- 1
  try <- 1
  d <- matrix(0, T, N)
  v <- matrix(0, T, N)
  while(i <= N){
    L3Y600 %>%
      filter(ID == idSubset[try]) %>%
      mutate(n = seq_along(v),
             vl = ifelse(n == 1, 0, lag(v)),
             vdiff = v - lag(v)) %>%
             filter(n > 1) %>% 
             select(delta, vdiff) -> tempData
    if(nrow(tempData) >= T){
      d[,i] <- tempData$delta[1:T]
      v[,i] <- tempData$vdiff[1:T]
      i <- i + 1
    }
    try <- try + 1
  }
  data <- list(N = N,
               T = T,
               v = v,
               d = d)

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
  fit <- stan('mixture.stan', data = data, iter = 300, chains = 1, control = list(adapt_delta = 0.95, max_treedepth = 20))
  
  M <- 500
  idSubset <- sample(noChange, M)
  data <- list()
  for(i in 1:M){
    L3Y600 %>%
      filter(ID == idSubset[i]) %>%
      mutate(n = seq_along(v),
             vl = ifelse(n == 1, 0, lag(v)),
             v = v - lag(v),
             d = delta - pi/2) %>%
      filter(n > 1) %>% 
      select(v , d) %>%
      as.matrix() -> data[[i]]
  }
  #saveRDS(data, 'rmixdata.RDS')
  N <- 75
  dataSub <- list()
  for(i in 1:N){
    dataSub[[i]] <- data[[i]]
  }
  reps <- 20000
  mixDraws <- mixtureMCMC(dataSub, reps)

  mapGroup <- NULL
  for(i in 1:N){
    kDraws <- mixDraws[[i+1]]$k[(reps/2+1):reps]
    mapGroup <- rbind(mapGroup, data.frame(group = mean(kDraws), ID = idSubset[i]))
  }
  ggplot(mapGroup) + geom_histogram(aes(group)) + theme_bw()
  
  mixDraws[[1]]$mean1 %>% 
    cbind(iter = 1:reps) %>%
    as.data.frame() -> mean1
  mixDraws[[1]]$mean2 %>%
    cbind(iter = 1:reps) %>%
    as.data.frame()  -> mean2
  
  colnames(mean1) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter')
  colnames(mean2) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter')
  
  mean1 %>%
    rbind(mean2) %>%
    cbind(method = rep(c('mu1', 'mu2'), rep(reps, 2))) %>%
    gather(var, draw, -iter, -method) %>%
    filter(iter > 5000) %>%
    ggplot() + geom_line(aes(iter, draw)) + 
    facet_grid(var ~ method, scales = 'free') + theme_bw()
  
  mean1 %>%
    filter(iter > 5000) %>%
    select(-iter) %>%
    GGally::ggpairs() + 
    theme_bw() + labs(title = 'mu_1')
  
  mean2 %>%
    filter(iter > 5000) %>%
    select(-iter) %>%
    GGally::ggpairs() + 
    theme_bw() + labs(title = 'mu_2')
}