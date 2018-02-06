########## TO DO #############
# Write more things.
# Check lane change car results (ongoing on slurm)
# Built the Neural Network (find why loss goes to nan)

library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(GGally)

setup {
id <- readRDS('carsID.RDS')

read.csv('carsAug.csv') %>%
  select(ID, relX, dist, relv, reldelta) %>%
  ungroup() -> carsAug
colnames(carsAug) <- c('ID', 'x', 'y', 'v', 'delta')

carsAug %>% 
  group_by(ID) %>%
  filter(min(v) != 0 & !ID %in% id$splinesID) %>%
  .$ID %>%
  unique() -> noStop

read.csv('carsChanged.csv') %>%
  select(ID, relX, dist, relv, reldelta) %>%
  ungroup() -> carsChanged
colnames(carsChanged) <- c('ID', 'x', 'y', 'v', 'delta')

carsChanged %>% 
  group_by(ID) %>%
  filter(min(v) != 0) %>%
  .$ID %>%
  unique() -> noStopChanged

dataChanged <- list()
for(i in 1:500){
  carsChanged %>%
    filter(ID == noStopChanged[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           a = v - lag(v),
           d = delta - pi/2,
           xl = ifelse(n == 1, 0, lag(x)),
           dx = x - xl,
           yl = ifelse(n == 1, 0, lag(y)),
           dy = y - yl) %>%
    filter(n > 1 & n <= 501) %>% 
    select(a , d, vl, dx, dy) %>%
    as.matrix() -> dataChanged[[i]]
}
dataChanged[[501]] <- noStopChanged
saveRDS(dataChanged, 'fcDataChanged.RDS')
}

hierMixtureModel{

set.seed(1)
N <- 2000
idSubset <- sample(noStop, N)
data <- list()
for(i in 1:N){
  carsAug %>%
    filter(ID == idSubset[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           a = v - lag(v),
           d = delta - pi/2) %>%
    filter(n > 1 & n <= 501) %>% 
    select(a , d) %>%
    as.matrix() -> data[[i]]
}
saveRDS(data, 'MCMCData.RDS')

reps <- 80000
K <- 6
thin <- 10
burn <- 0.9

draws <- list(list())
hyper <- list()
for(k in 1:K){
  hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(c(5, 5, 10, 10, 10, 10))), v = 6, scale = diag(1, 6))
  draws[[1]][[k]] <- list(mean = c(-5, -5, 0, 0, 0, 0), varInv = diag(10, 6))
}
draw[[1]]$pi <- rep(1/K, K)
hyper$alpha <- rep(1, K)
for(i in 1:N){
  draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:K, 1), pi = rep(1/K, K))
}



mixDraws <- mixtureMCMC(data, reps, draws, hyper, thin, K, 'gaussian', 0.01)

saveRDS(mixDraws, 'mixN2000K6.RDS')

mapGroup <- NULL
for(i in 1:N){
  kDraws <- mixDraws$draws[[i+1]]$k[(burn*reps/thin+1):(reps/thin)]
  mode <- which.max(table(c(1:K, kDraws)))
  mapGroup <- rbind(mapGroup, data.frame(group = c(mean(kDraws), mode), 
                                         method = c('mean', 'mode'),
                                         ID =  idSubset[i]))
}
ggplot(mapGroup) + geom_histogram(aes(group), binwidth = 0.1) + theme_bw() + facet_wrap(~method)

muK <- NULL
for(i in 1:K){
  mixDraws$draws[[1]][[i]]$mean[2:(reps/thin),] %>%
    cbind(iter = seq(2*thin, reps, thin)) %>%
    as.data.frame() %>%
    mutate(group = i)  -> temp
  colnames(temp) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter', 'group')
  muK <- rbind(muK, temp)
}

muK %>%
  gather(var, draw, -iter, -group) %>%
  mutate(var = factor(var, levels = c('log_sigSq_eps', 'phi1', 'phi2', 'log_sigSq_eta', 'gamma1', 'gamma2'))) %>%
  filter(iter > 2000 & iter %% 20 == 0 ) %>%
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_grid(var ~ group, scales = 'free') + theme_bw() + labs(title = 'mean') -> p1


sdK <- NULL
for(i in 1:K){
  mixDraws$draws[[1]][[i]]$varInv[,,2:(reps/thin)] %>%
    apply(3, function(x) sqrt(diag(solve(x)))) %>%
    t() %>%
    as.data.frame() %>%
    cbind(iter = seq(2*thin, reps, thin)) %>%
    mutate(group = i) -> temp
  colnames(temp) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter', 'group')
  sdK <- rbind(sdK, temp)
}

sdK %>%
  gather(var, draw, -iter, -group) %>%
  mutate(var = factor(var, levels = c('log_sigSq_eps', 'phi1', 'phi2', 'log_sigSq_eta', 'gamma1', 'gamma2'))) %>%
  filter(iter > 2000 & iter %% 20 == 0 ) %>%
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_grid(var ~ group, scales = 'free') + theme_bw() + labs(title = 'standard deviation') -> p2

gridExtra::grid.arrange(p1, p2, ncol = 2)

for(i in 1:K){
  mixDraws$draws[[1]][[i]]$varInv[,,(reps/(2*thin)+1):(reps/thin)] %>% 
    apply(3, function(x) cov2cor(solve(x))) %>%
    t() %>%
    colMeans() %>%
    matrix(6) -> mat
  colnames(mat) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2')
  rownames(mat) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2') 
  corrplot::corrplot(mat) 
}



densities <- NULL
support <- data.frame(seq(exp(-13), exp(-4.5), length.out = 1000),
                      seq(exp(-15), exp(-7.6), length.out = 1000),
                      seq(-0.5, 2, length.out = 1000),
                      seq(-1, 0.8, length.out = 1000),
                      seq(0, 2, length.out = 1000),
                      seq(-1, 0.5, length.out = 1000))
vars <- c('sigSq_eps', 'sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2')
index <- (burn*reps/thin+1):(reps/thin)
for(k in 1:K){
  mat <- matrix(0, 1000, 6)
  for(i in seq_along(index)){
    meanvec <- mixDraws$draws[[1]][[k]]$mean[i,]
    sdvec <- sqrt(diag(solve(mixDraws$draws[[1]][[k]]$varInv[,,i])))
    w <- 0
    for(j in 1:N){
      w <- w + mixDraws$draws[[j+1]]$pi[i, k] / N
    }
    for(j in 1:2){
      mat[,j] <- mat[,j] + w * dlnorm(support[,j], meanvec[j], sdvec[j]) / length(index) 
    } 
    for(j in 3:6){
      mat[,j] <- mat[,j] + w * dnorm(support[,j], meanvec[j], sdvec[j]) / length(index) 
    }
  }
  densities <- rbind(densities,
                     data.frame(dens = c(mat),
                                support = unlist(c(support)),
                                var = rep(vars, rep(1000, 6)),
                                group = k))
}




densities %>%
  mutate(var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2',
                                      'sigSq_eta', 'gamma1', 'gamma2'))) %>%
  group_by(var, support) %>%
  summarise(dens = sum(dens)) -> densMixMod

  ggplot(densMixMod) + geom_line(aes(support, dens)) +
  geom_line(data=densities, aes(support, dens, colour = factor(group))) +
  facet_wrap(~var, scales = 'free', ncol = 6) + 
  labs(title = 'Hierarchical Model Prior', x = NULL, y = NULL) + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) -> p1

  means <- rep(0, 6*6)
  varinv <- matrix(0, 6*6, 6)
  linv <- matrix(0, 6*6, 6)
  pi <- rep(0, 6)
  for(i in 4001:8000){
    for(k in 1:6){
      means[(k-1)*6 + 1:6] <- means[(k-1)*6 + 1:6] + mixDraws$draws[[1]][[k]]$mean[i,] / 4000
      uinv <- chol(mixDraws$draws[[1]][[k]]$varInv[,,i])
      linv[(k-1)*6 + 1:6,] <-  linv[(k-1)*6 + 1:6,] + uinv / 4000
      varinv[(k-1)*6 + 1:6,] <-  varinv[(k-1)*6 + 1:6,] + t(uinv) %*% uinv / 4000
      for(n in 1:2000){
        pi[k] <- pi[k] + mixDraws$draws[[1+n]]$pi[i, k] / (2000 * 4000)
      }
    }
  }
  dets <- numeric(6)
  for(k in 1:6){
    dets[k] <- det(linv[(k-1)*6 + 1:6, ])  
  }
  priorMix <- list(mean = means, linv = linv, varInv = varinv, dets = dets, weights = pi)
  startingLam <- c(means, rep(c(diag(0.5, 6)), 6), rep(1, 6))

}

otherMCMCModels{
  data <- readRDS('MCMCData.RDS')
  reps <- 50000
  
  draws <- list()
  hyper <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(10, 6)), v = 6, scale = diag(1, 6))
  draws[[1]] <- list(mean = c(-5, -5, 0, 0, 0, 0), varInv = diag(10, 6))
  for(i in 1:N){
    draws[[i+1]] <- c(-5, -5, 0, -0.1, 0.15, 0.05)
  }

  noMixDraws <- hierNoMixMCMC(data, reps, draws, hyper, thin = 10)
  saveRDS(noMixDraws, 'noMixN2000.RDS')
  noMixDraws <- readRDS('noMixN2000.RDS')
  
  noMixDraws$draws[[1]]$mean %>%
    as.data.frame() %>%
    cbind(iter = 1:(reps/thin)) %>%
    mutate(V1 = exp(V1), V2 = exp(V2)) %>%
    rename(sigSq_eps = V1, sigSq_eta = V2, phi1 = V3, phi2 = V4, gamma1 = V5, gamma2 = V6) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
    filter(iter > 1000) %>%
    ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free')
  
  noMixMean <- rep(0, 6)
  noMixVar <- matrix(0, 6, 6)
  
  support <- data.frame(seq(exp(-13), exp(-4.5), length.out = 1000),
                        seq(exp(-15), exp(-7.6), length.out = 1000),
                        seq(-0.5, 2, length.out = 1000),
                        seq(-1, 0.8, length.out = 1000),
                        seq(0, 2, length.out = 1000),
                        seq(-1, 0.5, length.out = 1000))
  
  densities <- matrix(0, 1000, 6)
  for(i in (reps/(2*thin)+1):(reps/thin)){
    noMixMean <- noMixMean + noMixDraws$draws[[1]]$mean[i,] / (reps / (2 * thin))
    noMixVar <- noMixVar + solve(noMixDraws$draws[[1]]$varInv[,,i]) / (reps / (2 * thin))
    for(k in 1:2){
      densities[,k] <- densities[,k] + dlnorm(support[,k], noMixDraws$draws[[1]]$mean[i, k], sqrt(solve(noMixDraws$draws[[1]]$varInv[,,i])[k,k])) / (reps / (2 * thin))
    }
    for(k in 3:6){
      densities[,k] <- densities[,k] + dnorm(support[,k], noMixDraws$draws[[1]]$mean[i, k], sqrt(solve(noMixDraws$draws[[1]]$varInv[,,i])[k,k])) / (reps / (2 * thin))
    }
  }
  noMixLam <- c(noMixMean, chol(noMixVar))
  densities <- as.data.frame(densities)
  colnames(densities) <- c('sigSq_eps', 'sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2')
  densities %>% 
    gather(var, dens) %>%
    cbind(support = unlist(c(support))) %>%
    mutate(var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
    ggplot() + geom_line(aes(support, dens)) + facet_wrap(~var, scales = 'free')
    
  
  draws <- c(-5, -5, 0, 0, 0, 0)
  hyper <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(10, 6)))
  hyper$var <- diag(solve(hyper$varInv))
  noHierDraws <- noHierMCMC(data, reps, draws, hyper, thin = 10, stepsize = 0.01)
  saveRDS(noHierDraws, 'noHierN2000.RDS')
  noHierDraws <- readRDS('noHierN2000.RDS')
  
  noHierDraws$draws %>%
    as.data.frame() %>%
    cbind(iter = 1:(reps/thin)) %>%
    mutate(V1 = exp(V1), V2 = exp(V2)) %>%
    rename(sigSq_eps = V1, sigSq_eta = V2, phi1 = V3, phi2 = V4, gamma1 = V5, gamma2 = V6) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
    filter(iter > 100) %>%
    ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free')
  
  noHierMean <- rep(0, 6)
  noHierVar <- var(noHierDraws$draws[(reps/(2*thin)+1):(reps/thin),])
  for(i in (reps/(2*thin)+1):(reps/thin)){
    noHierMean <- noHierMean + noHierDraws$draws[i,] / (reps / (2 * thin))
  }
  noHierLam <- c(noHierMean, chol(noHierVar))
  
}

VBSingleCar{
  
lags <- 2
dim <- 2 + 2 * lags
stepsize <- 10

n <- stepsize + lags
idUpdate <- id$fullID[!id$fullID %in% id$idSubset]
idU <- sample(idUpdate, 1)

carsAug %>%
  filter(ID == idU) %>%
  mutate(n = seq_along(v),
         vl = ifelse(n == 1, 0, lag(v)),
         a = v - lag(v),
         d = delta - pi/2) %>%
  filter(n > 1 & n <= 501) %>% 
  select(a , d) %>%
  as.matrix() -> dataUpdate

starting <- seq(1, nrow(dataUpdate), stepsize)
starting <- starting[-length(starting)]

# Fit distribution to MCMC - MVN Mixutre
mean1 <- rep(0, 6)
cov1 <- matrix(0, 6, 6)
mean2 <- rep(0, 6)
cov2 <- matrix(0, 6, 6)
w <- 0
for(i in (reps/(2*thin)+1):(reps/thin)){
  mean1 <- mean1 + mixDraws$draws[[1]][[1]]$mean[i,] / (reps / (2*thin))
  cov1 <- cov1 + solve(mixDraws$draws[[1]][[1]]$varInv[,,i]) / (reps / (2*thin))
  mean2 <- mean2 + mixDraws$draws[[1]][[2]]$mean[i,] / (reps / (2*thin))
  cov2 <- cov2 + solve(mixDraws$draws[[1]][[2]]$varInv[,,i]) / (reps / (2*thin))
  for(j in 1:N){
    w <- w + mixDraws$draws[[j+1]]$pi[i, 1] / (N * reps / (2 * thin))
  }
}
uinv1 <- solve(chol(cov1))
uinv2 <- solve(chol(cov2))
lambda <- matrix(c(mean1, mean2, uinv1, uinv2, log(w/(1-w))), ncol = 1)

### Check if it should be Linv or Uinv ###

for(t in seq_along(starting)){
  mean <- lambda[1:(2*dim)]
  u1 <- matrix(lambda[(2*dim+1):(2*dim + dim^2)], dim)
  linv1 <- solve(t(u1))
  u2 <- matrix(lambda[(dim*(2+dim)+1):(length(lambda)-1)], dim)
  linv2 <- solve(t(u2))
  linv <- rbind(linv1, linv2)
  w <- 1 / (1 + exp(-tail(lambda, 1)))
  det1 <- abs((2*pi)^(-3) * prod(diag(linv1)))
  det2 <- abs((2*pi)^(-3) * prod(diag(linv2)))
  priorComp <- c(w, det1, det2)
  dat <- dataUpdate[starting[t]:min(nrow(data), starting[t]+n-1),]
  if(is.matrix(dat)){
    if(nrow(dat) > lags){
      lambda <- carsVB(dat, lambda, dimTheta = 2*dim, model = arUpdaterMix,
                       mean = mean, Linv = linv, lags = lags, priorComp = priorComp)$lambda
    }
  } 
}


support <- data.frame(seq(-7, -1, length.out = 1000),
                      seq(-11, -6, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000),
                      seq(-2, 2, length.out = 1000))
mean1 <- lambda[1:6]
mean2 <- lambda[7:12]
u1 <- matrix(lambda[13:48], 6)
u2 <- matrix(lambda[49:84], 6)
sd1 <- abs(diag(u1))
sd2 <- abs(diag(u2))
w <- 1 / (1 + exp(-lambda[85]))

denVB <- data.frame()
for(i in 1:6){
  den1 <- dnorm(support[, i], mean1[i], sd1[i])
  den2 <- dnorm(support[, i], mean2[i], sd2[i])
  denVB <- rbind(denVB, data.frame(
    support = support[,i],
    dens = w * den1 + (1-w) * den2,
    var = vars[i]))
}
ggplot(denVB) + geom_line(aes(support, dens)) + facet_wrap(~var, scales = 'free')
}

IndepMCMC{
  
N <- 2000
posMeans <- NULL
vars <- c('sigSq_eps', 'sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2') 
for(i in 1:N){
  carsAug %>%
    filter(ID == id$idSubset[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           a = v - lag(v),
           d = delta - pi/2) %>%
    filter(n > 1 & n <= 501) %>% 
    select(a , d) %>%
    as.matrix() -> data
  
  draw <- c(-5, -5, 0, 0, 0, 0)
  hyper <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(10, 6)))
  
  MCMC <- singleMCMCallMH(data, 10000, draw, hyper)$draws
  mean <- c(colMeans(exp(MCMC[5001:10000,1:2])), colMeans(MCMC[5001:10000, 3:6]))
  
  posMeans <- rbind(posMeans, data.frame(mean = mean, var = vars, ID = id$idSubset[i]))
  if(i %% 50 == 0){
    print(i)
  }
}

# Fill in densMixMod and p1 from the Hierarchical Model results above (crtl-f 'densMixMod')

posMeans %>%
  mutate(var = as.character(var),
         var = ifelse(var == 'log_sigSq_eps', 'sigSq_eps', 
                      ifelse(var == 'log_sigSq_eta', 'sigSq_eta', var)),
         var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
  filter((var == 'sigSq_eps' & mean > exp(-13) & mean < exp(-4.5)) |
         (var == 'sigSq_eta' & mean > exp(-15) & mean < exp(-7.6)) | 
         (var == 'phi1' & mean > -0.5 & mean < 2) |
         (var == 'phi2' & mean > -1 & mean < 0.8) | 
         (var == 'gamma1' & mean > 0 & mean < 2) | 
         (var == 'gamma2' & mean > -1 & mean < 0.5)) %>%
  ggplot() + geom_density(aes(mean)) + 
  geom_blank(data = densMixMod, aes(x = support)) + 
  facet_wrap(~var, scales = 'free', ncol = 6) +
  labs(title = 'Single Model Means (Kernel Density Estimate)', x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))-> p2

gridExtra::grid.arrange(p2, p1, ncol = 1)
}

sliceSampler{
  set.seed(1)
  N <- 100
  idSubset <- sample(noStop, N)
  data <- list()
  for(i in 1:N){
    carsAug %>%
      filter(ID == idSubset[i]) %>%
      mutate(n = seq_along(v),
             vl = ifelse(n == 1, 0, lag(v)),
             a = v - lag(v),
             d = delta - pi/2) %>%
      filter(n > 1 & n <= 501) %>% 
      select(a , d) %>%
      as.matrix() -> data[[i]]
  }
  reps <- 100000
  maxK <- 20
  thin <- 10
  burn <- 0.9
  
  draws <- list(list())
  hyper <- list()
  for(k in 1:maxK){
    hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(c(5, 5, 10, 10, 10, 10))), v = 6, scale = diag(1, 6))
    draws[[1]][[k]] <- list(mean = c(-5, -5, 0, 0, 0, 0), varInv = diag(10, 6))
  }
  draws[[1]]$pi <- rep(1/maxK, maxK)
  hyper$M <- 1
  for(i in 1:N){
    draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:5, 1))
  }
  sliceDraws <- sliceSampler(data, reps, draws, hyper, thin, 10, 'gaussian', 0.01)
}

forecastResults {

results <- NULL
for(i in seq_along(id$idfc)){
  try(assign('results',
             rbind(results,
                   read.csv(paste0('eval/car', id$idfc[i], '.csv')))))
  if((i %% 100 == 0) | i == length(id$idfc)){
    print(i)
  }
}

mPerFoot <- 0.3048
results %>% 
  mutate(logscore = logscore + log(mPerFoot^2),
         mapDist = mapDist * mPerFoot,
         consDist = consDist * mPerFoot) -> results

results %>%
  group_by(method, prior, h) %>%
  summarise(meanDist = mean(mapDist)) %>%
  mutate(model = paste(method, prior)) %>%
  ggplot() + geom_line(aes(h, meanDist, colour = model)) + 
  geom_line(data = results %>% group_by(h) %>% summarise(meanDist = mean(consDist)), aes(h, meanDist)) + 
  labs(x = 'Forecast Horizon (seconds)', y = 'Mean Euclidean Error (metres)') + 
  theme_bw() + 
  scale_x_continuous(labels = c(1, 2, 3), breaks = c(10, 20, 30))

results %>% 
  mutate(group = factor(ceiling(S / 30)),
         method = ifelse(method == 'VB-Stream', 'VB-Update', as.character(method))) %>% 
  ggplot() + geom_boxplot(aes(x = group, y = logscore, colour = method)) + 
  facet_wrap(~prior, ncol = 1) + 
  ylim(-10, 5) +
  labs(x = 'T', y = 'Predictive Logscore') + 
  scale_x_discrete(labels = c('10-30', '40-60', '70-90', '100-120', '130-150', 
                              '160-180', '190-210', '220-240', '250-270', '280-300'))  

results %>%
  filter(h %in% c(10, 20, 30) & method == 'VB-Stream' & prior == 'Finite Mixture') %>% 
  mutate(horizon = paste(h / 10, ifelse(h == 10, 'second', 'seconds'), 'ahead')) %>%  
  group_by(S, horizon) %>%
  summarise(`Mixture/VB Update Model Predictive MAP` = mean(mapDist),
            `Constant Velocity and Angle` = mean(consDist)) %>%
  ungroup() %>%
  gather(Predictor, meanDist, -horizon, -S) %>%
  ggplot() + geom_line(aes(x = S, y = meanDist, colour = Predictor)) + 
  facet_wrap(~horizon) + 
  labs(x = 'T', y = 'Mean Euclidean Error (metres)') + 
  theme_bw() +  
  theme(legend.position = 'bottom') 

results %>%
  filter(h == 10) %>%
  group_by(S, method, prior) %>%
  summarise(meanDist = mean(mapDist)) %>%
  ggplot() + geom_line(aes(S, meanDist, colour = method)) + facet_wrap(~prior)

results %>% 
  group_by(method, prior) %>% 
  filter(is.finite(logscore)) %>% 
  summarise(med = median(logscore, na.rm = TRUE)) %>% 
  arrange(med)

results %>% 
  filter(h == 30 & method == 'VB-Stream' & prior == 'Finite Mixture') %>% 
  summarise(meanMAP = mean(mapDist), meanConst = mean(consDist)) 

results %>%
  group_by(method, prior) %>%
  summarise(meanError = mean(mapDist))  %>%
  arrange(meanError)
  
mean(results$consDist)
 
results[!is.na(results$xCDF),] %>%
  ggplot() + geom_density(aes(xCDF, colour = prior)) + facet_wrap(~method, ncol = 1)
 
}

noGapsSampler {
  set.seed(1)
  N <- 500
  idSubset <- sample(id$idSubset, N)
  data <- list()
  for(i in 1:N){
    carsAug %>%
      filter(ID == idSubset[i]) %>%
      mutate(n = seq_along(v),
             vl = ifelse(n == 1, 0, lag(v)),
             a = v - lag(v),
             d = delta - pi/2) %>%
      filter(n > 1 & n <= 501) %>% 
      select(a , d) %>%
      as.matrix() -> data[[i]]
  }
  reps <- 25000
  thin <- 10
  startK <- 15
  
  draws <- list(list())
  hyper <- list(mean = c(-7, -7, rep(0, 4)),
                varInv = solve(diag(c(5, 5, 1, 1, 1, 1))),
                df = 6, 
                scale = diag(1, 6),
                alpha = 1)
  for(k in 1:startK){
    draws[[1]][[k]] <- list(mean = c(-5, -5, 0, 0, 0, 0),
                            varInv = diag(10, 6))
  }
  for(i in 1:N){
    draws[[i+1]] <- list(theta = c(mvtnorm::rmvnorm(1, hyper$mean, solve(hyper$varInv))))
  }
  noGapDraws <- NoGaps(data, reps, draws, hyper, thin, startK, 0.01)
  data <- readRDS('MCMCData.RDS')
  
  reps <- 75000
  thin <- 10
  startK <- 10
  
  hyper <- list(mean = c(-5, -5, 0, 0, 0, 0),
                varInv = diag(c(0.1, 0.1, 1, 1, 1, 1)),
                var = diag(c(10, 10, 1, 1, 1, 1)),
                alpha = 1)
  
  draws <- list()
  for(i in 1:startK){
    draws[[i]] <- rnorm(6, hyper$mean, sqrt(diag(hyper$var)))
  }
  
  noGapDraws <- NoGaps2(data, reps, draws, hyper, thin, startK, 0.01)
  
}

neuralNetwork {
  
library(tensorflow)
library(keras)
xset <- data.frame()
yset <- data.frame()
for(i in 1:20){
  carsAug %>%
    filter(ID == id$idSubset[i]) %>%
    mutate(n = seq_along(v),
           d = delta - pi/2) %>%
    filter(n > 1 & n <= 501) %>% 
    ungroup() %>%
    select(x, y, v, d) -> carI
  T <- nrow(carI)
  for(t in 11:(T - 30)){
    carsub <- carI[(t-10):(t+30),]
    xset <- rbind(xset, unlist(carsub[1:10,]))  
    yset <- rbind(yset, unlist(carsub[11:40, 1:2]))

  }
  if(i %% 5 == 0){
    print(i)
  }
}

K <- backend()
euclideanLoss <- function(y_true, y_pred){
  K$mean(K$sqrt(K$square(y_true[1:30] - y_pred[1:30]) + K$square(y_true[31:60] - y_pred[31:60])) + 1e-8)
}

train <- sample(1:nrow(xset), 0.8 * nrow(xset))
x_train <- xset[train,] %>% as.matrix()
x_test <- xset[-train,] %>% as.matrix()
y_train <- yset[train,] %>% as.matrix()
y_test <- yset[-train,] %>% as.matrix()


model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 64, input_shape = c(40)) %>% 
  layer_activation('relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = 64) %>% 
  layer_activation('relu') %>%
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = 64) %>% 
  layer_activation('relu') %>%
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = 60) %>%
  layer_activation('linear')

model %>% compile(
  optimizer = 'adam',
  loss = euclideanLoss
)

model %>% fit(x_train, y_train, epochs=100, batch_size=128)
model %>% predict(x_test) -> y_pred

error <- numeric(0)
for(i in 1:30){
  for(j in 1:nrow(y_test)){
    error <- c(error, sqrt((y_pred[j, i] - y_test[j, i])^2 + (y_pred[j, i+30] - y_test[j, i+30])^2))
  }
}
mean(error)

}


