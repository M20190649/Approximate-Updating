########## TO DO #############
# Write more things.
# Build the Neural Network (find why loss goes to nan)

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
  hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(c(10, 10, 10, 10, 10, 10))), v = 6, scale = diag(1, 6))
  draws[[1]][[k]] <- list(mean = c(-5, -5, 0, 0, 0, 0), varInv = diag(10, 6))
}
draw[[1]]$pi <- rep(1/K, K)
hyper$alpha <- rep(1, K)
for(i in 1:N){
  draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:K, 1), pi = rep(1/K, K))
}



mixDraws <- mixtureMCMC(data, reps, draws, hyper, thin, K, 'gaussian', 0.01)

saveRDS(mixDraws, 'mixN2000K6.RDS')
mixDraws <- readRDS(file = 'mixN2000K6.RDS')

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
  mutate(var = case_when(
    var == 'sigSq_eps' ~ 'sigma[epsilon]^{2}',
    var == 'phi1' ~ 'phi[1]',
    var == 'phi2' ~ 'phi[2]',
    var == 'sigSq_eta' ~ 'sigma[eta]^{2}',
    var == 'gamma1' ~ 'gamma[1]',
    var == 'gamma2' ~ 'gamma[2]',
    TRUE ~ as.character(var)
  ),
  var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                      'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) -> densities
densities %>% 
  group_by(var, support) %>%
  summarise(dens = sum(dens)) -> densMixMod

  ggplot(densMixMod) + geom_line(aes(support, dens)) +
  geom_line(data=densities, aes(support, dens, colour = factor(group))) +
  facet_wrap(~var, scales = 'free', ncol = 6, labeller = label_parsed) + 
  labs(x = NULL, y = NULL) + 
  theme_bw() + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))  -> p1
  

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
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `phi[1]` = V3, `phi[2]` = V4, `gamma[1]` = V5, `gamma[2]` = V6) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                 'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) %>% 
    filter(iter > 100) %>%
    ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free', labeller = label_parsed)
  
  noHierDraws$draws %>%
    as.data.frame() %>%
    cbind(iter = 1:(reps/thin)) %>%
    mutate(V1 = exp(V1), V2 = exp(V2)) %>%
    rename(`sigma[epsilon]^{2}` = V1, `sigma[eta]^{2}` = V2, `phi[1]` = V3, `phi[2]` = V4, `gamma[1]` = V5, `gamma[2]` = V6) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                        'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) %>% 
    filter(iter > 1000) %>%
    ggplot() + geom_density(aes(draw)) + facet_wrap(~var, scales = 'free', ncol = 6, labeller = label_parsed) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
    labs(x = NULL, y = NULL)
  
  
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
vars <- c('sigma[epsilon]^{2}', 'sigma[eta]^{2}', 'phi[1]', 'phi[2]', 'gamma[1]', 'gamma[2]') 
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
saveRDS(posMeans, 'posMeans.RDS')

# Fill in densMixMod and p1 from the Hierarchical Model results above (crtl-f 'densMixMod')

posMeans %>%
  mutate(var = as.character(var),
         #var = ifelse(var == 'log_sigSq_eps', 'sigSq_eps', 
        #              ifelse(var == 'log_sigSq_eta', 'sigSq_eta', var)),
         var = factor(var, levels = c('sigma[epsilon]^{2}', 'phi[1]', 'phi[2]',
                                       'sigma[eta]^{2}', 'gamma[1]', 'gamma[2]'))) %>%
  filter((var == 'sigma[epsilon]^{2}' & mean > exp(-13) & mean < exp(-4.5)) |
         (var == 'sigma[eta]^{2}' & mean > exp(-15) & mean < exp(-7.6)) | 
         (var == 'phi[1]' & mean > -0.5 & mean < 2) |
         (var == 'phi[2]' & mean > -1 & mean < 0.8) | 
         (var == 'gamma[1]' & mean > 0 & mean < 2) | 
         (var == 'gamma[2]' & mean > -1 & mean < 0.5)) %>%
  ggplot() + geom_density(aes(mean)) + 
  geom_blank(data = densMixMod, aes(x = support)) + 
  facet_wrap(~var, scales = 'free', ncol = 6, labeller = label_parsed) +
  labs(x = NULL, y = NULL) +
  theme_bw() + 
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

# Read in raw results, transform from ft to metres
# Missing car 20771, car 20783, car 20813,
library(readr)
results <- NULL
for(i in seq_along(id$idfc)){
  if(i %in% c(687, 690, 699)){
    next
  }
  results <- rbind(results,
                   read_csv(paste0('forecasts/car', id$idfc[i], '.csv'), col_types = cols()))
  if((i %% 100 == 0) | i == length(id$idfc)){
    print(i)
  }
}
for(i in 1:500){
  results <- rbind(results,
                   read_csv(paste0('forecasts/car', noStopChanged[i], '.csv'), col_types = cols()))
  if(i %% 100 == 0){
    print(i)
  }
}

mPerFoot <- 0.3048
results %>% 
  mutate(logscore = logscore + log(mPerFoot^2),
         dist = dist * mPerFoot) -> results
write.csv(results, 'results.csv', row.names = FALSE)
results <- readr::read_csv('results.csv')

# Average prediction Error

results[results$h %in% c(10, 20, 30) &
          !results$model %in% c('Single Hierarchy MCMC', 'Single Hierarchy VB-Standard', 'Single Hierarchy VB-Updating'), ] %>%
  mutate(horizon = paste(h / 10, ifelse(h == 10, 'second', 'seconds'), 'ahead'),
         model = sub('Non-Informative', 'IH', model),
         model = sub('Finite Mixture', 'CH', model),
         model = sub('VB-Standard', 'VB', model),
         model = sub('VB-Updating', 'UVB', model),
         model = factor(model, levels = c('Naive 1', 'Naive 2', 'Naive 3', 'Naive 4', 'Naive 5', 
                                          'Naive 6', 'Naive 7', 'Naive 8', 'Naive 9', 'Homogenous MCMC', 
                                          'IH MCMC', 'IH VB', 'IH UVB',
                                          'CH MCMC', 'CH VB', 'CH UVB'))) %>%  
  group_by(S, horizon, model) %>%
  group_by(model, h) %>%
  summarise(mean = mean(dist)) %>% 
  spread(model, mean)

# Plot updating vs standard logscores

results[!is.na(results$logscore) & 
        results$model %in% c('Finite Mixture VB-Standard', 'Finite Mixture VB-Updating',
                             'Non-Informative VB-Standard', 'Non-Informative VB-Updating') & 
         results$h == 30,] %>%
  select(model, logscore, S, id) %>%
  mutate(model = ifelse(model == 'Non-Informative VB-Standard', 'Independent VB-Standard', 
                        ifelse(model == 'Non-Informative VB-Updating', 'Independent VB-Updating', model))) %>% 
  spread(model, logscore) %>%
  mutate(Clustered = `Finite Mixture VB-Updating` - `Finite Mixture VB-Standard`,
         Independent = `Independent VB-Updating` - `Independent VB-Standard`,
         S = ceiling((S - 100) / 50)) %>%
  filter(S > 0) %>%
  select(S, id, Clustered, Independent) %>%
  gather(Model, Difference, -S, -id) %>%
  ggplot() + geom_boxplot(aes(factor(S), Difference)) + 
  facet_wrap(~Model, ncol = 2) + 
  labs(x = 'T (100 milliseconds)', y = 'Difference In Logscore (UVB - VB)') + 
  theme_bw() +
  ylim(-0.5, 0.5) + 
  scale_x_discrete(labels = c('110-150', '160-200', '210-250', '260-300', '310-350', '360-400', '410-450'))

# Plot point estimate errors by model

results[results$h %in% c(10, 20, 30) &
        !results$model %in% c('Single Hierarchy MCMC', 'Single Hierarchy VB-Standard', 'Single Hierarchy VB-Updating'), ] %>%
  mutate(horizon = paste(h / 10, ifelse(h == 10, 'second', 'seconds'), 'ahead'),
         model = sub('Non-Informative', 'IH', model),
         model = sub('Finite Mixture', 'CH', model),
         model = sub('VB-Standard', 'VB', model),
         model = sub('VB-Updating', 'UVB', model),
         model = factor(model, levels = c('Naive 1', 'Naive 2', 'Naive 3', 'Naive 4', 'Naive 5', 
                                          'Naive 6', 'Naive 7', 'Naive 8', 'Naive 9', 'Homogenous MCMC', 
                                          'IH MCMC', 'IH VB', 'IH UVB',
                                          'CH MCMC', 'CH VB', 'CH UVB'))) %>%  
  rename(Model = model) %>%
  group_by(S, horizon, Model) %>%
  summarise(meanError = mean(dist)) %>%
  ungroup() %>%
  ggplot() + geom_line(aes(x = S, y = meanError, colour = Model)) + 
  facet_wrap(~horizon) + 
  labs(x = 'T', y = 'Mean Euclidean Error (metres)') + 
  theme_bw() +  
  theme(legend.position = 'bottom') 

# As above but exclude the naive models

results[results$h %in% c(10, 20, 30) &
        results$model %in% c('Homogenous MCMC', 'Non-Informative MCMC', 'Non-Informative VB-Standard', 'Non-Informative VB-Updating',
                             'Finite Mixture MCMC', 'Finite Mixture VB-Standard', 'Finite Mixture VB-Updating'), ] %>%
  mutate(horizon = paste(h / 10, ifelse(h == 10, 'second', 'seconds'), 'ahead'),
         model = sub('Non-Informative', 'IH', model),
         model = sub('Finite Mixture', 'CH', model),
         model = sub('VB-Standard', 'VB', model),
         model = sub('VB-Updating', 'UVB', model),
         model = factor(model, levels = c('Homogenous MCMC', 
                                          'IH MCMC', 'IH VB', 'IH UVB',
                                          'CH MCMC', 'CH VB', 'CH UVB'))) %>%  
  rename(Model = model) %>%
  group_by(S, horizon, Model) %>%
  summarise(meanError = mean(dist)) %>%
  ungroup() %>%
  ggplot() + geom_line(aes(x = S, y = meanError, colour = Model)) + 
  facet_wrap(~horizon) + 
  labs(x = 'T', y = 'Mean Euclidean Error (metres)') + 
  theme_bw() +  
  theme(legend.position = 'bottom') 

# Read in individual variances, plot against homogenous variances, and then split results by variance quintiles

sigma <- data.frame()
for(i in seq_along(id$idfc)){
  sigma <- rbind(sigma,
                 readr::read_csv(paste0('sigma/car', id$idfc[i], '.csv'), col_types = cols()))
}
for(i in 1:500){
  sigma <- rbind(sigma,
                 readr::read_csv(paste0('sigma/car', noStopChanged[i], '.csv'), col_types = cols()))
}
write.csv(sigma, 'sigma.csv')
sigma <- readr::read_csv('sigma.csv')

homogDraws <- readRDS('noHierN2000.RDS')$draws
homogMean = colMeans(exp(homogDraws[1001:5000, 1:2]))
homogMeanDf <- data.frame(variable = c('Acceleration', 'Angle'), 
                        posMean = colMeans(exp(homogDraws[1001:5000, 1:2])))


sigma %>%
  filter(method == 'MCMC' & prior == 'Finite Mixture') %>%
  mutate(prior = 'Hierarchy',
         variable = ifelse(variable == 'a', 'Acceleration', 'Angle')) %>%
  group_by(variable) %>%
  filter(posMean < quantile(posMean, 0.99)) %>%
  ggplot() + geom_density(aes(posMean)) +
  geom_vline(data = homogMeanDf, aes(xintercept = posMean), colour = 'red') +
  facet_wrap(~variable, scales = 'free', ncol = 2) + 
  theme_bw() + 
  labs(x = 'Variance Posterior Mean', y = NULL)

sigma %>%
  filter(method == 'MCMC' & prior == 'Finite Mixture' & variable == 'a') %>%
  .$posMean %>%
  quantile(seq(0.2, 1, 0.2)) -> aQuant

sigma %>%
  filter(method == 'MCMC' & prior == 'Finite Mixture' & variable == 'd') %>%
  .$posMean %>%
  quantile(seq(0.2, 1, 0.2)) -> dQuant


sigma %>%
  select(-X1) %>%
  filter(method == 'MCMC' & prior == 'Finite Mixture') %>%
  spread(variable, posMean) %>%
  mutate(aVar = ifelse(a < aQuant[1], 1, ifelse(a < aQuant[2], 2, ifelse(a < aQuant[3], 3, ifelse(a < aQuant[4], 4, 5)))),
         dVar = ifelse(d < dQuant[1], 1, ifelse(d < dQuant[2], 2, ifelse(d < dQuant[3], 3, ifelse(d < dQuant[4], 4, 5)))),
         aVar = factor(aVar),
         dVar = factor(dVar)) %>%
  select(id, aVar, dVar) -> varianceCategoryMix

results[!is.na(results$logscore) &
          results$h == 30 & 
          results$id %in% unique(varianceCategoryMix$id) & 
          results$model %in% c('Finite Mixture VB-Updating', 'Homogenous MCMC') & 
          results$S <= 400, ] %>%
  left_join(varianceCategoryMix) %>%
  group_by(model, S, h, aVar, dVar) %>%
  summarise(ls = median(logscore)) %>%
  ungroup() %>%
  spread(model, ls) %>%
  mutate(diff = `Finite Mixture VB-Updating` - `Homogenous MCMC`, 
         Model = 'Hierarchy') %>%
  select(S, aVar, dVar, diff, Model) -> hierarchyVarSplit

hierarchyVarSplit %>%
  ggplot() + geom_line(aes(S, diff)) +
  geom_hline(aes(yintercept = 0), colour = 'red') + facet_grid(aVar ~ dVar) + 
  labs(y = 'Median Difference in h = 30 Logscore (CH - Homogenous)', x = 'T (100 milliseconds)') + 
  theme_bw() + 
  theme(legend.position = 'bottom')

sigma %>%
  filter(method == 'MCMC' & prior == 'Non-Informative' & variable == 'a') %>%
  .$posMean %>%
  quantile(seq(0.2, 1, 0.2)) -> aQuantInd

sigma %>%
  filter(method == 'MCMC' & prior == 'Non-Informative' & variable == 'd') %>%
  .$posMean %>%
  quantile(seq(0.2, 1, 0.2)) -> dQuantInd


sigma %>%
  filter(method == 'MCMC' & prior == 'Non-Informative') %>%
  spread(variable, posMean) %>%
  mutate(aVar = ifelse(a < aQuant[1], 1, ifelse(a < aQuant[2], 2, ifelse(a < aQuant[3], 3, ifelse(a < aQuant[4], 4, 5)))),
         dVar = ifelse(d < dQuant[1], 1, ifelse(d < dQuant[2], 2, ifelse(d < dQuant[3], 3, ifelse(d < dQuant[4], 4, 5)))),
         aVar = factor(aVar),
         dVar = factor(dVar)) %>%
  select(id, aVar, dVar) -> varianceCategoryInd

results[!is.na(results$logscore) &
          results$h == 30 & 
          results$id %in% unique(varianceCategoryInd$id) & 
          results$model %in% c('Non-Informative VB-Updating', 'Homogenous MCMC') & 
          results$S <= 400, ] %>%
  left_join(varianceCategoryInd) %>%
  group_by(model, S, h, aVar, dVar) %>%
  summarise(ls = median(logscore)) %>%
  ungroup() %>%
  spread(model, ls) %>%
  mutate(diff = `Non-Informative VB-Updating` - `Homogenous MCMC`, 
         Model = 'Independent') %>%
  select(S, aVar, dVar, diff, Model) -> independentVarSplit

independentVarSplit %>%
  rbind(hierarchyVarSplit) %>%
  ggplot() + geom_line(aes(S, diff, colour = Model)) + 
  geom_hline(aes(yintercept = 0), colour = 'black') + facet_grid(aVar ~ dVar) + 
  labs(y = 'Median Difference in h = 30 Logscore (Heterogenous - Homogenous)', x = 'T (100 milliseconds)') + 
  theme_bw() + 
  theme(legend.position = 'bottom')
 

# Read in timings
timing <- data.frame()
for(i in 1:500){
  tmp <- readr::read_csv(paste0('timing/car', id$idfc[i], '.csv'), col_types = cols())
  timing <- rbind(timing, tmp)
}

timing %>%
  filter(prior == 'Hierarchy') %>% 
  group_by(method, T, TmS) %>%
  summarise(meanTime = mean(time)) %>%
  ungroup() %>% 
  mutate(`T - S` = factor(TmS)) -> meanTimings

  ggplot() + geom_line(data = filter(meanTimings, method == 'Standard'), aes(T, meanTime)) + 
    geom_line(data = filter(meanTimings, method == 'Updating'), aes(T, meanTime, colour = `T - S`)) + 
    labs(x = 'T (100 milliseconds)', y = 'Mean Time to Converge (seconds)') + 
    theme_bw() + 
    ylim(0, 25)


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

recurrentNeuralNetwork {
  
library(tensorflow)
library(keras)
sourceCpp('buildNNData.cpp')


data <- list()
for(i in 1:2000){
  carsAug[carsAug$ID == id$idSubset[i],] -> carI
  carI[2:min(501, nrow(carI)),] %>%
    mutate(d = delta - pi/2) %>%
    select(x, y, v, d) -> carI
  data[[i]] <- carI
}

K <- backend()
euclideanLoss <- function(y_true, y_pred){
  mPerFoot * K$mean(
                    K$sqrt(K$square(K$sum(y_true[,,1]) - K$sum(y_pred[,,1])) +
                           K$square(K$sum(y_true[,,2]) - K$sum(y_pred[,,2]))
                           + 1e-8
                    )
             )
}

eucLossVD <- function(y_true, y_pred){
  mPerFoot * K$mean(
                    K$sqrt(K$square(K$sum(exp(y_true[,1:10,1]) * cos(y_true[,1:10,2] + pi/2)) - 
                                    K$sum(exp(y_pred[,1:10,1]) * cos(y_pred[,1:10,2] + pi/2))) + 
                           K$square(K$sum(exp(y_true[,1:10,1]) * sin(y_true[,1:10,2] + pi/2)) - 
                                    K$sum(exp(y_pred[,1:10,1]) * sin(y_pred[,1:10,2] + pi/2))) + 1e-8
                           )
                    )
}

generator <- function(data, numCars, batch_size = 128){
  i <- 1
  function(){
    samples <- array(0, dim = c(batch_size, 10, 4))
    targets <- array(0, dim = c(batch_size, 10, 2))
    dataSub <- data[[i]]
    if(i >= numCars){
      i <<- 1
    } else {
      i <<- i + 1
    }
    t <- sample(30:(nrow(dataSub) - 5), batch_size)
    for(j in 1:batch_size){
      for(k in 1:4){
        samples[j, , k] <- dataSub[(-9:0) + t[j], k]
      }
      samples[j,,3] <- log(samples[j,,3])
      targets[j,,1] <- log(dataSub[(1:10) + t[j], 3])
      targets[j,,2] <- dataSub[(1:10) + t[j], 4]
    }
    
    list(samples, targets)
  }
}
trainSet <- sample(1:2000, 1600)
train_gen <- generator(data[trainSet], numCars = 1600)
test_gen <- generator(data[-trainSet], numCars = 400)



rnn <- keras_model_sequential()
rnn %>%
  layer_gru(units = 128, 
            dropout = 0.1, 
            recurrent_dropout = 0.5,
            return_sequences = FALSE,
            input_shape = list(NULL, 4)) %>% 
  layer_repeat_vector(10) %>%
  layer_gru(units = 128,
            dropout = 0.1,
            recurrent_dropout = 0.5,
            return_sequences = TRUE) %>%
  time_distributed(layer_dense(units = 128, 
                               activation = "relu")) %>%
  layer_dropout(rate = 0.5) %>%
  time_distributed(layer_dense(units = 2))

rnn %>% compile(
  optimizer = 'adam',
  loss = eucLossVD
)

history <- rnn %>% fit_generator(
  train_gen,
  steps_per_epoch = 50,
  epochs = 25,
  validation_data = test_gen,
  validation_steps = 25
)
plot(history)

nn <- keras_model_sequential() 
nn %>% 
  layer_dense(256,
              input_shape = 40,
              activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(256,
              activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(256,
              activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(20)

nnLoss <- function(y_true, y_pred) {
  mPerFoot * K$mean(K$sqrt(K$square(K$sum(y_true[1:10]) - K$sum(y_pred[1:10])) +
                           K$square(K$sum(y_true[11:20]) - K$sum(y_pred[11:20])) +
                            1e-8
                    )
              )
}
nnGen <- function(data, numCars, batch_size = 128){
  i <- 1
  function(){
    samples <- array(0, dim = c(batch_size, 40))
    targets <- array(0, dim = c(batch_size, 20))
    dataSub <- data[[i]]
    if(i >= numCars){
      i <<- 1
    } else {
      i <<- i + 1
    }
    t <- sample(10:(nrow(dataSub) - 10), batch_size)
    for(j in 1:batch_size){
      for(k in 1:2){
        samples[j, (k-1)*10+1:10] <- dataSub[(-9:0) + t[j], k] -  dataSub[(-10:-1) + t[j], k]
      }
      for(k in 3:4){
        samples[j, (k-1)*10+1:10] <- dataSub[(-9:0) + t[j], k] 
      }
      targets[j,1:10] <- dataSub[(1:10) + t[j], 1] - dataSub[(0:9) + t[j], 1]
      targets[j,11:20] <- dataSub[(1:10) + t[j], 2] - dataSub[(0:9) + t[j], 2]
    }
    
    list(samples, targets)
  }
}
trainSet <- sample(1:2000, 1600)
train_gen <- nnGen(data[trainSet], numCars = 1600)
test_gen <- nnGen(data[-trainSet], numCars = 400)

nn %>% compile(
  optimizer = 'adam',
  loss = 'mse'
)

history <- nn %>% fit_generator(
  train_gen,
  steps_per_epoch = 50,
  epochs = 25,
  validation_data = test_gen,
  validation_steps = 25
)
plot(history)



a <- test_gen()
pred <- predict(nn, a[[1]])
p <- rbind(exp(pred[1,,1]) * cos(pred[1,,2] + pi/2), exp(pred[1,,1]) * sin(pred[1,,2] + pi/2), 
      exp(a[[2]][1,,1]) * cos(a[[2]][1,,2] + pi/2), exp(a[[2]][1,,1]) * sin(a[[2]][1,,2] + pi/2))
apply(p, 2, function(x) sqrt((x[1] - x[3])^2 + (x[2] - x[4])^2)) %>% cumsum() %>% mean()

dataFC <- list()
for(i in 1:873){
  carsAug[carsAug$ID == id$idfc[i],] -> carI
  carI[2:min(501, nrow(carI)),] %>%
    mutate(d = delta - pi/2) %>%
    select(x, y, v, d) -> carI
  dataFC[[i]] <- as.matrix(carI)
}
for(i in 1:500){
  carsAug[carsAug$ID == noStopChanged[i],] -> carI
  carI[2:min(501, nrow(carI)),] %>%
    mutate(d = delta - pi/2) %>%
    select(x, y, v, d) -> carI
  dataFC[[i+873]] <- as.matrix(carI)
}




predictionData <- buildTestX(dataFC)
RNNResults <- data.frame()
for(i in 1:1373){
  posForecast <- predict(rnn, predictionData[[i]])
  if(i <= 873){
    carID <- id$idfc[i]
  } else {
    carID <- noStopChanged[i]
  }
  for(s in 1:30){
    true <- dataFC[[i]][10*s + 0:30, 1:2]
    pred <- posForecast[s, ]
    error <- mPerFoot * sqrt((true[2:31,1] - true[1:30, 1] - pred[1:30])^2 + (true[2:31,2] - true[1:30, 2] -  pred[31:60])^2)
    RNNResults <- rbind(RNNResults,
                        data.frame(h = 1:30,
                                   S = s,
                                   id = carID,
                                   model = 'RNN',
                                   error = error))
  }
}

}

