########## TO DO #############
# Replace K with a Direchlet Process & Slice Sampler in the MCMC
# Write VB algorithm with Graves' mixture derivatives - implement in arUpdaterMix
# Run Slurm code

# Consider new online implementations: 
# add new data every 100ms or batches slightly more often 
# Swap fit till converge, add more data, to add new data, switch prior asap.

# For the supervisors:
# Write more things.


library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(GGally)
source('mixtureMCMC.R')
sourceCpp('arvdMCMC.cpp')
sourceCpp('hierarchical.cpp')
id <- readRDS('carsID.RDS')
data <- readRDS('MCMCData.RDS')

read.csv('carsAug.csv') %>%
  select(ID, relX, dist, relv, reldelta) %>%
  ungroup() -> carsAug
colnames(carsAug) <- c('ID', 'x', 'y', 'v', 'delta')

carsAug %>% 
  group_by(ID) %>%
  filter(min(v) != 0 & !ID %in% id$splinesID) %>%
  .$ID %>%
  unique() -> noStop

HierMixMod{

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
hyper$alpha <- rep(1, K)
for(i in 1:N){
  draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:2, 1), pi = rep(1/K, K))
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
for(k in 1:K){
  mat <- matrix(0, 1000, 6)
  for(i in (burn*reps/thin+1):(reps/thin)){
    meanvec <- mixDraws$draws[[1]][[k]]$mean[i,]
    sdvec <- sqrt(diag(solve(mixDraws$draws[[1]][[k]]$varInv[,,i])))
    w <- 0
    for(j in 1:N){
      w <- w + mixDraws$draws[[j+1]]$pi[i, k] / N
    }
    for(j in 1:2){
      mat[,j] <- mat[,j] + w * dlnorm(support[,j], meanvec[j], sdvec[j]) / (reps / (2 * thin)) 
    } 
    for(j in 3:6){
      mat[,j] <- mat[,j] + w * dnorm(support[,j], meanvec[j], sdvec[j]) / (reps / (2 * thin)) 
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
  summarise(dens = sum(dens)) %>%
  ggplot() + geom_line(aes(support, dens)) +
  geom_line(data=densities, aes(support, dens, colour = factor(group))) +
  facet_wrap(~var, scales = 'free', ncol = 1) + 
  labs(title = 'Hierarchical Model') + 
  theme(legend.position = 'none')
}

OtherMods{
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
  
  noMixDraws$draws[[1]]$mean %>%
    as.data.frame() %>%
    cbind(iter = 1:(reps/thin)) %>%
    mutate(V1 = exp(V1), V2 = exp(V2)) %>%
    rename(sigSq_eps = V1, sigSq_eta = V2, phi1 = V3, phi2 = V4, gamma1 = V5, gamma2 = V6) %>%
    gather(var, draw, -iter) %>%
    mutate(var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
    filter(iter > 100) %>%
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
  colnames(densities) <- vars
  densities %>% 
    gather(var, dens) %>%
    cbind(support = c(support)) %>%
    mutate(var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
    ggplot() + geom_line(aes(support, dens)) + facet_wrap(~var, scales = 'free')
    
  
  draws <- c(-5, -5, 0, 0, 0, 0)
  hyper <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(10, 6)))
  hyper$var <- diag(solve(hyper$varInv))
  noHierDraws <- noHierMCMC(data, reps, draws, hyper, thin = 10, stepsize = 0.01)
  saveRDS(noHierDraws, 'noHierN2000.RDS')
  
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
  noHierLam <- c(noHierMean, chol(noHierVar))
  for(i in (reps/(2*thin)+1):(reps/thin)){
    noHierMean <- noHierMean + noHierDraws$draws[i,] / (reps / (2 * thin))
  }
  
}

VB{
  
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
posAC <- NULL
vars <- c('sigSq_eps', 'sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2') 
measures <- c('ac1a', 'ac2a', 'ac1d', 'ac2d', 'siga', 'sigd')
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
  
  x <- rlnorm(1000000, -5, sqrt(5))
  a <- mean(x)^2 / var(x) + 2
  b <- mean(x) * (mean(x)^2 / var(x) + 1)
  hyper <- c(a, b, a, b, 0, 10, 0, 10, 0, 10, 0, 10)

  
  MCMC <- ARVD_MCMC(data, hyper, 10000, 2)
  mean <- colMeans(MCMC[5001:10000,])
  ac <- rep(0, 6)
  for(l in 5001:10000){
    ac[1] <- ac[1] + MCMC[l, 3] / (1 - MCMC[l, 4])
    ac[2] <- ac[2] + (MCMC[l, 3]^2 / (1 - MCMC[l, 4]) + MCMC[l, 4])
    ac[3] <- ac[3] + MCMC[l, 5] / (1 - MCMC[l, 6])
    ac[4] <- ac[4] + (MCMC[l, 5]^2 / (1 - MCMC[l, 6]) + MCMC[l, 6])
    ac[5] <- ac[5] + sqrt((1-MCMC[l, 4]) * MCMC[l, 1] / ((1+MCMC[l, 4])*(1-MCMC[l, 3]-MCMC[l, 4])*(1+MCMC[l, 3]-MCMC[l, 4])))
    ac[6] <- ac[6] + sqrt((1-MCMC[l, 6]) * MCMC[l, 2] / ((1+MCMC[l, 6])*(1-MCMC[l, 5]-MCMC[l, 6])*(1+MCMC[l, 5]-MCMC[l, 6])))
  }
  if(any(is.na(ac))){
    print(paste('NA', i))
  }
  posMeans <- rbind(posMeans, data.frame(mean = mean, var = vars, ID = id$idSubset[i]))
  posAC <- rbind(posAC, data.frame(ac = ac/5000, var = measures, ID = id$idSubset[i]))
  if(i %% 50 == 0){
    print(i)
  }
}

posAC %>%
  spread(var, ac) %>%
  mutate(pac2a = (ac2a - ac1a^2) / (1 - ac1a^2),
         pac2d = (ac2d - ac1d^2) / (1 - ac1d^2)) %>%
  select(ID, ac1a, pac2a, ac1d, pac2d, siga, sigd) %>%
  rename(pac1a = ac1a, pac1d = ac1d) -> autocorrels

autocorrels %>%
  filter(siga < 0.5) %>%
  ggpairs(columns = c(6, 2, 3, 7, 4, 5))



posMeans %>%
  mutate(var = as.character(var),
         var = ifelse(var == 'log_sigSq_eps', 'sigSq_eps', 
                      ifelse(var == 'log_sigSq_eta', 'sigSq_eta', var)),
         var = factor(var, levels = c('sigSq_eps', 'phi1', 'phi2', 'sigSq_eta', 'gamma1', 'gamma2'))) %>%
  ggplot() + geom_density(aes(mean)) + 
  facet_wrap(~var, scales = 'free', ncol = 1) +
  labs(title = 'Single Model Means (Kernel Density Estimate)')
}

sliceSampler{
  set.seed(1)
  N <- 200
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
  reps <- 10000
  K <- 10
  thin <- 1
  draws <- list(list())
  hyper <- list()
  for(k in 1:K){
    hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(10, 6)), v = 6, scale = diag(1, 6))
    draws[[1]][[k]] <- list(mean = c(-5, -5, 0, 0, 0, 0), varInv = diag(10, 6))
  }
  draws[[1]]$pi <- rep(1/K, K)
  draws[[1]]$alpha <- 5
  for(i in 1:N){
    draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:K, 1))
  }
    
  sliceDraws <- sliceSampler(data, reps, draws, hyper, thin, K, 'gaussian', 0.01)
}

forecasts{
  # This should be parallelised and put onto slurm.
  increment <- FALSE
  S <- 10
  maxT <- 300
  K <- 6
  sSeq <- seq(S, maxT, S)
  results <- data.frame()
  methods <- c('None', 'No Hier', 'Single', 'Mixture')
  prior <- list()
  starting <- matrix(c(-5, -5, rep(0, 4), c(chol(diag(0.5, 6)))), ncol = 1)
  prior[[1]] <- matrix(c(-5, -5, rep(0, 4), c(chol(diag(c(5, 5, 10, 10, 10, 10))))), ncol = 1)
  prior[[2]] <- matrix(noHierLam, ncol = 1)
  prior[[3]] <- matrix(noMixLam, ncol = 1)
  fit <- prior
  
  datafc <- list()
  for(i in seq_along(id$idfc)){
    carsAug %>%
    filter(ID == id$idfc[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           a = v - lag(v),
           d = delta - pi/2) %>%
    filter(n > 1 & n <= 501) %>% 
    select(a , d) %>%
    as.matrix() -> datafc[[i]]
  }
  saveRDS(datafc, 'ForecastData.RDS')
  
  
  for(i in seq_along(id$idfc)){
    
    # Extract Data
    datafc[[i]] -> data
    
    # Set forcast supports
    asup <- seq(0.8*min(data[,1]), 1.25*max(data[,1]), length.out=1000)
    dsup <- seq(0.8*min(data[,2]), 1.25*max(data[,2]), length.out=1000)
      
    # Incrementally add data to VB fits
    for(s in seq_along(sSeq)){
      if(sSeq[s] > nrow(data)){
        break
      }
      if(s == 1 | !increment){
        dat <- data[1:sSeq[s],]
      } else {
        dat <- data[(sSeq[s-1]+1):sSeq[s],]
      }
      # Update posterior approximations - Or re-estimate new ones from scratch
      if(increment){
        fit <- fitCarMods(dat, fit, increment, NULL)
      } else {
        fit <- fitCarMods(dat, prior, increment, starting)
      }

      # Extract Lower Triangular Matrices
      L <- NULL
      for(k in 1:3){
        L <- rbind(L, t(matrix(fit[[k]][7:42], 6)))
      }
      means <- cbind(fit[[1]][1:6],
                     fit[[2]][1:6],
                     fit[[3]][1:6]) 
      #sigma4 <- list()
      #for(k in 1:K){
      #  u <- matrix(fit[[4]][1:36 + 6*K + 36*(k-1)], 6)
      #  sigma4[[k]] <- t(u) %*% u
      #}
      #weights <- fit[[4]][42*K+1:K]
      
      densities <- evalFcDens(data[-1:S + sSeq[s], ], means, L, 1000, S, asup, dsup)
      
      # Grab logscores for each method, h, and variable.
      for(k in 1:3){
        for(h in 1:S){
          aindex <- min(which(asup > data[sSeq[s]+h,1]))
          alogscore <- log(densities[(h-1)*1000 + aindex, 1, k])
          dindex <- min(which(dsup > data[sSeq[s]+h,2]))
          dlogscore <- log(densities[(h-1)*1000 + dindex, 2, k])
          # Attach results
          results <- rbind(results, 
                           data.frame(logscore = c(alogscore, dlogscore),
                                      variable = c('a', 'd'),
                                      method = methods[k],
                                      S = sSeq[s],
                                      h = h,
                                      id = id$idfc[i]))
        }
      }
    }
    print(i)
  }
  
  
  results %>%
    group_by(S, method, variable) %>%
    summarise(ls = mean(logscore, na.rm = TRUE)) %>%
    ggplot() + geom_line(aes(S, ls, colour = method)) + facet_wrap(~variable, scales = 'free')
  
}


