library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
source('mixtureMCMC.R')
splinesID <- readRDS('splinesID.RDS')

read.csv('carsAug.csv') %>%
  select(ID, relX, dist, relv, reldelta) %>%
  ungroup() -> carsAug
colnames(carsAug) <- c('ID', 'x', 'y', 'v', 'delta')

carsAug %>% 
  group_by(ID) %>%
  filter(min(v) != 0 & !ID %in% splinesID) %>%
  .$ID %>%
  unique() -> noStop

set.seed(1)
N <- 2500
idSubset <- sample(noStop, N)
data <- list()
for(i in 1:N){
  carsAug %>%
    filter(ID == idSubset[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           v = v - lag(v),
           d = delta - pi/2) %>%
    filter(n > 1 & n < 501) %>% 
    select(v , d) %>%
    as.matrix() -> data[[i]]
}
saveRDS(data, 'mixmod.RDS')

reps <- 1000
K <- 2


draws <- list(list(list(mean = c(-6, -6, 0.5, 0, 0.3, 0.1), varInv = diag(10, 6)),
                   list(mean = c(-4, -4, -0.5, -0.2, 0, 0), varInv = diag(10, 6))))
hyper <- list()
for(k in 1:K){
  hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(5, 6)),  v = 6, scale = diag(0.5, 6))
}
for(i in 1:N){
  draws[[i+1]] <- list(theta = c(-5, -5, 0, -0.1, 0.15, 0.05), k = sample(1:2, 1))
}


mixDraws <- mixtureMCMC(data, reps, draws, hyper, 10, K, 'gaussian')


mapGroup <- NULL
for(i in 1:N){
  kDraws <- mixDraws[[i+1]]$k[(reps/2+1):reps]
  mapGroup <- rbind(mapGroup, data.frame(group = mean(kDraws), ID = idSubset[i]))
}
ggplot(mapGroup) + geom_histogram(aes(group)) + theme_bw()

muK <- NULL
for(i in 1:K){
  mixDraws[[1]][[i]]$mean %>%
    cbind(iter = 1:reps) %>%
    as.data.frame() %>%
    mutate(group = i)  -> temp
  colnames(temp) <- c('log_sigSq_eps', 'log_sigSq_eta', 'phi1', 'phi2', 'gamma1', 'gamma2', 'iter', 'group')
  muK <- rbind(muK, temp)
}


muK %>%
  gather(var, draw, -iter, -group) %>%
  mutate(var = factor(var, levels = c('log_sigSq_eps', 'phi1', 'phi2', 'log_sigSq_eta', 'gamma1', 'gamma2'))) %>%
  filter(iter > 2000 & iter %% 50 == 0 ) %>%
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_grid(var ~ group, scales = 'free') + theme_bw()

