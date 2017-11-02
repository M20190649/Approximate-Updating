library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
source('mixtureMCMC.R')

read.csv('carsAug.csv') %>%
  select(ID, relX, dist, relv, reldelta)  -> cars
colnames(cars) <- c('ID', 'x', 'y', 'v', 'delta')

cars %>% 
  group_by(ID) %>%
  filter(min(v) > 0.1) %>%
  .$ID %>%
  unique() -> noChange


cars %>% 
  filter(ID < 40) %>%
  group_by(ID) %>% 
  mutate(n = seq_along(v),
         vl = ifelse(n == 1, 0, lag(v)),
         vd = v - vl) %>% filter(n > 1) %>% 
  ggplot() + geom_path(aes(n, v + ID, group = ID))


N <- 250
idSubset <- sample(noChange, N)
data <- list()
for(i in 1:N){
  cars %>%
    filter(ID == idSubset[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           v = v - lag(v),
           d = delta - pi/2) %>%
    filter(n > 1 & n < 501) %>% 
    select(v , d) %>%
    as.matrix() -> data[[i]]
}

reps <- 20000
K <- 2

mixDraws <- mixtureMCMC(data, reps, K, 'gaussian')


mapGroup <- NULL
for(i in 1:N){
  kDraws <- mixDraws[[i+1]]$k[(reps/2+1):reps]
  mode <- which.max(table(kDraws))
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
  filter(iter > 10000) %>%
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_grid(var ~ group, scales = 'free') + theme_bw()

