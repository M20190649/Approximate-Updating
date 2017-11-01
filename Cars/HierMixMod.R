library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
source('mixtureMCMC.R')

Y600 <- read.csv('Y600.csv')

Y600 %>% 
  group_by(ID) %>%
  filter(min(v) > 0.1 & changed == FALSE) %>%
  .$ID %>%
  unique() -> noChange

N <- 150
idSubset <- sample(noChange, N)
data <- list()
for(i in 1:N){
  Y600 %>%
    filter(ID == idSubset[i]) %>%
    mutate(n = seq_along(v),
           vl = ifelse(n == 1, 0, lag(v)),
           v = v - lag(v),
           d = delta - pi/2) %>%
    filter(n > 1) %>% 
    select(v , d) %>%
    as.matrix() -> data[[i]]
}

reps <- 20000
mixDraws <- mixtureMCMC(data, reps)


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
#  filter(iter > 4000) %>%
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