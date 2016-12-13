library(readr)
library(ggplot2)
library(GGally)
library(tidyr)
library(moments)
library(MASS)
library(mixtools)


draws <- read_csv("exp_sm_MCMC.csv")
ggpairs(draws[seq(1, 10000, 10),])
drawsl = gather(draws, variable, draw) 
ggplot(drawsl, aes(draw)) + geom_density() + facet_wrap(~variable, ncol = 5, scales = "free")

calc_moments = function(x){
  require(moments)
  return(c(mean(x), var(x), skewness(x), kurtosis(x)))
}

moments = apply(draws, 2, calc_moments)
rownames(moments) = c("Mean", "Variance", "Skewness", "Kurtosis")
moments
#Use IG for sigmasq and Gaussian for initial states - same as conditionals
#IG method of moment parameters here


a1 = fitdistr(draws$Alpha, "beta", list(shape1 = 1.5, shape2 = 1.5))
- 2 * a1$loglik + 2 * log(10000)
a2 = normalmixEM(draws$Alpha, k = 2)
- 2 * a2$loglik + 6 * log(10000)
a3 = normalmixEM(draws$Alpha, k = 3)
- 2 * a3$loglik + 9 * log(10000)
a4 = normalmixEM(draws$Alpha, k = 4)
- 2 * a4$loglik + 12 * log(10000)
a5 = normalmixEM(draws$Alpha, k = 5)
- 2 * a5$loglik + 15 * log(10000)

d1 = fitdistr(draws$Delta, "beta", list(shape1 = 0.5, shape2 = 0.5))
-2 * d1$loglik + 2 * log(10000)
d2 = normalmixEM(draws$Delta, k = 2)
-2 * d2$loglik + 6 * log(10000)


#Choose two component normals for both
support = seq(0, 1, 0.001)
alphadens1 = dbeta(support, a1$estimate[1], a1$estimate[2])
alphadens2 = a2$lambda[1] * dnorm(support, a2$mu[1], a2$sigma[1]) + a2$lambda[2] * dnorm(support, a2$mu[2], a2$sigma[2])
ggplot() + geom_density(aes(draws$Alpha)) + geom_line(aes(support, alphadens1), colour = "red") + 
  geom_line(aes(support, alphadens2), colour = "blue")
deltadens1 = dbeta(support, d1$estimate[1], d1$estimate[2])
deltadens2 = d2$lambda[1] * dnorm(support, d2$mu[1], d2$sigma[1]) + d2$lambda[2] * dnorm(support, d2$mu[2], d2$sigma[2])
ggplot() + geom_density(aes(draws$Delta)) + geom_line(aes(support, deltadens1), colour = "red") + 
  geom_line(aes(support, deltadens2), colour = "blue") 

#Next find copulas
library(VineCopula)
pobdraws = pobs(draws)
Vine = RVineStructureSelect(pobdraws, indeptest = TRUE)
saveRDS(Vine, "Vine.rds")
