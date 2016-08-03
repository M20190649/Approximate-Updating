##AR(1) Model
library(Rcpp)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RcppArmadillo)
library(microbenchmark)
library(MASS)
library(parallel)
sourceCpp("AR1c.cpp")

T <- 50
J <- 10
sigma <- 2
mu <- 2
rho <- 0.5

y <- genAR1(T+J, c(mu, rho, sigma))


X <- as.matrix(cbind(rep(1,59), y[1:59]))
InvXX <- solve(t(X)%*%X)
b <- InvXX%*%t(X)%*%y[2:60]                   
v <- T+J-3

#exact posterior
sigsqexact <- t(y[2:60]-X%*%b)%*%(y[2:60]-X%*%b)/(T+J-3)
sigmaexact <- 1/sqrt(gen_gamma(2000, v/2, 1/(sigsqexact*v/2)))
sigtrans <- transformsigma(sigmaexact, InvXX)
betaexact <- t(apply(sigtrans, 3, gen_mvrnorm, n=1, mu=b))
exact <- data.frame(mu = betaexact[,1], psi = betaexact[,2], sigma = sigmaexact, method="Exact")



#The ABC Algorithm

index.bottom.N = function(x, N){
  o = order(x, na.last=FALSE, decreasing=TRUE)
  o.length = length(o)
  o[((o.length-N+1):o.length)]
}

abcpar <- function(theta, cores, n, T, yss) { ##Parallelisation returns a 20% performance boost, 4 cores seems optimal
  theta.list <- list()
  n <- nrow(theta)
  for(i in 1:cores) { #split theta into a list, one for each core to reduce parallelisation overhead
    theta.list[[i]] <- theta[((i-1)*(n/cores)+1):(i*(n/cores)),]
  }
  distl <- unlist(mclapply(theta.list, ABC_Loop, n/cores, T, yss=yss, mc.cores=cores)) #run the mclapply on each element on list
  return(distl)
} 
#ABC_Loop function does all the work


yss <- sum_stats(y[51:60]) #Calculate summarry statistics of the updated data
InvXXprior <- solve(t(X[1:49,])%*%X[1:49,]) #run through exact inference on the original data to get a prior
bprior <- InvXXprior%*%t(X[1:49,])%*%y[2:50]
sigsqprior <- as.numeric(t(y[2:50]-X[1:49,]%*%bprior)%*%(y[2:50]-X[1:49,]%*%bprior)/(T-3))
sigmaprior <- 1/sqrt(gen_gamma(200000, (T-3)/2, 1/(sigsqprior*(T-3)/2)))
sigtransp <- transformsigma(sigmaprior, InvXXprior) ## Sigma * inverse(X'X)
betaprior <- t(apply(sigtransp, 3, gen_mvrnorm, n=1, mu=bprior)) ##Draw mu from the sigma
theta <- cbind(betaprior, sigmaprior)
dist <- abcpar(theta, 4, n, J, yss) ##The ABC Algorithm
accept <- theta[index.bottom.N(dist, 2000),]
abc <- data.frame(mu = accept[,1], psi = accept[,2], sigma = accept[,3], method = "ABC")



##the prior
prior <- data.frame(mu = theta[1:2000,1], psi = theta[1:2000,2], sigma = theta[1:2000,3], method = "Prior")

#setting up the variational algorithm
k = 2
sigsqhat <- var(y[2:60])*59/(59-k)
shape = (T+J)/2
scale = (T+J-1-k)*sigsqhat/(2*(1-k/(T+J)))
lambda = solve(XX)*shape/scale
betavar <- mvrnorm(2000, b, lambda)
sigmavar <- 1/sqrt(rgamma(2000, shape, rate=scale))
variational <- data.frame(mu = betavar[,1], psi = betavar[,2], sigma = sigmavar, method="Variational")


AR1 <- rbind(prior, exact, variational, abc)
#AR1 <- mutate(AR1, mu.actual = mu/(1-psi))
AR1l <- gather(AR1, theta, density, -method)

ggplot(data=AR1l, aes(x=density, colour=method)) + geom_density() + facet_wrap(~theta, scales="free")

