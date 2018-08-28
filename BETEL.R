library(tidyverse)

N <- 250
z <- rnorm(N, 0.5, 1)
alpha <- 0
beta <- 1
y <- alpha + beta * z
for(i in 1:N){
  u <- runif(1)
  if(u < 0.5){
    y[i] <- y[i] + rnorm(1, 0.75, 0.75)
  } else {
    y[i] <- y[i] + rnorm(1, -0.75, 1.25)
  }
}
nu <- 0.5 * (3 * 0.75 * 0.75^2 + 0.75^3) + 0.5 * (3 * (-0.75) * 1.25^2 + (-0.75)^3)
lam <- c(1, 2, 3)

hLinReg <- function(lam, g){
  mean(exp(lam %*% t(g)))
}

density <- function(y, z, theta, lambdaStart = rep(0.1, 3)){
  alpha <- theta[1]
  beta <- theta[2] 
  nu <- theta[3]
  
  g1 <- y - alpha - beta * z
  g2 <- z * g1
  g3 <- g1^3 - nu
  g <- cbind(g1, g2, g3)
  
  lambdaHat <- nlm(hLinReg, lambdaStart, g = g)$estimate
  

  exponent <- c(exp(lambdaHat %*% t(g)))
  
  list(dens = sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent)),
       lambdaHat = lambdaHat)
}



reps <- 10000
betelMCMC <- matrix(0, reps, 3)
draw <- c(0, 1, -1)
accept <- 0
dens <- density(y, z, draw)
oldDens <- dens$dens
lambdaStart <- dens$lambdaHat
lam <- matrix(0, reps, 3)

for(i in 1:reps){
  
  candidate <- draw + 0.1 * rt(3, 2.5)
  dens <- density(y, z, candidate, lambdaStart)
  canDens <- dens$dens
  alpha <- min(1, exp(canDens - oldDens))
  if(runif(1) < alpha){
    draw <- candidate
    lambdaStart <- dens$lambdaHat
    oldDens <- canDens
    accept <- accept + 1
  }
  betelMCMC[i, ] <- draw
  lam[i, ] <- lambdaStart
}

accept / reps

betelMCMC %>%
  as.tibble() %>%
  rename(alpha = V1, beta = V2, nu = V3) %>%
  mutate(iter = seq_along(alpha)) %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free')
  

betelMCMC %>%
  as.tibble() %>%
  rename(alpha = V1, beta = V2, nu = V3) %>%
  mutate(iter = seq_along(alpha)) %>%
  filter(iter > 2500) %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(draw), stat = 'density') + facet_wrap(~var, scales = 'free')


trueDensity <- function(y, z, theta){
  dens <- 0.5 * 1 / sqrt(2 * pi * 0.75^2) * exp(-(y - theta[1] - theta[2] * z - 0.75)^2 / (2 * 0.75^2)) +
    0.5 * 1 / sqrt(2 * pi * 1.25^2) *  exp(-(y - theta[1] - theta[2] * z + 0.75)^2 / (2 * 1.75^2))
  sum(log(dens))
}

reps <- 10000
trueMCMC <- matrix(0, reps, 2)
draw <- c(0, 1)
accept <- 0
oldDens <- trueDensity(y, z, draw)

for(i in 1:reps){
  candidate <- draw + 0.1 * rt(2, 2.5)
  canDens <- trueDensity(y, z, candidate)
  alpha <- min(1, exp(canDens - oldDens))
  if(runif(1) < alpha){
    draw <- candidate
    oldDens <- canDens
    accept <- accept + 1
  }
  trueMCMC[i, ] <- draw
}

trueMCMC %>%
  as.tibble() %>%
  rename(alpha = V1, beta = V2) %>%
  mutate(iter = seq_along(alpha)) %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free')

trueMCMC %>%
  as.tibble() %>%
  rename(alpha = V1, beta = V2) %>%
  mutate(iter = seq_along(alpha)) %>%
  filter(iter > 2500) %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(draw), stat = 'density') + facet_wrap(~var, scales = 'free')

MCMCdraws <-
  rbind(trueMCMC %>%
          as.tibble() %>%
          rename(alpha = V1, beta = V2) %>%
          mutate(iter = seq_along(alpha)) %>%
          filter(iter > 2500) %>%
          gather(var, draw, -iter) %>%
          mutate(method = 'MCMC-true'),
        betelMCMC %>%
          as.tibble() %>%
          rename(alpha = V1, beta = V2, nu = V3) %>%
          mutate(iter = seq_along(alpha)) %>%
          filter(iter > 2500) %>%
          gather(var, draw, -iter) %>%
          mutate(method = 'MCMC-empirical'))

ggplot(MCMCdraws) + geom_line(aes(draw, colour = method), stat = 'density') + 
  geom_vline(data = data.frame(nu = nu, method = 'MCMC-true', var = 'nu'), aes(xintercept = nu, colour = method)) + 
  facet_wrap(~var, scales = 'free')

reparamDeriv <- function(y, z, epsilon, lambda, lamStart){
  # Transform
  U <- matrix(lambda[4:12], 3)
  theta <- lambda[1:3] + t(U) %*% epsilon
  
  # Set of transform derivatives
  dtdl <- matrix(0, 12, 3)
  diag(dtdl[1:3, 1:3]) <- 1
  
  dtdl[4, 1] <- dtdl[7, 2] <- dtdl[10, 3] <- epsilon[1]
  dtdl[8, 2] <- dtdl[11, 3] <- epsilon[2]
  dtdl[12, 3] <- epsilon[3]
  
  djdl <- rep(0, 12)
  djdl[c(4, 8, 12)] <- 1 / lambda[c(4, 8, 12)]
  
  # Calculate dpdt components
  n <- length(y)
  alpha <- theta[1]
  beta <- theta[2] 
  nu <- theta[3]
  
  g1 <- y - alpha - beta * z
  g2 <- z * g1
  g3 <- g1^3 - nu
  g <- cbind(g1, g2, g3)
  minima <- nlm(hLinReg, rep(0, 3), g = g, hessian = TRUE)
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  # export numerical work to C++
  gradients <- matrixCalculations(y, z, alpha, beta, g, dh2dlam2, lambdaHat, exponent)
  dpdt <- gradients$grad
  logp <- gradients$val
  
  # Calculate the derivative of the  ELBO
  dELBO <- dtdl %*% dpdt + djdl
  
  list(grad = dELBO, val = logp)
}


fitBETEL <- function(y, z, lambda, model, S = 25, dimTheta = 3, zEpsilon = TRUE,
                  maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01){
  if(!is.matrix(lambda)){
    lambda <- matrix(lambda, ncol = 1)
  }
  dimLambda <- length(lambda)
  
  lambda <- cbind(lambda, matrix(0, dimLambda, maxIter))
  
  

  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- V <- numeric(dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    eval <- numeric(S)
    grad <- matrix(0, dimLambda, S)
    q <- numeric(S)
    unif <- matrix(runif(S*dimTheta), ncol = dimTheta)
    for(s in 1:S){
      epsilon <- qnorm(unif[s,])
      q <- sum(dnorm(epsilon, log = TRUE))
     
      logpj <- model(y, z, epsilon, as.matrix(lambda[,iter]), rep(0, 3))
      eval[s] <- logpj$val
      grad[,s] <- logpj$grad
    }
    
    eval[eval == -Inf] = NA
    gradient <- rowMeans(grad, na.rm = TRUE)
    gradientSq <- rowMeans(grad^2, na.rm = TRUE)
    LB[iter] <- mean(eval - q, na.rm=TRUE) 
    
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mstar <- M / (1 - beta1^iter)
    Vstar <- V / (1 - beta2^iter)
    update <- alpha * Mstar / (sqrt(Vstar) + e)
    if(any(is.na(update))){
      print('Break')
      break
    }
    lambda[,iter + 1] <- lambda[, iter] + update
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 10 == 0){
        print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', min(iter-1, maxIter), ' ELBO: ', meanLB))
  return(list(lambda=lambda, LB = LB[1:min(iter-1, maxIter)], iter = min(maxIter, iter-1)))
}

lambda <- c(0, 0, 0, diag(0.1, 3))
fit <- fitBETEL(y, z, lambda, reparamDeriv, alpha = 0.03, S = 25, maxIter = 1000, threshold = 0.01)


qplot(1:fit$iter, fit$LB[1:fit$iter], geom = 'line')


best <- which.max(fit$LB[1:fit$iter])

vbDens <- vbDensity(list(mean = fit$lambda[1:3, best], U = fit$lambda[4:12, best]),
                     rep('identity', 3),
                     c('alpha', 'beta', 'nu'))
vbDens$method <- 'VB-Reparam'


ggplot(MCMCdraws) + geom_line(aes(draw, colour = method), stat = 'density') +
  geom_vline(data = data.frame(nu = nu, method = 'MCMC-true', var = 'nu'), aes(xintercept = nu, colour = method)) + 
  geom_line(data = vbDens, aes(support, density, colour = method)) + 
  facet_wrap(~var, scales = 'free')




fit$lambda[,1:fit$iter] %>%
  t() %>%
  as.data.frame() %>%
  cbind(iter = 1:fit$iter) %>% 
  gather(var, lambda, -iter) %>%
  group_by(var) %>%
  filter(!all(lambda == 0)) %>%
  ungroup() %>%
  mutate(var = rep(c('alpha', 'beta', 'nu', 'L11', 'L21', 'L22', 'L31', 'L32', 'L33'), rep(fit$iter, 9)),
         var = factor(var, levels = c('alpha', 'beta', 'nu', 'L11', 'L21', 'L22', 'L31', 'L32', 'L33'))) %>%
  ggplot() + geom_line(aes(iter, lambda)) + facet_wrap(~var, scales = 'free') + theme_bw()



### Stein
Rcpp::sourceCpp('misc/particleVB.cpp')
steinDeriv <- function(y, z, theta, lamStart){
  
  # Calculate dpdt components
  n <- length(y)
  alpha <- theta[1]
  beta <- theta[2] 
  nu <- theta[3]
  
  g1 <- y - alpha - beta * z
  g2 <- z * g1
  g3 <- g1^3 - nu
  g <- cbind(g1, g2, g3)
  minima <- nlm(hLinReg, lamStart, g = g, hessian = TRUE)
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  # export numerical work to C++
  gradients <- matrixCalculations(y, z, alpha, beta, g, dh2dlam2, lambdaHat, exponent)
  dpdt <- gradients$grad
  logp <- gradients$val
  
  list(grad = dpdt, val = logp)
}

steinVBWrapper <- function(y, z, particles, N, h = 0.01, alpha = 1e-3, maxIter = 1000, threshold = 1e-3){
  
  iter <- 1
  LB <- rep(0, maxIter)
  diff <- threshold + 1
  M <- V <- rep(0, N)
  e <- 1e-8
  oldLB <- 0
  
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    
    phiHat <- matrix(0, nrow(particles), ncol(particles))
    logP <- matrix(0, N, ncol(particles))
    elbo <- rep(0, N)
    
    for(i in 1:N){
      stein <- steinDeriv(y, z, particles[i, ],  rep(0, 3))
      elbo[i] <- stein$val
      logP[i, ] <- stein$grad
    }
    
    for(i in 1:N){
      for(j in 1:N){
        kernel <- RBFKernel(particles[i, ], as.matrix(particles[j, ]), h)
        phiHat[i, ] <- phiHat[i, ] + 1/N * (kernel$val * logP[j, ] + kernel$grad)
      }
    }
    LB[iter] <- mean(elbo)
    
    particles <- particles + alpha * phiHat
    
    
    if(iter %% 5 == 0){
      meanLB <- mean(LB[iter:(iter- 4)])
      diff <- abs(meanLB - oldLB)
      oldLB <- meanLB
    }
    if(iter %% 25 == 0){
      print(paste0('Iteration: ', iter, ', ELBO: ', meanLB))
    }
    
    iter <- iter + 1
  }
  print(paste0('Converged at Iteration: ', iter - 1, ' at ELBO: ', meanLB))
  particles
}
N <- 100
particles <-matrix(rnorm(3 * N, 0, 0.5), ncol = 3) +  betelMCMC[sample(5001:10000, N), ]

initial <- data.frame(particle = c(particles), var = rep(c('alpha', 'beta', 'nu'), rep(N, 3)))


ggplot(MCMCdraws) + geom_line(aes(draw, colour = method), stat = 'density') +
  geom_vline(data = data.frame(nu = nu, method = 'MCMC-true', var = 'nu'), aes(xintercept = nu, colour = method)) + 
  geom_line(data = initial, aes(particle), stat = 'density') + 
  labs(colour = 'MCMC Likelihood') +
  facet_wrap(~var, scales = 'free')

h <- 0.1


fitStein <- steinVBWrapper(y, z, particles, N, maxIter = 1000, threshold = 0.05)

steinVB <- data.frame(particle = c(fitStein), var = rep(c('alpha', 'beta', 'nu'), rep(N, 3)), method = 'VB-Stein')
  

cbPalette <- c("#000000", "#E69F00", "#2684B9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
MCMCdraws %>%
  mutate(method = ifelse(method == 'MCMC-true', 'MCMC-parametric', method)) %>%
  ggplot() + geom_line(aes(draw, colour = method), stat = 'density') +
  geom_vline(data = data.frame(nu = nu, method = 'MCMC-parametric', var = 'nu'), aes(xintercept = nu, colour = method)) + 
  geom_line(data = steinVB, aes(particle, colour = method), stat = 'density') + 
  geom_line(data = vbDens, aes(support, density, colour = method)) +
  theme_bw() + 
  facet_wrap(~var, scales = 'free', label = label_parsed) + 
  scale_colour_manual(values=cbPalette) + 
  theme(legend.position = 'bottom') + 
  labs(x = 'support')




vbmu <- fit$lambda[1:3, best]
vbU <- matrix(fit$lambda[4:12, best], 3)

MCMCdraws %>%
  filter(iter %in% sample(2501:10000, 1000)) %>%
  spread(var, draw) %>%
  select(-iter) %>%
  rbind(mvtnorm::rmvnorm(1000, vbmu, t(vbU) %*% vbU) %>%
          as.data.frame() %>%
          rename(alpha = V1, beta = V2, nu = V3) %>%
          mutate(method = 'VB-Reparam')) %>%
  rbind(particles %>% 
          as.data.frame() %>%
          rename(alpha = V1, beta = V2, nu = V3) %>%
          mutate(method = 'VB-Stein'))-> sample

ggplot(sample) + geom_density2d(aes(alpha, beta, colour = method)) + theme_bw()
ggplot(sample) + geom_density2d(aes(alpha, nu, colour = method)) + theme_bw()
ggplot(sample) + geom_density2d(aes(beta, nu, colour = method)) + theme_bw()


