library(tidyverse)

carsVB <- function(data, lambda, S = 25, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, 
                   dimTheta = 4, model = arimaDeriv, ...){
  dimLambda <- length(lambda)
  sobol <- sobol_points(100+S, dimTheta)
  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- numeric(dimLambda)
  V <- numeric(dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    #if(any(is.na(lambda))){
    #  break
    #}
    grad <- matrix(0, dimLambda, S)
    eval <- numeric(S)
    q <- numeric(S)
    unif <- shuffle(sobol)
    epsilon <- qnorm(unif[101:(100+S), ])
    epsilon[epsilon < -3] = -3
    epsilon[epsilon > 3] = 3
    if(S == 1){
      logpj <- model(data, lambda, epsilon, ...)
      eval <- logpj$val
      grad <- logpj$grad
      q <- sum(dnorm(epsilon, log=TRUE))
      gradient <- grad
      gradientSq <- grad^2
      LB[iter] <- eval - q
    } else {
      for(s in 1:S){
        logpj <- model(data, lambda, epsilon[s,], ...)    
        eval[s] <- logpj$val
        grad[,s] <- logpj$grad
        q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
      }
      eval[eval == -Inf] = NA
      gradient <- rowMeans(grad, na.rm = TRUE)
      gradientSq <- rowMeans(grad^2, na.rm=TRUE)
      LB[iter] <- mean(eval - q, na.rm=TRUE) 
    }
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    if(any(is.na(alpha * Mst / sqrt(Vst + e)))){
      print('Break')
      break
    }
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 100 == 0){
      print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', min(iter-1, maxIter), ' ELBO: ', LB[min(iter-1, maxIter)]))
  return(list(lambda=lambda, LB = LB[1:min(iter-1, maxIter)], iter = min(maxIter, iter-1)))
}

vbDensity <- function(fit, transform, names, supports = NULL){
  n = length(transform)
  if(is.null(supports)){
    supports = as.list(rep(NA, n))
  }
  mu = fit$mean
  if(length(fit$U) == n^2){
    u = matrix(fit$U, n)
    sigma = sqrt(diag(t(u) %*% u))
  } else {
    sigma = fit$U
  }
  dens = data.frame()
  for(i in 1:n){
    if(transform[i] == 'exp'){
      if(is.na(supports[[i]][1])){
        mean = exp(mu[i] + 0.5 * sigma[i]^2)
        stdev = sqrt((exp(sigma[i]^2) - 1)*exp(2*mu[i]+sigma[i]^2))
        support = seq(max(1e-08, mean - 5*stdev), mean+5*stdev, length.out=1000)
      } else {
        support = supports[[i]]
      }
      density = dlnorm(support, mu[i], sigma[i])
    } else if (transform[i] == 'sigmoid') {
      if(is.na(supports[[i]][1])){
        sample = 1 / (1 + exp(-rnorm(1000, mu[i], sigma[i])))
        mean = mean(sample)
        stdev = sd(sample)
        support = seq(max(0.001, mean-5*stdev), min(0.999, mean+5*stdev), length.out=1000)
      } else {
        support = supports[[i]]
      }
      density = dnorm(log(support / (1-support)), mu[i], sigma[i]) / (support - support^2)
    } else if(transform[i] == 'identity') {
      if(is.na(supports[[i]][1])){
        support = seq(mu[i] - 5*sigma[i], mu[i] + 5*sigma[i], length.out=1000)
      } else {
        support = supports[[i]]
      }
      density = dnorm(support, mu[i], sigma[i])
    } else if(transform[i] == 'stretchedSigmoid'){
      if(is.na(supports[[i]][1])){
        sample = 2 / (1 + exp(-rnorm(1000, mu[i], sigma[i]))) - 1
        mean = mean(sample)
        stdev = sd(sample)
        support = seq(max(-0.999, mean-5*stdev), min(0.999, mean+5*stdev), length.out=1000)
      } else {
        support = supports[[i]]
      }
      density = dnorm(-log(2/(support+1)-1), mu[i], sigma[i]) * 2 / (2*(support+1) - (support+1)^2)
    }
    df = data.frame(support = support, density = density, var = names[i])
    dens = rbind(dens, df)
  }
  dens
}

compareModels <- function(origList, IDvec, j, lags, transform, names, heir = TRUE, ...) {
  size = 2 + 2 * lags
  L3Y600 %>%
    filter(ID == IDvec[j]) %>%
    select(v, delta) -> car
  
  data <- cbind(car$v[2:nrow(car)] - car$v[1:(nrow(car)-1)], car$delta[2:nrow(car)] - pi/2)
  mu <- rep(0, 2 + 2 * lags)
  sd <- rep(0.2, 2 + 2 * lags)
  lambda <- matrix(c(mu, diag(sd)), ncol=1)
  hyper <- c(2, 0.0002, 2, 0.00002, rep(c(0, 1), 2 * lags))
  
  
  fit <- carsVB(data, lambda, hyper=hyper, dimTheta=size, model = arDeriv, lags =lags, ...)$lambda
  
  if(heir){
    density1 <- vbDensity(fit = origList[[j+1]],
                          transform = transform,
                          names = names)
    density1$method <- 'heirarchical'
  } else {
    density1 <- vbDensity(fit = origList,
                          transform = transform,
                          names = names)
    density1$method <- 'updater'
  }
  
  density2 <- vbDensity(fit = list(mean = fit[1:size], U = fit[(size+1):(size*(size+1))]),
                        transform = transform,
                        names = names,
                        supports = density1 %>% 
                          select(var, support) %>% 
                          group_by(var) %>%
                          mutate(n = 1:n()) %>% 
                          spread(var, support) %>% 
                          select(-n) %>%
                          as.list())
  density2$method <- 'single'
  
  density1 %>%
    rbind(density2) %>%
    ggplot() + geom_line(aes(support, density)) + 
    facet_wrap(method~var, scales = 'free', ncol = size) + 
    theme_bw() + 
    theme(strip.background = element_blank()) + 
    labs(x = NULL, y = NULL, title = paste0('carID: ', IDvec[j])) -> plot
  print(plot)
  
  density1 %>%
    rbind(density2) %>%
    group_by(var, method) %>%
    summarise(map = support[which.max(density)])
}

drawTheta <- function(thetaHat, var){
  require(mvtnorm)
  B <- t(matrix(c(thetaHat[5], 0, 0, 0, thetaHat[6:7], 0, 0, thetaHat[8:10], 0, thetaHat[11:14]), 4)) %*% 
    matrix(c(thetaHat[5], 0, 0, 0, thetaHat[6:7], 0, 0, thetaHat[8:10], 0, thetaHat[11:14]), 4)
  mean <- thetaHat[1:4]
  
  draw <- NULL
  for(i in 1:1000){
    inv <- rnorm(10, var$mean, abs(var$sd))
    Uinv <- matrix(c(Linv[1], 0, 0, 0, Linv[2:3], 0, 0, Linv[4:6], 0, Linv[7:10]), 4, 4)
    Sig <- solve(t(Uinv) %*% Uinv)
    Var <- B + Sig
    draw <- rbind(draw, rmvnorm(1, mean, Var))
  }
  draw
}

updateVB <- function(data, lambda, stepsize = 10, lags = 2, model = arUpdater,  ...){
  dim <- 2 + 2 * lags
  n <- stepsize + lags
  starting <- seq(1, nrow(data), stepsize)
  for(t in seq_along(starting)){
    mean <- lambda[1:dim]
    u <- matrix(lambda[(dim+1):length(lambda)], dim)
    linv <- solve(t(u))
    dat <- data[starting[t]:min(nrow(data), starting[t]+n-1),]
    if(is.matrix(dat)){
      if(nrow(dat) > lags){
        lambda <- carsVB(dat, lambda, dimTheta = dim, model = model, mean = mean, Linv = linv, lags = lags, ...)$lambda
      }
    } 
  }
  lambda
}

updateVBMix <- function(data, lambda, stepsize = 10, lags = 2, K = 6){
  n <- stepsize + lags
  starting <- seq(1, nrow(data), stepsize)
  for(t in seq_along(starting)){
    mean <- prevFit[[4]][1:(dim*K)]
    linv <- NULL
    dets <- NULL
    for(k in 1:K){
      u <- matrix(lambda[6*K + 1:36 + (k-1)*36], 6)
      linv <- rbind(linv, solve(t(u)))
      dets <- c(dets, det(solve(u)))
    }
    weights <- lambda[K * 42 + 1:K]
    dat <- data[starting[t]:min(nrow(data), starting[t]+n-1),]
    if(is.matrix(dat)){
      if(nrow(dat) > lags){
        lambda <- carsVBMixScore(data, prevFit[[4]], lags = 2, 
                                 priorMix = list(mean = mean, linv = linv, dets = dets, weights = weights),
                                 S = 5, maxIter = 10000, threshold = 0.1)$lambda
      }
    } 
  }
  lambda
}

sampleCars <- function(cars, IDvector){
  cars %>%
    filter(ID %in% IDvector) %>%
    select(v, delta, ID, class) -> carsSub
  
  carsSub %>%
    group_by(ID) %>%
    mutate(n = seq_along(v)) %>%
    filter(n > 1) %>%
    summarise(n =n()) %>%
    .$n %>%
    cumsum() -> obsSum
  
  carsSub%>%
    group_by(ID) %>%
    mutate(n = seq_along(v),
           vlag = ifelse(n == 1, 0, lag(v)),
           vdiff = v - vlag,
           delta = delta - pi/2) %>%
    filter(n > 1) %>%
    ungroup() %>%
    select(vdiff, delta) %>%
    as.matrix() -> data
  
  carsSub %>% 
    group_by(ID) %>%
    summarise(class = head(class, 1)) %>%
    .$class -> class
  return(list(data=data, obsSum=obsSum, class = class))
}

fitCarMods <- function(data, prevFit, increment, starting){
  S <- nrow(data)
  results <- list()
  if(increment){
    for(k in 1:3){
      results[[k]] <- updateVB(data, prevFit[[k]], stepsize = S, lags = 2, model = arUpdater)$lambda
    }
    #results[[4]] <- updateVBMix(data, prevFit[[4]], stepsize = S, lags = 2, model = arUpdateMix)$lambda
  } else {
    for(k in 1:3){
      mean <- prevFit[[k]][1:6]
      u <- matrix(prevFit[[k]][7:42], 6)
      linv <- solve(t(u))
      
      if(k == 1){
        results[[k]] <- carsVB(data, starting, lags = 2, model = arUpdater, mean = mean, Linv = linv, dimTheta = 6)$lambda
      } else {
        results[[k]] <- carsVB(data, prevFit[[k]], lags = 2, model = arUpdater, mean = mean, Linv = linv, dimTheta = 6)$lambda
      }
    }
    #mean <- prevFit[[4]][1:(6*6)]
    #linv <- NULL
    #logdets <- NULL
    #for(k in 1:6){
    #  u <- matrix(prevFit[[4]][36 + 1:36 + (k-1)*36], 6)
    #  linv <- rbind(linv, solve(t(u)))
    #  logdets <- c(logdets, -log(det(u)))
    #}
    #weights <- prevFit[[4]][36*7 + 1:6]
    #results[[4]] <- carsVB(data, prevFit[[4]], lags = 2, model = arUpdateMix, mean = mean, Linv = linv, logdets = logdets, weights = weights)
  }
  results
}

