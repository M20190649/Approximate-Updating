carsVB <- function(data, lambda, S, maxIter, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, threshold = 0.01, 
                   dimTheta = 6, dimLambda = NULL, model = arimaDeriv, ...){
  if(is.null(dimLambda)){
    dimLambda = dimTheta * (dimTheta + 1)
  }
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
    for(s in 1:S){
      if(S == 1){
        logpj <- model(data, lambda, epsilon, ...)
        eval <- logpj$val
        grad <- logpj$grad
        logpj <- model(data, lambda, epsilon, ...)     
        q <- sum(dnorm(epsilon[s,], log=TRUE))
      } else {
        logpj <- model(data, lambda, epsilon[s,], ...)    
        eval[s] <- logpj$val
        grad[,s] <- logpj$grad
        q[s] <- sum(dnorm(epsilon[s,], log=TRUE))
      }
     
    }
    gradient <- rowMeans(grad, na.rm = TRUE)
    gradientSq <- rowMeans(grad^2, na.rm=TRUE)
    LB[iter] <- mean(eval - q, na.rm=TRUE) 
    if(any(is.na(gradient)) | any(is.na(gradientSq))) {
      break
    }
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mst <- M / (1 - beta1^iter)
    Vst <- V / (1 - beta2^iter)
    lambda <- lambda + alpha * Mst / sqrt(Vst + e)
    if(iter %% 5 == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-4)])
      diff <- abs(meanLB - oldMeanLB)
    } 
    if(iter %% 25 == 0){
      print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  return(lambda)
}

vbDensity <- function(fit, transform, names, supports = NULL){
  n = length(transform)
  if(is.null(supports)){
    supports = as.list(rep(NA, n))
  }
  mu = fit$mean
  if(length(fit$U) == n^2){
    u = matrix(fit$U, n)
  } else {
    u = NULL
    for(i in 1:n){
      u = c(u, fit$U[sum(0:(i-1))+1:i], rep(0, n - i))
    }
    u = matrix(u, n)
  }
  sigma = sqrt(diag(t(u) %*% u))
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

compareModels <- function(heirList, IDvec, j, transform = c(rep('exp', 2), rep('stretchedSigmoid', 2)), 
                          names = c('sigma^2[V]', 'sigma^2[D]', 'phi[V]', 'phi[D]')) {
  size = length(transform)
  L3Y600 %>%
    filter(ID == IDvec[j]) %>%
    select(v, delta) -> car
  
  data <- cbind(car$v[2:nrow(car)] - car$v[1:(nrow(car)-1)], car$delta[2:nrow(car)])
  mu <- rep(0, size)
  sd <- rep(1, size)
  lambda <- matrix(c(mu, diag(sd)), nrow=size*(size+1))
  hyper <- c(2, 0.0002, 2, 0.00002, rep(1, 2*(size-2)))
  
  if(size == 4){
    fit <- carsVB(data, lambda, hyper=hyper, S=5, maxIter=5000, alpha=0.01, beta1=0.9, beta2=0.99,
                  dimTheta=4, model = ar1Deriv, threshold=0.01)
  } else {
    fit <- carsVB(data, lambda, hyper=hyper, S=5, maxIter=5000, alpha=0.01, beta1=0.9, beta2=0.99,
                  dimTheta=6, model = var1Deriv, threshold=0.01)
  }
  
  heirDensity <- vbDensity(fit = heirList[[j+1]],
                           transform = transform,
                           names = names)
  heirDensity$method <- 'heirarchical'
  
  density <- vbDensity(fit = list(mean = fit[1:size], U = fit[(size+1):(size*(size+1))]),
                       transform = transform,
                       names = names,
                       supports = heirDensity %>% 
                         select(var, support) %>% 
                         group_by(var) %>%
                         mutate(n = 1:n()) %>% 
                         spread(var, support) %>% 
                         select(-n) %>%
                         as.list())
  density$method <- 'single'

  density %>%
    rbind(heirDensity) %>%
    ggplot() + geom_line(aes(support, density)) + 
    facet_wrap(method~var, scales = 'free', ncol = size) + 
    theme_bw() + 
    theme(strip.background = element_blank()) + 
    labs(x = NULL, y = NULL, title = paste0('carID: ', IDvec[j])) -> plot
  print(plot)
  
  density %>%
    rbind(heirDensity) %>%
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
