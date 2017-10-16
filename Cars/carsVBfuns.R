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
  mu = fit[1:n]
  if(length(fit) == n*(n+1)){
    u = matrix(fit[(n+1):length(fit)], n)
  } else {
    u = matrix(c(fit[n+1], 0, 0, 0, fit[n+2:3], 0, 0, fit[n+4:6], 0, fit[n+7:10]), n)
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

compareModels <- function(heirList, IDvec, j) {
  L3Y600 %>%
    filter(ID == IDvec[j]) %>%
    select(v, delta) -> car
  
  data <- cbind(car$v[2:nrow(car)] - car$v[1:(nrow(car)-1)], car$delta[2:nrow(car)])
  mu <- c(0, 0, 0, 0)
  sd <- c(1, 1, 1, 1)
  lambda <- matrix(c(mu, diag(sd)), nrow=20)
  hyper <- c(2, 0.0002, 2, 0.00002, 1, 1, 1, 1)
  
  fit <- carsVB(data, lambda, hyper=hyper, S=5, maxIter=5000, alpha=0.01, beta1=0.9, beta2=0.99,
                dimTheta=4, model = ar1Deriv, threshold=0.01)
  
  
  heirDensity <- vbDensity(fit = heirList[[j+1]],
                           transform = c(rep('exp', 2), rep('stretchedSigmoid', 2)),
                           names = c('sigma^2[V]', 'sigma^2[D]', 'phi[V]', 'phi[D]'))
  heirDensity$method <- 'heirarchical'
  
  density <- vbDensity(fit = fit,
                       transform = c(rep('exp', 2), rep('stretchedSigmoid', 2)),
                       names = c('sigma^2[V]', 'sigma^2[D]', 'phi[V]', 'phi[D]'),
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
    facet_wrap(method~var, scales = 'free', ncol =4) + 
    theme_bw() + 
    theme(strip.background = element_blank()) + 
    labs(x = NULL, y = NULL, title = paste0('carID: ', IDvec[j])) -> plot
  print(plot)
  
  density %>%
    rbind(heirDensity) %>%
    group_by(var, method) %>%
    summarise(map = support[which.max(density)])
}
