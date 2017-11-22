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
    iter <- iter + 1
  }
  return(list(lambda=lambda, LB = LB[1:min(iter-1, maxIter)], iter = min(maxIter, iter-1)))
}

fitCarMods <- function(data, prevFit, increment, starting){
  S <- nrow(data)
  results <- list()
  #if(increment){
  #  for(k in 1:3){
  #    results[[k]] <- updateVB(data, prevFit[[k]], stepsize = S, lags = 2, model = arUpdater)$lambda
  #  }
    #results[[4]] <- updateVBMix(data, prevFit[[4]], stepsize = S, lags = 2, model = arUpdateMix)$lambda
  #} else {
    for(k in 1:3){
      mean <- prevFit[[k]][1:6]
      u <- matrix(prevFit[[k]][7:42], 6)
      linv <- solve(t(u))
      results[[k]] <- carsVB(data, starting, lags = 2, model = arUpdater, mean = mean, Linv = linv, dimTheta = 6)$lambda
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
  #}
  results
}

singleMCMCallMH <- function(data, reps, draw, hyper, thin = 1, error = 'gaussian', stepsize = 0.01){
  # set up likelihood function and theta dimension
  if(error == 'gaussian'){
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  accept <- 0
  # set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- matrix(0, nSave, dim)
  # changing MH acceptance rate
  stepsizeCons <- 0.44 * (1 - 0.44)
  
  for(i in 2:reps){
    candidate <- draw
    candidate <- candidate + stepsize * rnorm(dim)
    canDens <- likelihood(data, candidate, hyper$mean, hyper$varInv)
    oldDens <- likelihood(data, draw, hyper$mean, hyper$varInv)
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
    } else {
      stepsize <- stepsize - c * 0.44 / (18 + i)
    }
    
    # save draws
    if(i %% thin == 0){
      saveDraws[i/thin,] <- draw
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}