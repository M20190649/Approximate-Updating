library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('mixtureMCMCMetHast.cpp')

mixtureMCMC <- function(data, reps, draw, hyper, thin = 1, K = 2, error = 'gaussian', stepsize = 0.01){
  N <- length(data)
  stepsize <- rep(stepsize, N)
  accept <- rep(0, N)
 
  
  #set up likelihood function and theta dimension
  if(error == 'gaussian'){
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  #set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- list(list())
  for(i in 1:K){
    saveDraws[[1]][[i]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
  }
  for(i in 1:N){
    saveDraws[[i+1]] <- list(theta = matrix(0, nSave, dim), k = rep(0, nSave), pi = matrix(0, nSave, K))
  }
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)  +  1 / (dim * 0.234 * (1 - 0.234))
  
  for(i in 2:reps){
    # timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # theta_i
    for(j in 1:N){
      candidate <- draw[[j+1]]$theta +  stepsize[j] * rnorm(dim)
      group <- draw[[j+1]]$k
      canDens <- likelihood(data[[j]], candidate, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv)
      oldDens <- likelihood(data[[j]], draw[[j+1]]$theta, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv)
   
      ratio <- exp(canDens - oldDens)
    
      c <- stepsize[j] * stepsizeCons
      if(runif(1) < ratio){
        accept[j] <- accept[j] + 1
        draw[[j+1]]$theta <- candidate
        stepsize[j] <- stepsize[j] + c * (1 - 0.234) / (28 + i)
      } else {
        stepsize[j] <- stepsize[j] - c * 0.234 / (28 + i)
      }
    }
    # k_i
    for(j in 1:N){
      p <- numeric(K)
      for(k in 1:K){
        p[k] <-  log(draw[[j+1]]$pi[k])  +  0.5 * log(det(draw[[1]][[k]]$varInv)) -
          0.5 * (draw[[j+1]]$theta - draw[[1]][[k]]$mean) %*% draw[[1]][[k]]$varInv %*%
          (draw[[j+1]]$theta - draw[[1]][[k]]$mean)
      }
      p <- p - max(p) 
      p <- exp(p) / sum(exp(p))
      draw[[j+1]]$k <- base::sample(1:K, 1, prob=p)
    }
    #pi_i
    for(j in 1:N){
      group <- rep(0, K)
      group[draw[[j+1]]$k] <- 1
      draw[[j+1]]$pi <- c(MCMCpack::rdirichlet(1, hyper$alpha + group))
    }
    # thetaHat_k
    sumK <- rep(0, K)
    sumTheta <- matrix(0, dim, K)
    for(j in 1:N){
      sumK[draw[[j+1]]$k] <- sumK[draw[[j+1]]$k] + 1
      sumTheta[,draw[[j+1]]$k] <-  sumTheta[,draw[[j+1]]$k] + draw[[j+1]]$theta
    }
    for(k in 1:K){
      var <- hyper[[k]]$varInv + sumK[k] * draw[[1]][[k]]$varInv
      var <- solve(var)
      mean <- var %*% hyper[[k]]$varInv %*% hyper[[k]]$mean  +  var %*% draw[[1]][[k]]$varInv %*%  sumTheta[,k]
    
      draw[[1]][[k]]$mean <- c(rmvnorm(1, mean, var))
      vardf <- hyper[[k]]$v 
      scaleMat <- hyper[[k]]$scale
      for(j in 1:N){
        if(draw[[j+1]]$k == k){
          vardf <- vardf + 1
          scaleMat <- scaleMat + outer(draw[[j+1]]$theta - draw[[1]][[k]]$mean, draw[[j+1]]$theta - draw[[1]][[k]]$mean)
        }
      }
      draw[[1]][[k]]$varInv <- rWishart(1, vardf, solve(scaleMat))[,,1]
    }
    # save draws
    if(i %% thin == 0){
      for(k in 1:K){
        saveDraws[[1]][[k]]$mean[i/thin,] <- draw[[1]][[k]]$mean
        saveDraws[[1]][[k]]$varInv[,,i/thin] <- draw[[1]][[k]]$varInv
      }
      for(j in 1:N){
        saveDraws[[j+1]]$theta[i/thin, ] <- draw[[j+1]]$theta
        saveDraws[[j+1]]$k[i/thin] <- draw[[j+1]]$k
        saveDraws[[j+1]]$pi[i/thin, ] <- draw[[j+1]]$pi
      }
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}

sliceSampler <- function(data, reps, draw, hyper, thin = 1, K = 20, error = 'gaussian', stepsize = 0.01){
  N <- length(data)
  stepsize <- rep(stepsize, N)
  accept <- rep(0, N)
  
  #set up likelihood function and theta dimension
  if(error == 'gaussian'){
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  #set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- list(list())
  for(i in 1:K){
    saveDraws[[1]][[i]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
  }
  saveDraws[[1]]$pi <- matrix(0, nSave, K)
  for(i in 1:N){
    saveDraws[[i+1]] <- list(theta = matrix(0, nSave, dim), k = rep(0, nSave))
  }
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2*3.141598) * exp(alpha^2/2) / (2 * alpha)  +  1 / (dim * 0.234 * (1 - 0.234))
  for(i in 2:reps){
    # timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }

    # thetaHat_k
    sumK <- rep(0, K)
    sumTheta <- matrix(0, dim, K)
    for(j in 1:N){
      sumK[draw[[j+1]]$k] <- sumK[draw[[j+1]]$k] + 1
      sumTheta[,draw[[j+1]]$k] <-  sumTheta[,draw[[j+1]]$k] + draw[[j+1]]$theta
    }
    for(k in 1:K){
      if(sumK[k] > 0){
        var <- hyper[[k]]$varInv + sumK[k] * draw[[1]][[k]]$varInv
        var <- solve(var)
        mean <- var %*% hyper[[k]]$varInv %*% hyper[[k]]$mean  +  var %*% draw[[1]][[k]]$varInv %*%  sumTheta[,k]
        
        draw[[1]][[k]]$mean <- c(mvtnorm::rmvnorm(1, mean, var))
        vardf <- hyper[[k]]$v 
        scaleMat <- hyper[[k]]$scale
        for(j in 1:N){
          if(draw[[j+1]]$k == k){
            vardf <- vardf + 1
            scaleMat <- scaleMat + outer(draw[[j+1]]$theta - draw[[1]][[k]]$mean, draw[[j+1]]$theta - draw[[1]][[k]]$mean)
          }
        }
        draw[[1]][[k]]$varInv <- rWishart(1, vardf, solve(scaleMat))[,,1]
      } else {
        draw[[1]][[k]]$mean <- c(mvtnorm::rmvnorm(1, hyper[[k]]$mean, solve(hyper[[k]]$varInv)))
        draw[[1]][[k]]$varInv <- rWishart(1, hyper[[k]]$v, solve(hyper[[k]]$scale))[,,1]
      }
    }
    # v
    v <- rep(0, K)
    for(k in 1:(K-1)){
      v[k] <- rbeta(1, 1 + sumK[k],  hyper$M  + sum(sumK[(k+1):K]))
    }
    v[K] <- 1
    # pi
    draw[[1]]$pi <- v * cumprod(c(1, 1-v[1:(K-1)]))
    # u
    #u <- runif(N)
    #for(j in 1:N){
    # u[j] <- u[j] * draw[[1]]$pi[draw[[j+1]]$k] 
    #}
    # k_i
    for(j in 1:N){
      p <- rep(0, K)
      for(k in 1:K){
        p[k] <-  log(draw[[1]]$pi[k])  +  0.5 * log(det(draw[[1]][[k]]$varInv)) -
          0.5 * (draw[[j+1]]$theta - draw[[1]][[k]]$mean) %*% draw[[1]][[k]]$varInv %*%
          (draw[[j+1]]$theta - draw[[1]][[k]]$mean)
      }
      p <- p - max(p) 
      p <- exp(p) / sum(exp(p))
      draw[[j+1]]$k <- base::sample(1:K, 1, prob=p)
    }
    
    # theta_i
    for(j in 1:N){
      candidate <- draw[[j+1]]$theta +  stepsize[j] * rnorm(dim)
      group <- draw[[j+1]]$k
      canDens <- likelihood(data[[j]], candidate, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv)
      oldDens <- likelihood(data[[j]], draw[[j+1]]$theta, draw[[1]][[group]]$mean, draw[[1]][[group]]$varInv)
      
      ratio <- exp(canDens - oldDens)
      
      c <- stepsize[j] * stepsizeCons
      if(runif(1) < ratio){
        accept[j] <- accept[j] + 1
        draw[[j+1]]$theta <- candidate
        stepsize[j] <- stepsize[j] + c * (1 - 0.234) / (28 + i)
      } else {
        stepsize[j] <- stepsize[j] - c * 0.234 / (28 + i)
      }
    }
    
  
    # save draws
    if(i %% thin == 0){
      for(k in 1:K){
        saveDraws[[1]][[k]]$mean[i/thin,] <- draw[[1]][[k]]$mean
        saveDraws[[1]][[k]]$varInv[,,i/thin] <- draw[[1]][[k]]$varInv
      }
      for(j in 1:N){
        saveDraws[[j+1]]$theta[i/thin, ] <- draw[[j+1]]$theta
        saveDraws[[j+1]]$k[i/thin] <- draw[[j+1]]$k
      }
      saveDraws[[1]]$pi[i/thin,] <- draw[[1]]$pi
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}

hierNoMixMCMC <-function(data, reps, draw, hyper, thin = 1, error = 'gaussian', stepsize = 0.01){
  N <- length(data)
  stepsize <- rep(stepsize, N)
  accept <- rep(0, N)
  
  #set up likelihood function and theta dimension
  if(error == 'gaussian'){
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  #set up storage for saved draws
  nSave <- floor(reps / thin)
  saveDraws <- list()
  saveDraws[[1]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
  for(i in 1:N){
    saveDraws[[i+1]] <- matrix(0, nSave, dim)
  }
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)  +  1 / (dim * 0.234 * (1 - 0.234))
  
  for(i in 2:reps){
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # theta_i
    for(j in 1:N){
      candidate <- draw[[j+1]]  +  stepsize[j] * rnorm(dim)
      canDens <- likelihood(data[[j]], candidate, draw[[1]]$mean, draw[[1]]$varInv)
      oldDens <- likelihood(data[[j]], draw[[j+1]], draw[[1]]$mean, draw[[1]]$varInv)
      
      ratio <- exp(canDens - oldDens)
      
      c <- stepsize[j] * stepsizeCons
      if(runif(1) < ratio){
        accept[j] <- accept[j] + 1
        draw[[j+1]] <- candidate
        stepsize[j] <- stepsize[j] + c * (1 - 0.234) / (28 + i)
      } else {
        stepsize[j] <- stepsize[j] - c * 0.234 / (28 + i)
      }
    }

    sumTheta <- rep(0, dim)
    for(j in 1:N){
      sumTheta <-  sumTheta + draw[[j+1]]
    }
    var <- hyper$varInv + N * draw[[1]]$varInv
    var <- solve(var)
    mean <- var %*% hyper$varInv %*% hyper$mean  +  var %*% draw[[1]]$varInv %*%  sumTheta
      
    draw[[1]]$mean <- c(rmvnorm(1, mean, var))
    vardf <- hyper$v 
    scaleMat <- hyper$scale
    for(j in 1:N){
      vardf <- vardf + 1
      scaleMat <- scaleMat + outer(draw[[j+1]] - draw[[1]]$mean, draw[[j+1]] - draw[[1]]$mean)
    }
    draw[[1]]$varInv <- rWishart(1, vardf, solve(scaleMat))[,,1]
    
    # save draws
    if(i %% thin == 0){
      saveDraws[[1]]$mean[i/thin,] <- draw[[1]]$mean
      saveDraws[[1]]$varInv[,,i/thin] <- draw[[1]]$varInv
      for(j in 1:N){
        saveDraws[[j+1]][i/thin, ] <- draw[[j+1]]
      }
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
}
  
noHierMCMC <- function(data, reps, draw, hyper, thin = 1, error = 'gaussian', stepsize = 0.01){
  N <- length(data)
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
  
  # sums of data
  a11 <- a22 <- a12 <- a01 <- a02 <- 0
  d11 <- d22 <- d12 <- d01 <- d02 <- 0
  for(i in 1:N){
    T <- nrow(data[[i]])
    for(t in 3:T){
      a11 <- a11 + data[[i]][t-1, 1]^2
      a22 <- a22 + data[[i]][t-2, 1]^2
      a12 <- a12 + data[[i]][t-1, 1] * data[[i]][t-2, 1]
      a01 <- a01 + data[[i]][t-1, 1] * data[[i]][t, 1]
      a02 <- a02 + data[[i]][t-2, 1] * data[[i]][t, 1]
      d11 <- d11 + data[[i]][t-1, 2]^2
      d22 <- d22 + data[[i]][t-2, 2]^2
      d12 <- d12 + data[[i]][t-1, 2] * data[[i]][t-2, 2]
      d01 <- d01 + data[[i]][t-1, 2] * data[[i]][t, 2]
      d02 <- d02 + data[[i]][t-2, 2] * data[[i]][t, 2]
    }
  }
  # Metropolis Hastings variables
  MHvar <- 1:2
  if(error == 't'){
    MHvar <- c(MHvar, 7:8)
  }
  
  for(i in 2:reps){
    # timing
    if(i == 50){
      startTime <- Sys.time()
    } else if(i == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    }
    # log_sigma_squared and nu (if error ~ t)
    candidate <- draw
    candidate[MHvar] <- candidate[MHvar] + stepsize * rnorm(dim-4)
    canDens <- 0
    oldDens <- 0
    for(j in 1:N){
      canDens <- canDens + likelihood(data[[j]], candidate, hyper$mean, hyper$varInv)
      oldDens <- oldDens + likelihood(data[[j]], draw, hyper$mean, hyper$varInv)
    }
    ratio <- exp(canDens - oldDens)
    c <- stepsize / stepsizeCons
    if(runif(1) < ratio){
      accept <- accept + 1
      draw <- candidate
      stepsize <- stepsize + c * (1 - 0.44) / (18 + i)
    } else {
      stepsize <- stepsize - c * 0.44 / (18 + i)
    }
    # phi 1
    meanNumer <- hyper$var[3] * (a01 - draw[4] * a12)  +  exp(draw[1]) * hyper$mean[3]
    meanDenom <- hyper$var[3] * a11  +  exp(draw[1])
    var <- exp(draw[1]) * hyper$var[3] / meanDenom
    draw[3] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # phi 2
    meanNumer <- hyper$var[4] * (a02 - draw[3] * a12)  +  exp(draw[1]) * hyper$mean[4]
    meanDenom <- hyper$var[4] * a22  +  exp(draw[1])
    var <- exp(draw[1]) * hyper$var[4] / meanDenom
    draw[4] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # gamma 1
    meanNumer <- hyper$var[5] * (d01 - draw[6] * d12)  +  exp(draw[2]) * hyper$mean[5]
    meanDenom <- hyper$var[5] * d11  +  exp(draw[2])
    var <- exp(draw[2]) * hyper$var[5] / meanDenom
    draw[5] <- rnorm(1, meanNumer / meanDenom, sqrt(var))
    
    # gamma 2
    meanNumer <- hyper$var[6] * (d02 - draw[5] * d12)  +  exp(draw[2]) * hyper$mean[6]
    meanDenom <- hyper$var[6] * d22  +  exp(draw[2])
    var <- exp(draw[2]) * hyper$var[6] / meanDenom
    draw[6] <- rnorm(1, meanNumer / meanDenom, sqrt(var))

    # save draws
    if(i %% thin == 0){
      saveDraws[i/thin,] <- draw
    }
    # print progress
    if(i %% 1000 == 0){
      mins <-  (reps - i) * timePerIter[1] / 60
      if(mins > 180){
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins / 60, 2), ' hours.'))
      } else {
        print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round(mins, 2), ' minutes.'))
      }
    }
  }
  list(draws = saveDraws, accept = accept, steps = stepsize)
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

NoGaps <- function(data, reps, draw, hyper, thin = 1, k = 10, stepsize = 0.01){
  N <- length(data)
  stepsize <- rep(stepsize, N)
  accept <- rep(0, N)

  #set up storage for saved draws
  nSave <- floor(reps / thin)
  save <- list(list())
  dim <- 6
  for(i in 1:k){
    save[[1]][[i]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
  }
  for(i in 1:N){
    save[[i+1]] <- list(theta = matrix(0, nSave, dim), s = rep(0, nSave))
  }
  
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/6) * sqrt(2*3.141598) * exp(alpha^2/2) / (2 * alpha)  +  1 / (6 * 0.234 * (1 - 0.234))
  
  # initialise groups randomly
  s <- sample(1:k, N, replace = TRUE)
  n <- table(s)
  
  for(iter in 2:reps){
    # timing
    if(iter == 50){
      startTime <- Sys.time()
    } else if(iter == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    } else if(iter %% 100 == 50){
      timePerIter <- (Sys.time() - startTime) / (iter - 50);
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      } else if(attr(timePerIter, 'units') == 'hours'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 3600
      }
    }
    
    # Resample s indices for each i
    # First sample a k+1'th mu / precision for any new groups
    draw[[1]][[k+1]] <- list(mean = c(mvtnorm::rmvnorm(1, hyper$mean, solve(hyper$varInv))),
                             varInv = rWishart(1, hyper$df, hyper$scale)[,,1])
    save[[1]][[k+1]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
    

    for(j in 1:N){
      skip = FALSE
      # mu_j is the only element of the cluster
      if(n[s[j]] == 1){
        prob <- k / (k+1)
        if(runif(1) > prob){
          # Keep s_j the same with above probability
          # Otherwise remove group, shuffle labels
          if(s[j] < k){
            # If the group already was k, then we keep it the same,
            # Otherwise shift everything above it down by one and put the now empty cluster in the k'th spot
            temp <- draw[[1]][[s[j]]]
            tempSave <- save[[1]][[s[j]]]
            
            draw[[1]][s[j]:(k-1)] <- draw[[1]][(s[j]+1):k]
            save[[1]][s[j]:(k-1)] <- save[[1]][(s[j]+1):k]
            
            draw[[1]][[k]] <- temp
            save[[1]][[k]] <- tempSave
            
            s[s > s[j]] <- s[s > s[j]] - 1
            s[j] <- k
            n <- table(s)
          }
          # remove the group
          k <- k - 1
        } else {
          # repeating the one element cluster
          skip = TRUE
        }
      }
      # If we didn't chose to keep a one element cluster the same we must resample s via MacEachern Muller 3.2
      if(!skip){
        n[s[j]] <- n[s[j]]  - 1
        probS <- rep(0, k+1)
        for(i in 1:k){
          probS[i] <- log(n[i])  +  nlogDensity(data[[j]], draw[[j+1]]$theta, draw[[1]][[i]]$mean, draw[[1]][[i]]$varInv)
        }
        probS[k+1] <- log(hyper$alpha / (k + 1))  +   nlogDensity(data[[j]], draw[[j+1]]$theta, draw[[1]][[k+1]]$mean, draw[[1]][[k+1]]$varInv)
        probS <- probS - max(probS)
        probS <- exp(probS) / sum(exp(probS))
        s[j] <- sample(1:(k+1), 1, prob = probS)
        if(s[j] == k+1){
          # add a new group, draw a new k+1'th group
          n <- c(n, 1)
          k <- k + 1
          draw[[1]][[k+1]] <- list(mean = c(mvtnorm::rmvnorm(1, hyper$mean, solve(hyper$varInv))),
                                   varInv = rWishart(1, hyper$df, hyper$scale)[,,1])
          save[[1]][[k+1]] <- list(mean = matrix(0, nSave, dim), varInv = array(0, dim = c(dim, dim, nSave)))
          
        } else {
          # add to the count of thew newly assigned group
          n[s[j]] <- n[s[j]] + 1
        }
      }
      n <- table(s)
    }
    draw[[1]] <- draw[[1]][1:k]
    save[[1]] <- save[[1]][1:k]
    
    # draw new values for each mu / varInv pair
    sumTheta <- matrix(0, 6, k)
    for(j in 1:N){
      sumTheta[,s[j]] <-  sumTheta[,s[j]] + draw[[j+1]]$theta
    }
    for(group in 1:k){
      var <- hyper$varInv + n[group] * draw[[1]][[group]]$varInv
      var <- solve(var)
      mean <- var %*% hyper$varInv %*% hyper$mean  +  var %*%  draw[[1]][[group]]$varInv %*%  sumTheta[,group]
      
      draw[[1]][[group]]$mean <- c(mvtnorm::rmvnorm(1, mean, var))
      vardf <- hyper$df 
      scaleMat <- hyper$scale
      for(j in 1:N){
        if(s[j] == group){
          vardf <- vardf + 1
          scaleMat <- scaleMat + outer(draw[[j+1]]$theta - draw[[1]][[group]]$mean, draw[[j+1]]$theta - draw[[1]][[group]]$mean)
        }
      }
      draw[[1]][[group]]$varInv <- rWishart(1, vardf, solve(scaleMat))[,,1]
    }
    
    # draw new values for each theta_i
    for(j in 1:N){
      candidate <- draw[[j+1]]$theta +  stepsize[j] * rnorm(dim)
      canDens <- nlogDensity(data[[j]], candidate, draw[[1]][[s[j]]]$mean, draw[[1]][[s[j]]]$varInv)
      oldDens <- nlogDensity(data[[j]], draw[[j+1]]$theta, draw[[1]][[s[j]]]$mean, draw[[1]][[s[j]]]$varInv)
      
      ratio <- exp(canDens - oldDens)
      
      c <- stepsize[j] * stepsizeCons
      if(runif(1) < ratio){
        accept[j] <- accept[j] + 1
        draw[[j+1]]$theta <- candidate
        stepsize[j] <- stepsize[j] + c * (1 - 0.234) / (28 + i)
      } else {
        stepsize[j] <- stepsize[j] - c * 0.234 / (28 + i)
      }
    }
    
    # save draws
    if(iter %% thin == 0){
      for(group in 1:k){
        save[[1]][[group]]$mean[iter/thin,] <- draw[[1]][[group]]$mean
        save[[1]][[group]]$varInv[,,iter/thin] <- draw[[1]][[group]]$varInv
      }
      for(j in 1:N){
        save[[j+1]]$theta[iter/thin, ] <- draw[[j+1]]$theta
        save[[j+1]]$s[iter/thin] <- s[j]
      }
    }
    # print progress
    if(iter %% 1000 == 0){
      time <- as.numeric((reps - iter) * timePerIter[1])
      hours <- floor(time / 3600)
      mins <- floor(time %% 3600 / 60)
      secs <- floor(time %% 60)
      
      print(paste0('Iteration: ', iter, '/', reps, '. Time Remaining: ', hours, ':', mins, ':', secs))
      
    }
  }
  
  save
}
  
NoGaps2 <- function(data, reps, draw, hyper, thin = 1, k = 10, stepsizeStart = 0.01){
  N <- length(data)
  stepsize <- rep(stepsizeStart, k)
  accept <- rep(0, k)
  
  #set up storage for saved draws
  nSave <- floor(reps / thin)
  saveS <- matrix(0, nSave, N)
  saveT <- list()
  for(i in 1:k){
    saveT[[i]] <- matrix(0, nSave, 6)
  }
  
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/6) * sqrt(2*3.141598) * exp(alpha^2/2) / (2 * alpha)  +  1 / (6 * 0.234 * (1 - 0.234))
  
  # initialise groups randomly
  s <- sample(1:k, N, replace = TRUE)
  n <- table(s)
  
  for(iter in 2:reps){
    # timing
    if(iter == 50){
      startTime <- Sys.time()
    } else if(iter == 150){
      timePerIter <- (Sys.time() - startTime) / 100
      class(timePerIter) <- 'numeric'
      print(paste0('Estimated Finishing Time: ', Sys.time() + timePerIter * (reps - 150)))
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      }
    } else if(iter %% 100 == 50){
      timePerIter <- (Sys.time() - startTime) / (iter - 50);
      if(attr(timePerIter, 'units') == 'mins'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 60
      } else if(attr(timePerIter, 'units') == 'hours'){
        attr(timePerIter, 'units') = 'secs'
        timePerIter <- timePerIter * 3600
      }
    }
    
    # Resample s indices for each i
    # First sample a k+1'th mu / precision for any new groups
    draw[[k+1]] <-  c(mvtnorm::rmvnorm(1, hyper$mean, hyper$var))
    saveT[[k+1]] <- matrix(0, nSave, 6)
    accept <- c(accept, 0)
    stepsize <- c(stepsize, stepsizeStart)
    for(j in 1:N){
      skip = FALSE
      # theta_j is the only element of the cluster
      if(n[s[j]] == 1){
        prob <- k / (k+1)
        u <- runif(1)
        if(u > prob){
          # Keep s_j the same with above probability
          # Otherwise remove group, shuffle labels
          if(s[j] < k){
            # If the group already was k, then we keep it the same,
            # Otherwise shift everything above it down by one and put the now empty cluster in the k'th spot
            temp <- draw[[s[j]]]
            tempSave <- saveT[[s[j]]]
            ss <- stepsize[s[j]]
            acc <- accept[s[j]]
            
            draw[s[j]:(k-1)] <- draw[(s[j]+1):k]
            saveT[s[j]:(k-1)] <- saveT[(s[j]+1):k]
            stepsize[s[j]:(k-1)] <- stepsize[(s[j]+1):k]
            accept[s[j]:(k-1)] <- accept[(s[j]+1):k]
            
            draw[[k]] <- temp
            saveT[[k]] <- tempSave
            stepsize[k] <- ss
            accept[k] <- acc
            
            s[s > s[j]] <- s[s > s[j]] - 1
            s[j] <- k
          }
          # remove the group
          k <- k - 1
          n <- table(s)
        } else {
          # repeating the one element cluster
          skip = TRUE
        }
      }
      # If we didn't chose to keep a one element cluster the same we must resample s via MacEachern Muller 3.2
      if(!skip){
        n[s[j]] <- n[s[j]]  - 1
        probS <- rep(0, k+1)
        for(i in 1:k){
          probS[i] <- log(n[i]) + nlikelihood(data[[j]], draw[[i]])
        }
        probS[k+1] <- log(hyper$alpha / (k + 1)) + nlikelihood(data[[j]], draw[[k+1]])
        probS <- probS - max(probS)
        probS <- exp(probS) / sum(exp(probS))
        s[j] <- sample(1:(k+1), 1, prob = probS)
        if(s[j] == k+1){
          # add a new group, draw a new k+1'th group
          n <- c(n[1:k], 1)
          k <- k + 1
          draw[[k+1]] <- c(mvtnorm::rmvnorm(1, hyper$mean, hyper$var))
          saveT[[k+1]] <- matrix(0, nSave, dim)
          accept <- c(accept, 0)
          stepsize <- c(stepsize, stepsizeStart)
        } else {
          # add to the count of thew newly assigned group
          n[s[j]] <- n[s[j]] + 1
        }
      }
    }
<<<<<<< HEAD
    n <- n[1:k]
    saveT <- saveT[1:k]
=======
    n <- table(s)
>>>>>>> c5acdb096e8bb10943437842d30d097d91be66ac
    accept <- accept[1:k]
    stepsize <- stepsize[1:k]
    draw <- draw[1:k]

    # draw new values for each theta_i
    MH <- NoGapsMH(data, draw, stepsize, accept, stepsizeCons, hyper$mean, hyper$varInv, s, iter)
    draw <- MH$draws
    stepsize <- MH$stepsize
    accept <- MH$accept
    
    # save draws
    if(iter %% thin == 0){
      for(group in 1:k){
        saveT[[group]][iter/thin,] <- draw[[group]]
      }
      saveS[iter/thin,] <- s
    }
    # print progress
    if(iter %% 1000 == 0){
      time <- as.numeric((reps - iter) * timePerIter[1])
      hours <- floor(time / 3600)
      mins <- floor(time %% 3600 / 60)
      secs <- floor(time %% 60)

      print(paste0('Iteration: ', iter, '/', reps, '. Time Remaining: ', hours, ':', mins, ':', secs))
      
    }
  }
  
  list(saveT, saveS, accept)
}  
  