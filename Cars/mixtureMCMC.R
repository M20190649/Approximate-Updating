library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('mixtureMCMCMetHast.cpp')

mixtureMCMC <- function(data, reps, draw, hyper, thin = 1, K = 2, error = 'gaussian'){
  N <- length(data)
  stepsize <- rep(0.01, N)
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
    saveDraws[[i+1]] <- list(theta = matrix(0, nSave, dim), k = rep(0, nSave))
  }
  # changing MH acceptance rate
  alpha <- - qnorm(0.234/2)
  stepsizeCons <- (1 - 1/dim) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)  +  1 / (dim * 0.234 * (1 - 0.234))
  
  for(i in 2:reps){
    # timing
    if(i == 2){
      startTime <- Sys.time()
    } else if(i < 5000 & i %% 100 == 2){
      timePerIter <- (Sys.time() - startTime)/(i - 2)
      class(timePerIter) <- 'numeric'
      if(attr(timePerIter, 'units') == 'mins'){
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
        p[k] <-  0.5 * log(det(draw[[1]][[k]]$varInv)) - 0.5 * (draw[[j+1]]$theta - draw[[1]][[k]]$mean) %*% 
          draw[[1]][[k]]$varInv %*% (draw[[j+1]]$theta - draw[[1]][[k]]$mean)
      }
      
      if(sum(exp(p)) == 0){
        draw[[j+1]]$k <- which.max(p)
      } else {
        prob <- exp(p) / sum(exp(p))
        draw[[j+1]]$k <- base::sample(1:K, 1, prob=prob)
      }
    }
    # thetaHat_k
    for(k in 1:K){
      ki <- 0
      for(j in 1:N){
        if(draw[[j+1]]$k == k){
          ki <- ki + 1
        }
      }
      var <- hyper[[k]]$varInv + ki * draw[[1]][[k]]$varInv
      var <- solve(var)
      mean <- var %*% hyper[[k]]$varInv %*% hyper[[k]]$mean
      
      for(j in 1:N){
        if(draw[[j+1]]$k == k){
          mean <- mean + var %*% draw[[1]][[k]]$varInv %*% draw[[j+1]]$theta
        }
      }
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
