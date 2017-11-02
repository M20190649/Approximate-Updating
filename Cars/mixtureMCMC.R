library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('mixtureMCMCMetHast.cpp')

mixtureMCMC <- function(data, reps, K = 2, error = 'gaussian'){
  N <- length(data)
  stepsize <- rep(0.01, N)
  accept <- rep(0, N)

  if(error == 'gaussian'){
    draws <- list(list(list(mean = matrix(c(-6, -6, 0.5, 0, 0.3, 0.1), reps, 6, byrow = TRUE), varInv = array(0, dim = c(6, 6, reps))),
                       list(mean = matrix(c(-4, -4, -0.5, -0.2, 0, 0), reps, 6, byrow = TRUE), varInv = array(0, dim = c(6, 6, reps)))))
    hyper <- list()
    for(k in 1:K){
      diag(draws[[1]][[k]]$varInv[,,1]) = 10
      hyper[[k]] <- list(mean = c(-5, -5, rep(0, 4)), varInv = solve(diag(5, 6)),  v = 6, scale = diag(0.5, 6))
    }
    
    diag(draws[[1]][[2]]$varInv[,,1]) = 10
    diag(draws[[1]][[3]]$varInv[,,1]) = 10
    for(i in 1:N){
      draws[[i+1]] <- list(theta = matrix(c(-5, -5, 0, -0.1, 0.15, 0.05), reps, 6, byrow = TRUE), k = sample(1:K, reps, replace=TRUE))
    }
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    hyper <- list(list(mean = c(-8, -8, rep(0, 6)), varInv = solve(diag(5, 8)),  v = 8, scale = diag(0.5, 8)),
                  list(mean = c(-8, -8, rep(0, 6)), varInv = solve(diag(5, 8)),  v = 8, scale = diag(0.5, 8)))
    
    draws <- list(list(list(mean = matrix(c(-6, -5, 0.2, 0.1, 0.2, 0.1, 2, 2), reps, 8, byrow = TRUE), varInv = array(0, dim = c(8, 8, reps))),
                       list(mean = matrix(c(-4, -2, -0.5, -0.5, 0, 0, 2, 2), reps, 8, byrow = TRUE), varInv = array(0, dim = c(8, 8, reps)))))
    diag(draws[[1]][[1]]$varInv[,,1]) = 10
    diag(draws[[1]][[2]]$varInv[,,1]) = 10
    for(i in 1:N){
      draws[[i+1]] <- list(theta = matrix(c(-5, -3.5, 0.15, -0.2, 0.1, 0.05, 2, 2), reps, 8, byrow = TRUE), k = sample(1:2, reps, replace=TRUE))
    }
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  
  for(i in 2:reps){
    # theta_i
    if(i == 2){
      startTime <- Sys.time()
    } else if(i == 102){
      timePerIter <- (Sys.time() - startTime)/100
      class(timePerIter) <- 'numeric'
    }
    for(j in 1:N){
      candidate <- draws[[j+1]]$theta[i-1, ] +  stepsize[j] * rnorm(dim)
      group <- draws[[j+1]]$k[i-1]
      canDens <- likelihood(data[[j]], candidate, draws[[1]][[group]]$mean[i-1,], draws[[1]][[group]]$varInv[,,i-1])
      oldDens <- likelihood(data[[j]], draws[[j+1]]$theta[i-1,], draws[[1]][[group]]$mean[i-1,], draws[[1]][[group]]$varInv[,,i-1])
   
      ratio <- exp(canDens - oldDens)
      alpha <- - qnorm(0.234/2)
      c <- stepsize[j] * ((1 - 1/dim) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)
                          + 1 / (dim * 0.234 * (1 - 0.234)))
      if(runif(1) < ratio){
        accept[j] <- accept[j] + 1
        draws[[j+1]]$theta[i,] <- candidate
        stepsize[j] <- stepsize[j] + c * (1 - 0.234) / (28 + i)
      } else {
        draws[[j+1]]$theta[i, ] <- draws[[j+1]]$theta[i-1,]
        stepsize[j] <- stepsize[j] - c * 0.234 / (28 + i)
      }
    }
    # k_i
    for(j in 1:N){
      p <- numeric(K)
      for(k in 1:K){
        p[k] <-  0.5 * log(det(draws[[1]][[k]]$varInv[,,i-1])) - 0.5 * (draws[[j+1]]$theta[i,] - draws[[1]][[k]]$mean[i-1,]) %*% 
          draws[[1]][[k]]$varInv[,,i-1] %*% (draws[[j+1]]$theta[i,] - draws[[1]][[k]]$mean[i-1,])
      }
      
      if(sum(exp(p)) == 0){
        draws[[j+1]]$k[i] <- which.max(p)
      } else {
        prob <- exp(p) / sum(exp(p))
        draws[[j+1]]$k[i] <- base::sample(1:K, 1, prob=prob)
      }
    }
    
    # thetaHat_k
    for(k in 1:K){
      ki <- 0
      for(j in 1:N){
        if(draws[[j+1]]$k[i] == k){
          ki <- ki + 1
        }
      }
      var <- hyper[[k]]$varInv + ki * draws[[1]][[k]]$varInv[,,i-1]
      var <- solve(var)
      mean <- var %*% hyper[[k]]$varInv %*% hyper[[k]]$mean
      
      for(j in 1:N){
        if(draws[[j+1]]$k[i] == k){
          mean <- mean + var %*% draws[[1]][[k]]$varInv[,,i-1] %*% draws[[j+1]]$theta[i,]
        }
      }
      draws[[1]][[k]]$mean[i,] <- rmvnorm(1, mean, var)
      vardf <- hyper[[k]]$v 
      scaleMat <- hyper[[k]]$scale
      for(j in 1:N){
        if(draws[[j+1]]$k[i] == k){
          vardf <- vardf + 1
          scaleMat <- scaleMat + (draws[[j+1]]$theta[i, ] - draws[[1]][[k]]$mean[i, ]) %*% t((draws[[j+1]]$theta[i, ] - draws[[1]][[k]]$mean[i, ]))
        }
      }
      draws[[1]][[k]]$varInv[,,i] <- rWishart(1, vardf, scaleMat)[,,1]
      
    }
    
    if(i %% 1000 == 0){
      print(paste0('Iteration: ', i, '. Est. Time Remaining: ', round((reps - i) * timePerIter[1] / 60, 2), ' minutes.'))
    }
  }
  print(accept / reps)
  draws
}
