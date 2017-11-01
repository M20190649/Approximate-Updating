library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('mixtureMCMCMetHast.cpp')

mixtureMCMC <- function(data, reps, error = 'gaussian'){
  N <- length(data)
  stepsize <- rep(0.01, N)
  accept <- rep(0, N)

  if(error == 'gaussian'){
    hyper <- list(mean1 = c(-8, -8, rep(0, 4)), var1Inv = solve(diag(5, 6)),  v1 = 6, scale1 = diag(1, 6), 
                  mean2 = c(-8, -8, rep(0, 4)), var2Inv = solve(diag(5, 6)),  v2 = 6, scale2 = diag(1, 6))
    
    draws <- list(list(mean1 = matrix(c(-6.7, -8.5, 1.15, -0.5, 1.05, -0.35), reps, 6, byrow = TRUE), var1Inv = array(0, dim = c(6, 6, reps)),
                       mean2 = matrix(c(-7.3, -10, 1.2, -0.65, 1, -0.3), reps, 6, byrow = TRUE), var2Inv = array(0, dim = c(6, 6, reps))))
    diag(draws[[1]]$var1Inv[,,1]) = 10
    diag(draws[[1]]$var2Inv[,,1]) = 10
    for(i in 1:N){
      draws[[i+1]] <- list(theta = matrix(c(-7, -9.25, 1.175, -0.575, 1.025, -0.325), reps, 6, byrow = TRUE), k = sample(1:2, reps, replace=TRUE))
    }
    likelihood <- nlogDensity
    dim <- 6
  } else if(error == 't') {
    hyper <- list(mean1 = c(-8, -8, rep(0, 6)), var1Inv = solve(diag(5, 8)),  v1 = 8, scale1 = diag(0.5, 8), 
                  mean2 = c(-8, -8, rep(0, 6)), var2Inv = solve(diag(5, 8)),  v2 = 8, scale2 = diag(0.5, 8))
    
    draws <- list(list(mean1 = matrix(c(-6.7, -8.5, 1.15, -0.5, 1.05, -0.35, 1, 1), reps, 8, byrow = TRUE), var1Inv = array(0, dim = c(8, 8, reps)),
                       mean2 = matrix(c(-7.3, -10, 1.2, -0.65, 1, -0.3, 1, 1), reps, 8, byrow = TRUE), var2Inv = array(0, dim = c(8, 8, reps))))
    diag(draws[[1]]$var1Inv[,,1]) = 10
    diag(draws[[1]]$var2Inv[,,1]) = 10
    for(i in 1:N){
      draws[[i+1]] <- list(theta = matrix(c(-7, -9.25, 1.175, -0.575, 1.025, -0.325, 1, 1), reps, 8, byrow = TRUE), k = sample(1:2, reps, replace=TRUE))
    }
    likelihood <- tlogDensity
    dim <- 8
  } else {
    stop('error must be gaussian or t')
  }
  
  for(i in 2:reps){
    # theta_i
    for(j in 1:N){
      candidate <- draws[[j+1]]$theta[i-1, ] +  stepsize[j] * rnorm(dim)
      if(draws[[j+1]]$k[i-1] == 1){
        canDens <- likelihood(data[[j]], candidate, draws[[1]]$mean1[i-1,], draws[[1]]$var1Inv[,,i-1])
        oldDens <- likelihood(data[[j]], draws[[j+1]]$theta[i-1,], draws[[1]]$mean1[i-1,], draws[[1]]$var1Inv[,,i-1])
      } else {
        canDens <- likelihood(data[[j]], candidate, draws[[1]]$mean2[i-1,], draws[[1]]$var2Inv[,,i-1])
        oldDens <- likelihood(data[[j]], draws[[j+1]]$theta[i-1,], draws[[1]]$mean2[i-1,], draws[[1]]$var2Inv[,,i-1])
      }
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
      p1 <- 0.5 * log(det(draws[[1]]$var1Inv[,,i-1])) - 0.5 * (draws[[j+1]]$theta[i,] - draws[[1]]$mean1[i-1,]) %*% 
        draws[[1]]$var1Inv[,,i-1] %*% (draws[[j+1]]$theta[i,] - draws[[1]]$mean1[i-1,])
      p2 <- 0.5 * log(det(draws[[1]]$var2Inv[,,i-1])) - 0.5 * (draws[[j+1]]$theta[i,] - draws[[1]]$mean2[i-1,]) %*%
        draws[[1]]$var2Inv[,,i-1] %*% (draws[[j+1]]$theta[i,] - draws[[1]]$mean2[i-1,])
      if(p1 > (p2 + 10)){
        draws[[j+1]]$k[i] <- 1
      } else if(p2 > (p1 + 10)){
        draws[[j+1]]$k[i] <- 2
      } else {
        draws[[j+1]]$k[i] <- rbinom(1, 1, exp(p2) / (exp(p1) + exp(p2))) + 1
      }
    }
    
    # thetaHat_1
    k <- 0
    for(j in 1:N){
      if(draws[[j+1]]$k[i] == 1){
        k <- k + 1
      }
    }
    var <- hyper$var1Inv + k * draws[[1]]$var1Inv[,,i-1]
    var <- solve(var)
    mean <- var %*% hyper$var1Inv %*% hyper$mean1
    for(j in 1:N){
      if(draws[[j+1]]$k[i] == 1){
        mean <- mean + var %*% draws[[1]]$var1Inv[,,i-1] %*% draws[[j+1]]$theta[i,]
      }
    }
    draws[[1]]$mean1[i,] <- rmvnorm(1, mean, var)
    
    vardf <- hyper$v1 
    scaleMat <- hyper$scale1
    for(j in 1:N){
      if(draws[[j+1]]$k[i] == 1){
        vardf <- vardf + 1
        scaleMat <- scaleMat + (draws[[j+1]]$theta[i, ] - draws[[1]]$mean1[i, ]) %*% t((draws[[j+1]]$theta[i, ] - draws[[1]]$mean1[i, ]))
      }
    }
    draws[[1]]$var1Inv[,,i] <- rWishart(1, vardf, scaleMat)[,,1]
    # thetaHat_2
    k <- 0
    for(j in 1:N){
      if(draws[[j+1]]$k[i] == 2){
        k <- k + 1
      }
    }
    var <- hyper$var2Inv + k * draws[[1]]$var2Inv[,,i-1]
    var <- solve(var)
    mean <- var %*% hyper$var2Inv %*% hyper$mean2
    for(j in 1:N){
      if(draws[[j+1]]$k[i] == 2){
        mean <- mean + var %*% draws[[1]]$var2Inv[,,i-1]  %*% draws[[j+1]]$theta[i,]
      }
    }
    draws[[1]]$mean2[i,] <- rmvnorm(1, mean, var)
    
    vardf <- hyper$v2 
    scaleMat <- hyper$scale2
    for(j in 1:N){
      if(draws[[j+1]]$k[i] == 2){
        vardf <- vardf + 1
        scaleMat <- scaleMat + (draws[[j+1]]$theta[i, ] - draws[[1]]$mean1[i, ]) %*% t((draws[[j+1]]$theta[i, ] - draws[[1]]$mean1[i, ]))
      }
    }
    draws[[1]]$var2Inv[,,i] <- rWishart(1, vardf, scaleMat)[,,1]
    
    if(i %% 500 == 0){
      print(i)
    }
  }
  print(accept / reps)
  draws
}
