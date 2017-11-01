library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('mixtureMCMCMetHast.cpp')

mixtureMCMC <- function(data, reps){

  hyper <- list(mean1 = c(-5, -5, rep(0, 4)), var1Inv = solve(diag(c(0.1, 0.1, rep(5, 4)))),  v1 = 10, scale1 = diag(6), 
                mean2 = c(-8, -8, rep(0, 4)), var2Inv = solve(diag(c(0.1, 0.1, rep(5, 4)))),  v2 = 10, scale2 = diag(6))

  draws <- list(list(mean1 = matrix(c(-3, -3, 1, -0.5, 1, -0.5), reps, 6, byrow = TRUE), var1Inv = array(0, dim = c(6, 6, reps)),
                     mean2 = matrix(c(-6, -6, 0.5, 0, 0.5, 0), reps, 6, byrow = TRUE),   var2Inv = array(0, dim = c(6, 6, reps))))
  diag(draws[[1]]$var1[,,1]) = 0.1
  diag(draws[[1]]$var2[,,1]) = 0.1
  for(i in 1:N){
    draws[[i+1]] <- list(theta = matrix(0, reps, 6), k = sample(1:2, reps, replace=TRUE))
  }
  N <- length(data)
  stepsize <- rep(0.01, N)
  accept = rep(0, N)
  
  for(i in 2:reps){
    # thetaHat_1
    k = 0
    for(j in 1:N){
      if(draws[[j+1]]$k[i-1] == 1){
        k = k + 1
      }
    }
    var = hyper$var1Inv + k * draws[[1]]$var1Inv[,,i-1]
    var = solve(var)
    mean = var %*% hyper$var1Inv %*% hyper$mean1
    for(j in 1:N){
      if(draws[[j+1]]$k[i-1] == 1){
        mean = mean + var %*% draws[[1]]$var1Inv[,,i-1]  %*% draws[[j+1]]$theta[i-1,]
      }
    }
    draws[[1]]$mean1[i,] <- rmvnorm(1, mean, var)
    
    vardf <- hyper$v1 
    scaleMat <- hyper$scale1
    for(j in 1:N){
      if(draws[[j+1]]$k[i-1] == 1){
        vardf = vardf + 1
        scaleMat = scaleMat + (draws[[j+1]]$theta[i-1, ] - draws[[1]]$mean1[i, ]) %*% t((draws[[j+1]]$theta[i-1, ] - draws[[1]]$mean1[i, ]))
      }
    }
    draws[[1]]$var1Inv[,,i] <- rWishart(1, vardf, scaleMat)[,,1]
    # thetaHat_2
    k = 0
    for(j in 1:N){
      if(draws[[j+1]]$k[i-1] == 2){
        k = k + 1
      }
    }
    var = hyper$var2Inv + k * draws[[1]]$var2Inv[,,i-1]
    var = solve(var)
    mean = var %*% hyper$var2Inv %*% hyper$mean2
    for(j in 1:N){
      if(draws[[j+1]]$k[i-1] == 2){
        mean = mean + var %*% draws[[1]]$var2Inv[,,i-1]  %*% draws[[j+1]]$theta[i-1,]
      }
    }
    draws[[1]]$mean2[i,] <- rmvnorm(1, mean, var)
    
    vardf <- hyper$v2 
    scaleMat <- hyper$scale2
    for(j in 1:N){
      if(draws[[j+1]]$k[i-1] == 2){
        vardf = vardf + 1
        scaleMat = scaleMat + (draws[[j+1]]$theta[i-1, ] - draws[[1]]$mean1[i, ]) %*% t((draws[[j+1]]$theta[i-1, ] - draws[[1]]$mean1[i, ]))
      }
    }
    draws[[1]]$var2Inv[,,i] <- rWishart(1, vardf, scaleMat)[,,1]
    
    # theta_i
    for(j in 1:N){
      candidate <- draws[[j+1]]$theta[i-1, ] +  stepsize[j] * rnorm(6)
      if(draws[[j+1]]$k[i-1] == 1){
        canDens <- logDensity(data[[j]], candidate, draws[[1]]$mean1[i,], draws[[1]]$var1Inv[,,i])
        oldDens <- logDensity(data[[j]], draws[[j+1]]$theta[i-1,], draws[[1]]$mean1[i,], draws[[1]]$var1Inv[,,i])
      } else {
        canDens <- logDensity(data[[j]], candidate, draws[[1]]$mean2[i,], draws[[1]]$var2Inv[,,i])
        oldDens <- logDensity(data[[j]], draws[[j+1]]$theta[i-1,], draws[[1]]$mean2[i,], draws[[1]]$var2Inv[,,i])
      }
      ratio = exp(canDens - oldDens)
      alpha <- - qnorm(0.234/2)
      c <- stepsize[j] * ((1 - 1/6) * sqrt(2*pi) * exp(alpha^2/2) / (2 * alpha)
                          + 1 / (6 * 0.234 * (1 - 0.234)))
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
      p1 <- 0.5 * log(det(draws[[1]]$var1Inv[,,i])) - 0.5 * (draws[[j+1]]$theta[i,] - draws[[1]]$mean1[i,]) %*% draws[[1]]$var1Inv[,,i] %*% (draws[[j+1]]$theta[i,] - draws[[1]]$mean1[i,])
      p2 <- 0.5 * log(det(draws[[1]]$var2Inv[,,i])) - 0.5 * (draws[[j+1]]$theta[i,] - draws[[1]]$mean2[i,]) %*% draws[[1]]$var2Inv[,,i] %*% (draws[[j+1]]$theta[i,] - draws[[1]]$mean2[i,])
      if(p1 > (p2 + 100)){
        draws[[j+1]]$k[i] <- 1
      } else if(p2 > (p1 + 100)){
        draws[[j+1]]$k[i] <- 2
      } else {
        draws[[j+1]]$k[i] <- rbinom(1, 1, exp(p2) / (exp(p1) + exp(p2))) + 1
      }
    }
    if(i %% 100 == 0){
      print(i)
    }
  }
  print(accept / reps)
  draws
}
