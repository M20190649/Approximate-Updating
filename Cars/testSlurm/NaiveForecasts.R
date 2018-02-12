rm(list=ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(repenv)

ids <- readRDS('carsID.RDS')

if(i <= 873){
  datafc <- readRDS('ForecastData.RDS')
  data <- datafc[[i]]
  id <- ids$idfc[i]
} else {
  datafc <- readRDS('fcDataChanged.RDS')
  data <- datafc[[i-873]]
  id <- datafc[[501]][i-873]
}

H <- 30
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
results <- data.frame()

# Incrementally add data
for(s in seq_along(sSeq)){
  if((sSeq[s] + H) > nrow(data)){
    break
  }
  # Predict naive movements
  dConst <- data[sSeq[s], 2]
  aConst <- data[sSeq[s] +1, 3] - data[sSeq[s] , 3]
  aAvg <- 0.1 * (data[sSeq[s] + 1, 3] - data[sSeq[s] -9, 3])
  dAvg <- mean(data[-9:0 + sSeq[s], 2])
  # Model 1: Past A, Const D
  # Model 2: Past A, Zero D
  # Model 3: Past A, Mean D
  # Model 4: Mean A, Const D
  # Model 5: Mean A, Zero D
  # Model 6: Mean A, Mean D
  # Model 7-9: Const V (Zero A)
  x0 <- rep(data[sSeq[s], 4], 9)
  y0 <- rep(data[sSeq[s], 5], 9)
  v0 <- rep(data[sSeq[s] + 1, 3], 9) 
  dist <- matrix(0, H, 9)
  for(h in 1:H){
    v0 <- v0 + c(rep(aConst, 3), rep(aAvg, 3), rep(0, 3))
    x0 <- x0 + v0 * rep(c(cos(dConst + pi/2), cos(pi/2), cos(pi/2 + dAvg)), 3)
    y0 <- y0 + v0 * rep(c(sin(dConst + pi/2), sin(pi/2), sin(pi/2 + dAvg)), 3)
    dist[h, ] <- sqrt((x0 - sum(data[sSeq[s] + 1:h, 4]))^2 + (y0 - sum(data[sSeq[s] + 1:h, 5]))^2)
  }
  dist <- as.data.frame(dist)
  colnames(dist) <- paste('Naive', 1:9)
  dist$h <- 1:H
  dist$S <- sSeq[s]
  dist$id <- id
  results <- rbind(results, dist)
}

write.csv(results, paste0('naive/car', id, '.csv'), row.names=FALSE)
