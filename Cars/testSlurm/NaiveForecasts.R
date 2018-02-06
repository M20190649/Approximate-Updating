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
  vConst <- data[sSeq[s] + 1, 3] # Velocity is lagged in the data
  vAvg <- mean(data[-8:1 + sSeq[s], 3 ])
  dAvg <- mean(data[-9:0 + sSeq[s], 2])
  # Model 1: Const V, Const D
  # Model 2: Const V, Zero D
  # Model 3: Const V, Mean D
  # Model 4: Mean V, Const D
  # Model 5: Mean V, Zero D
  # Model 6: Mean V, Mean D
  dxdy <- matrix(0, 2, 6)
  dxdy[,1] <- vConst * c(cos(pi/2 + dConst), sin(pi/2 + dConst))
  dxdy[,2] <- vConst * c(cos(pi/2), sin(pi/2))
  dxdy[,3] <- vConst * c(cos(pi/2 + dAvg), sin(pi/2) + dAvg)
  dxdy[,4] <- vAvg * c(cos(pi/2 + dConst), sin(pi/2 + dConst))
  dxdy[,5] <- vAvg * c(cos(pi/2), sin(pi/2))
  dxdy[,6] <- vAvg * c(cos(pi/2 + dAvg), sin(pi/2 + dAvg))
  
  dist <- matrix(0, H, 6)
  for(h in 1:H){
    for(j in 1:6){
      dist[h, j] <- sqrt((h*dxdy[1, j] - sum(data[sSeq[s] + 1:h, 4]))^2 + (h*dxdy[2, j] - sum(data[sSeq[s] + 1:h, 5]))^2)
    }
  }
  dist <- as.data.frame(dist)
  colnames(dist) <- paste('Naive', 1:6)
  dist$h <- 1:H
  dist$S <- sSeq[s]
  dist$id <- id
  results <- rbind(results, dist)
}

write.csv(results, paste0('naive/car', id, '.csv'), row.names=FALSE)