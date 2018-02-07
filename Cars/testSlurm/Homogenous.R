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

draws <- readRDS('noHierN2000.RDS')$draws
draws[,1:2] <- exp(draws[,1:2])

H <- 30
S <- 10
maxT <- 300

sSeq <- seq(S, maxT, S)
aLower <- min(data[,1])
if(aLower < 0){
  aLower <- 2 * aLower
} else {
  aLower <- aLower - 0.5
}
dLower <- min(data[,2])
if(dLower < 0){
  dLower <- 2 * dLower
} else {
  dLower <- 0.5 * dLower
}

asup <- seq(aLower, 2*max(data[,1]), length.out = 200)
dsup <- seq(dLower, 2*max(data[,2]), length.out=200)
grid <- cbind(asup, dsup)

results <- data.frame()

# Incrementally add data
for(s in seq_along(sSeq)){
  if((sSeq[s] + H) > nrow(data)){
    break
  }
  MCMCDens(data = data[(sSeq[s]-1):(sSeq[s] + H), ],
           N = 200,
           H = H,
           grid = grid,
           MCMCdraws = draws)
}

write.csv(results, paste0('homog/car', id, '.csv'), row.names=FALSE)