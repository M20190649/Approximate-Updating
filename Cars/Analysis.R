library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(GGally)

# Load data
cars1 = read.table('vehicle-trajectory-data/trajectories1.txt')

cars2 = read.table('vehicle-trajectory-data/trajectories2.txt')
cars2[,1] = cars2[,1] + max(cars1[,1])

cars3 = read.table('vehicle-trajectory-data/trajectories3.txt')
cars3[,1] = cars3[,1] + max(cars2[,1])

cars = rbind(cars1, cars2, cars3)
rm(cars1, cars2, cars3)

colnames(cars) = c('ID', 'frame', 'totalFrames', 'time', 'x', 'y', 
                   'globalX', 'globalY', 'length', 'width', 'class',
                   'veloc', 'accel', 'lane', 'proceeding', 'following', 
                   'spacing', 'headway')

#Operations on data

cars %>%
  mutate(time = time - min(time), 
         xmin = x - 0.5 * width,
         xmax = x + 0.5 * width) -> cars

cars %>%
  group_by(ID) %>%
  summarise(medLane = median(lane),
            changed = any(lane != medLane),
            enterExit = any(lane > 5)) %>%
  ungroup() %>%
  right_join(cars, by = "ID") -> cars

cars %>%
  group_by(ID) %>%
  summarise(endTime = max(time),
            startTime = min(time)) %>%
  left_join(cars, by = 'ID') %>%
  filter(time == startTime) %>%
  select(ID, lane) -> startLanes

startLanes %>%
  filter(lane == 1) -> startIn1

yQuantiles = quantile(cars$y, seq(0, 1, 0.0001))
cars$yQ = sapply(cars$y, function(x) max(which(yQuantiles <= x)))

cars %>%
  filter(lane < 6 & yQ < max(yQ)) %>%
  group_by(lane, yQ) %>%
  summarise(lowerMid = quantile(x, 0.025),
            upperMid = quantile(x, 0.975),
            medianMid = quantile(x, 0.5),
            lowerLeft = quantile(xmin, 0.025),
            medianLeft = quantile(xmin, 0.5),
            upperRight = quantile(xmax, 0.975),
            medianRight = quantile(xmax, 0.5)) %>%
  mutate(yQ = yQuantiles[yQ]) -> lanePath


length(unique(cars$ID))

# Various plots
cars %>% 
  group_by(ID) %>%
  tally() %>%
  ggplot() + geom_histogram(aes(n))

cars %>% 
  filter(ID %in% head(sort(unique(cars$ID)), 10)) %>%
  ggplot() + geom_path(aes(x, y, group = ID, colour = factor(ID))) + theme(legend.position = 'none')

cars %>%
  filter(ID %in% laneChange$ID) %>%
  filter(ID %in% head(sort(unique(.$ID)), 25)) %>%
  ggplot() + geom_path(aes(x, y, group = ID, colour = factor(lane))) + theme(legend.position = 'none')

cars %>%
  filter(ID %in% stayed$ID) %>%
  filter(ID %in% head(sort(unique(.$ID)), 500)) %>%
  ggplot() + geom_path(aes(x, y, group = ID, colour = factor(ID), alpha=0.1)) +
  geom_path(data = lanePath, aes(lowerMid, yQ, group=factor(lane))) +
  geom_path(data = lanePath, aes(upperMid, yQ, group=factor(lane))) +
  theme(legend.position = 'none')

cars %>%
  filter(ID %in% stayed$ID) %>%
  filter(ID %in% head(sort(unique(.$ID)), 100)) %>%
  ggplot() + geom_path(aes(xmin, y, group = ID, colour = factor(ID), alpha=0.1)) +
  geom_path(aes(xmax, y, group = ID, colour = factor(ID), alpha=0.1)) +
  geom_path(data = lanePath, aes(lowerOut, yQ, group=factor(lane))) +
  geom_path(data = lanePath, aes(upperOut, yQ, group=factor(lane))) +
  theme(legend.position = 'none')

cars %>% group_by(ID) %>%
  summarise(startTime = min(time)) %>%
  right_join(cars, by = 'ID') %>%
  filter(ID %in% startIn1$ID &
           enterExit == FALSE) %>%
  filter(ID %in% head(sort(unique(.$ID)), 15)) %>%
  mutate(xL = lanePath$medianUpper[yQ],
         xU = lanePath$medianUpper[yQ+1],
         yL = yQuantiles[yQ],
         yU = yQuantiles[yQ+1],
         xAdj = xL + (y - yL) * (xU - xL) / (yU - yL))  %>%
  ggplot() + geom_path(aes(xmax, time - startTime, group = ID, colour = changed, alpha = 0.01)) #+ 
geom_path(data = filter(lanePath, lane == 1), aes(medianUpper, yQ)) + 
  theme(legend.position = 'none')

grid = seq(-5, 12, length.out = 1000)
dens = dnorm(grid, 4.77086, 1.26847)
df = data.frame(grid, dens) 

cars %>%
  group_by(ID) %>%
  summarise(startTime = min(time)) %>%
  right_join(cars, by = 'ID') %>%
  filter(time == startTime &
           lane == 1) %>%
  ggplot() + geom_histogram(aes(x, y = ..density..)) + geom_line(data=df, aes(grid, dens))



# Generate some fake data
rskew = function(n=1, mu, sigmaSq, nu, delta){
  out = rep(0, n)
  for(i in 1:n){
    w = rgamma(1, nu/2, nu/2)
    z = abs(rnorm(1, 0, sqrt(1/w)))
    out[i] = rnorm(1, mu + delta*z, sqrt(sigmaSq / w))
  }
  out
}
#sourceCpp('carsPMCMC.cpp')
sourceCpp('cars.cpp')
set.seed(34)
T = 500
sigSqX = 0.01
sigSqE = 0.01
nu = 15
phi = 0.9

x = rep(0, T)
vx = rep(0, T)
x[1] = rnorm(1, 4.77, 1.26)
vx[1] = rnorm(1, 0, sqrt(sigSqE/(1-phi^2)))

for(t in 2:T){
  vx[t] = phi*vx[t-1] + sqrt(sigSqE) * rt(1, nu)
  x[t] = x[t-1] + vx[t] + sqrt(sigSqX) * rnorm(1)
}

dx = x[2:T] - x[1:(T-1)]
ggplot() + geom_path(aes(x, 1:T))
#ggplot() + geom_point(aes(vx[2:T], dx))

mean = c(-3, -3, 2.5, 0.5, 0)
var = rep(0.1, 5)
lambda = cbind(mean, diag(var))

VB = VB_Cars(x = x,
             lambdaIn = lambda, 
             hyperParams = c(15, 0.1, 15, 0.02, 2, 0.05, 0.5, 0.1),
             alpha = 0.05,
             S = 1,
             N = 10,
             maxIter = 5000,
             threshold = 0.01,
             thresholdIS = 0.8)
lambda

Sigma = t(VB$U) %*% VB$U





MCMC = PMCMC(x, 10000, 50, c(5, 0.5, 5, 0.5, 2, 0.05, 0, 5), c(0.05, 0.1, 0.1), 0, 0.5)

theta = as.data.frame(MCMC$theta[5001:10000,])
colnames(theta) = c('SigSqX', 'SigSqE', 'Nu', 'Delta')
ggpairs(theta)

# MCMC for realsies
cars %>%
  group_by(ID) %>%
  summarise(endTime = max(time),
            startTime = min(time)) %>%
  ungroup() %>%
  left_join(cars, by = 'ID') %>%
  filter(time <= startTime + 100) %>%
  group_by(ID) %>%
  summarise(veloc = x[2] - x[1]) %>%
  .$veloc -> initialVelocities

t = MASS::fitdistr(initialVelocities, 't')
grid = seq(-0.6, 0.6, length.out=1000)
density = dt((grid - t$estimate[1])/t$estimate[2], t$estimate[3])/(t$estimate[2])
ggplot() + geom_histogram(aes(initialVelocities, y=..density..), binwidth=0.005) + geom_line(aes(grid, density))

cars %>%
  filter(ID == 10) -> car10Path

ggplot(car10Path) + geom_path(aes(x, y))

MCMCcar10 = list(theta = matrix(0, nrow=0, ncol=4), vx = matrix(0, nrow=0, ncol=400))
for(i in 1:25){
  if(i == 1){
    MCMCtemp = PMCMC(x = car10Path$x[1:400], 
                     initialTheta = c(0.1, 0.1, 10, 0),
                     reps = 1000,
                     N = 200,
                     hyperParams = c(5, 0.5, 5, 0.5, 2, 0.05, 0, 5),
                     stepSize = c(0.001, 0.001, 2.5, 0.01),
                     velocMean = t$estimate[1], 
                     velocSd = t$estimate[2], 
                     df = t$estimate[3], 
                     Normal = FALSE)
  } else {
    MCMCtemp = PMCMC(x = car10Path$x[1:400], 
                     initialTheta = MCMCcar10$theta[nrow(MCMCcar10$theta), ],
                     reps = 1000,
                     N = 200,
                     hyperParams = c(5, 0.5, 5, 0.5, 2, 0.05, 0, 5),
                     stepSize = c(0.001, 0.001, 2.5, 0.01),
                     velocMean = t$estimate[1], 
                     velocSd = t$estimate[2], 
                     df = t$estimate[3], 
                     Normal = FALSE)
  }
  MCMCcar10$theta = rbind(MCMCcar10$theta, MCMCtemp$theta)
  MCMCcar10$vx = rbind(MCMCcar10$vx, MCMCtemp$vx)
  save(MCMCcar10, file = 'car10.RDA')
  print(i)
}

MCMCcar10$theta %>%
  as.data.frame() %>%
  cbind(t = 1:25000) %>%
  gather(param, value, -t) %>%
  ggplot() + geom_line(aes(t, value)) + facet_wrap(~param, ncol = 1, scales = 'free')


theta10 = as.data.frame(MCMCcar10$theta[1001:25000, ])
colnames(theta10) = c('SigSqX', 'SigSqE', 'Nu', 'Delta')
ggpairs(theta10)

car10Path = filter(cars, ID==10)
predictedX = matrix(0, 436, 1000)
predictedVX = matrix(0, 436, 1000)

for(t in 2:436){
  for(j in 1:1000){
    draw = sample(1001:25000, 1)
    if(t <= 336){
      predictedVX[t, j] = MCMCcar10$vx[draw, t-1] + rskew(1, 0, MCMCcar10$theta[draw, 2], MCMCcar10$theta[draw, 3], MCMCcar10$theta[draw, 4])
      predictedX[t, j] = predictedVX[t, j] + car10Path$x[t-1] + rnorm(1, 0, sqrt(MCMCcar10$theta[draw, 1]))
    } else {
      predictedVX[t, j] = predictedVX[t-1, j] + rskew(1, 0, MCMCcar10$theta[draw, 2], MCMCcar10$theta[draw, 3], MCMCcar10$theta[draw, 4]) 
      predictedX[t, j] = predictedVX[t, j] + predictedX[t-1, j]  + rnorm(1, 0, sqrt(MCMCcar10$theta[draw, 1]))
    }
  }
}

predictedX %>%
  apply(1, function(x) c(quantile(x, 0.025), quantile(x, 0.25), quantile(x, 0.5), quantile(x, 0.75), quantile(x, 0.975))) %>%
  t() %>%
  as.data.frame() -> predictions
colnames(predictions) = c('l95', 'l50', 'median', 'u50', 'u95') 
predictions$y = car10Path$y
predictions[4:436,] %>% gather(stat, value, -y) -> predictions
ggplot() + geom_path(data=car2Path, aes(x, y)) +
  geom_path(data=predictions, aes(value, y, group=stat), colour = 'red')


