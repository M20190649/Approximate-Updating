library(tidyverse)

# Load data

dataset = '101'

cars1 = read.table(paste0(dataset, '/vehicle-trajectory-data/trajectories1.txt'))

cars2 = read.table(paste0(dataset, '/vehicle-trajectory-data/trajectories2.txt'))
cars2[,1] = cars2[,1] + 10000
cars2[,15] = cars2[,15] + 10000

cars3 = read.table(paste0(dataset, '/vehicle-trajectory-data/trajectories3.txt'))
cars3[,1] = cars3[,1] + 20000
cars3[,15] = cars3[,15] + 20000

cars = rbind(cars1, cars2, cars3)
rm(cars1, cars2, cars3)

colnames(cars) = c('ID', 'frame', 'totalFrames', 'time', 'x', 'y', 
                   'globalX', 'globalY', 'length', 'width', 'class',
                   'veloc', 'accel', 'lane', 'proceeding', 'following', 
                   'spacing', 'headway')

cars %>%
  filter(ID %in% head(unique(.$ID), 10)) %>%
  ggplot() + geom_path(aes(x, y, group = ID)) +
  labs(y = NULL, x = 'Lane') + 
  theme_bw() + 
  scale_x_continuous(labels = c('Lane 1', 'Lane 2', 'Lane 3', 'Lane 4', 'Lane 5', 'Entry Lane'),
                     breaks = c(6, 18, 30, 40, 50, 65))

#Operations on data
cars %>%
  mutate(time = time - min(time)) -> cars
cars %>%
  group_by(ID) %>%
  summarise(medLane = median(lane),
            changed = any(lane != medLane),
            enterExit = any(lane > 5)) %>%
  ungroup() %>%
  right_join(cars, by = "ID") -> cars

cars %>% 
  group_by(ID) %>%
  mutate(n = seq_along(time)) %>%
  filter(n == 1) %>%
  select(ID, lane) %>%
  rename(startLane = lane) %>%
  left_join(cars) -> cars

cars %>%
  group_by(ID) %>%
  filter(changed == FALSE & lane < 6) %>% #& y > 500 & y < 1700) %>%
  mutate(n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         delta = atan2(y - ylag, x - xlag),
         dist = cumsum(v)) -> carsAug

cars %>%
  group_by(ID) %>%
  filter(changed == TRUE & enterExit == FALSE) %>%
  mutate(n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         delta = atan2(y - ylag, x - xlag),
         dist = cumsum(v)) -> carsChanged

degree <- 50

carsAug %>%
  group_by(lane) %>%
  filter(ID %in% head(unique(ID), 100)) %>%
  .$ID %>%
  unique() -> splinesID

pred <- NULL
predChanged <- NULL
for(i in 1:5){
  carsAug %>%
    filter(lane == i  & ID %in% splinesID) -> carSubset
  modX <- smooth.spline(carSubset$dist, carSubset$x, df = degree)
  modY <- smooth.spline(carSubset$dist, carSubset$y, df = degree)
  carsAug %>%
    filter(lane == i) -> ys
  df <- data.frame(xfit = predict(modX, ys$dist)$y,
                   yfit = predict(modY, ys$dist)$y,
                   ID = ys$ID,
                   time = ys$time)
  pred <- rbind(pred, df)
  carsChanged %>% 
    filter(startLane == i) -> ys
  df <- data.frame(xfit = predict(modX, ys$dist)$y,
                   yfit = predict(modY, ys$dist)$y,
                   ID = ys$ID,
                   time = ys$time)
  predChanged <- rbind(predChanged, df)
}

carsAug %>% 
  left_join(pred) %>%
  group_by(ID) %>%
  mutate(relX = x - xfit,#sign(x - xfit) * sqrt((x-xfit)^2 + (y-yfit)^2),
         n = seq_along(time),
         lagRelX = ifelse(n == 1, 0, lag(relX)),
         lagDist = ifelse(n == 1, 0, lag(dist)),
         relv = sqrt((relX - lagRelX)^2 + (dist - lagDist)^2),
         reldelta = atan2(dist - lagDist, relX - lagRelX),
         class = factor(class, levels = 1:3, labels = c('bike', 'car', 'truck'))) %>%
  filter(n > 1) %>%# & y > 600 & y < 1600) %>%
  select(ID, changed, enterExit, frame, x, y, relX, lagRelX, lagDist, dist, lane, delta, v, relv, reldelta, proceeding, time, class) -> carsAug

write.csv(carsAug, 'carsAug.csv', row.names = FALSE)

carsChanged %>% 
  left_join(predChanged) %>%
  group_by(ID) %>%
  mutate(relX = x - xfit,#sign(x - xfit) * sqrt((x-xfit)^2 + (y-yfit)^2),
         n = seq_along(time),
         lagRelX = ifelse(n == 1, 0, lag(relX)),
         lagDist = ifelse(n == 1, 0, lag(dist)),
         relv = sqrt((relX - lagRelX)^2 + (dist - lagDist)^2),
         reldelta = atan2(dist - lagDist, relX - lagRelX),
         class = factor(class, levels = 1:3, labels = c('bike', 'car', 'truck'))) %>%
  filter(n > 1) %>%# & y > 600 & y < 1600) %>%
  select(ID, changed, enterExit, frame, x, y, relX, lagRelX, lagDist, dist, lane, delta, v, relv, reldelta, proceeding, time, class) -> carsChanged

write.csv(carsChanged, 'carsChanged.csv', row.names = FALSE)
saveRDS(splinesID, 'splinesID.RDS')





