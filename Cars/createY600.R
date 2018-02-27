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


#Operations on data
cars %>%
  group_by(ID) %>%
  summarise(medLane = median(lane),
            changed = any(lane != medLane),
            enterExit = any(lane > 5)) %>%
  ungroup() %>%
  right_join(cars, by = "ID") -> cars

cars %>%
  group_by(ID) %>%
  filter(y > 500 & y < 1700 & changed == FALSE & lane < 6) %>%
  mutate(n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         dist = cumsum(v)) -> Y500
degree <- 10

pred <- NULL
for(i in 1:5){
  cars %>%
    group_by(ID) %>%
    filter(lane == i & changed == FALSE & y > 500 & y < 1700) %>%
    filter(ID %in% head(unique(.$ID), 200)) %>%
    mutate(n = seq_along(frame), 
           xlag = ifelse(n == 1, 0, lag(x)), 
           ylag = ifelse(n == 1, 0, lag(y)),
           v = sqrt((x-xlag)^2 + (y-ylag)^2),
           dist = cumsum(v)) %>%
    ungroup() -> carSubset
  modX <- lm(x ~ poly(dist, degree, raw = TRUE), data=carSubset)
  modY <- lm(y ~ poly(dist, degree, raw = TRUE), data=carSubset)
  Y500 %>%
    filter(lane == i) -> ys
  df <- data.frame(xfit = predict(modX, ys),
                   yfit = predict(modY, ys),
                   ID = ys$ID,
                   time = ys$time)
  pred <- rbind(pred, df)
}

Y500 %>% 
  left_join(pred) %>%
  mutate(relX = sign(x - xfit) * sqrt((x-xfit)^2 + (y-yfit)^2),
         lagRelX = ifelse(n == 1, 0, lag(relX)),
         lagDist = ifelse(n == 1, 0, lag(dist)),
         v = sqrt((relX - lagRelX)^2 + (dist - lagDist)^2),
         delta = atan2(dist - lagDist, relX - lagRelX),
         class = factor(class, levels = 1:3, labels = c('bike', 'car', 'truck'))) %>%
  filter(n > 1 & y > 600 & y < 1600) %>%
  select(ID, changed, enterExit, frame, x, y, relX, dist, lane, delta, v, proceeding, time, class) -> Y600

Y600 %>%
  filter(ID != 10874) -> Y600

write.csv(Y600, 'Y600.csv', row.names = FALSE)

Y600 %>%
  filter(ID %in% head(unique(.$ID), 50)) %>%
  ggplot() + geom_path(aes(relX, dist, group = ID, colour = class))

Y600 %>%
  group_by(ID) %>%
  summarise(class = head(class, 1)) %>%
  .$class %>%
  table()


