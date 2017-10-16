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
  mutate(n = seq_along(frame)) %>%
  filter(n == 1) %>%
  summarise(startLane = lane) %>% 
  ungroup() %>%
  right_join(cars, by = "ID") -> cars



degree <- 20
cars %>%
  group_by(ID) %>%
  filter(startLane == 3 & changed == FALSE & y > 600 & y < 1600) %>%
  filter(ID %in% head(unique(.$ID), 150)) %>%
  mutate(n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         dist = cumsum(v)) %>%
  ungroup() %>%
  cbind(fitX = lm(x ~ poly(dist, degree, raw = TRUE), data=.)$fitted) %>%
  cbind(fitY = lm(y ~ poly(dist, degree, raw = TRUE), data=.)$fitted) %>% 
  mutate(relX = sign(x-fitX)*sqrt((x-fitX)^2 + (y-fitY)^2)) -> carSubset

carSubset %>%
  ggplot() + geom_path(aes(relX, dist, group = ID)) -> p1

carSubset %>%
  ggplot() + geom_path(aes(x, y, group=ID)) + geom_path(aes(fitX, fitY, group = ID), colour = 'red') -> p2

gridExtra::grid.arrange(p1, p2, nrow=1)

fitXMod <- lm(x ~ poly(dist, degree, raw = TRUE), data=carSubset)

cars %>%
  group_by(ID) %>%
  filter(y > 600 & y < 1600) %>%
  mutate(n = seq_along(frame), 
         xlag = ifelse(n == 1, 0, lag(x)), 
         ylag = ifelse(n == 1, 0, lag(y)),
         v = sqrt((x-xlag)^2 + (y-ylag)^2),
         dist = cumsum(v)) -> Y600

Y600$fitX <- predict(fitXMod, Y600)

Y600 %>%
  mutate(relX = x - fitX, #sign(x-fitX)*sqrt((x-fitX)^2 + (y-fitY)^2),
         lagRelX = ifelse(n == 1, 0, lag(relX)),
         lagDist = ifelse(n == 1, 0, lag(dist)),
         v = sqrt((relX - lagRelX)^2 + (dist - lagDist)^2),
         delta = atan2(dist - lagDist, relX - lagRelX),
         class = factor(class, levels = 1:3, labels = c('bike', 'car', 'truck'))) %>%
  filter(n > 1) %>%
  select(ID, changed, enterExit, frame, x, y, relX, dist, lane, delta, v, proceeding, time, class) -> carsAugRel

carsAugRel %>%
  group_by(ID) %>%
  mutate(n = seq_along(frame)) %>%
  filter(n == 1) %>%
  select(ID, lane) %>%
  rename(lane600 = lane) %>%
  right_join(carsAugRel) %>%
  filter(lane600 == 3) -> L3Y600

L3Y600 %>%
  ggplot() + geom_path(aes(relX, y, group = ID, colour = class))

augment = aug(as.matrix(select(L3Y600, -class)), as.matrix(select(carsAugRel, -class)))
augment = as.data.frame(augment)
colnames(augment) = c('Vdiff', 'xydiff', 'ttCol')

L3Y600 = data.frame(L3Y600, augment)

ttChange = timeToChange(as.matrix(select(L3Y600, -class)))
tsChange = timeSinceChange(as.matrix(select(L3Y600, -class)))
L3Y600$ttChange = ttChange[,1]
L3Y600$tsChange = tsChange[,1]

L3Y600 %>%
  group_by(ID) %>%
  summarise(changed = mean(lane) != 3) %>%
  right_join(select(L3Y600, -changed), by = 'ID') -> L3Y600

write.csv(L3Y600, 'L3Y600.csv', row.names=FALSE)
