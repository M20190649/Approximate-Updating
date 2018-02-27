carsAug <- readr::read_csv('carsAug.csv')
library(tidyverse)

carsAug %>%
  group_by(ID) %>%
  mutate(n = seq_along(ID),
         lv = ifelse(n == 1, relv, lag(relv)),
         timeD = time - time[which.max(abs(reldelta - pi/2))],
         timeA = time - time[which.max(abs(relv - lv))]) -> timeToAction

timeToAction %>%
  filter(timeD >= -1000 & timeD <= 500) %>%
  group_by(ID) %>%
  mutate(percentDelta = (reldelta - pi/2) / (max(reldelta) - pi/2)) %>%
  filter(abs(percentDelta) <= 1) %>%
  ungroup() %>%
  group_by(timeD) %>%
  summarise(mean = mean(percentDelta)) %>%
  ggplot() + geom_line(aes(timeD, 100 * mean)) + 
  labs(x = 'Time to Maximum Absolute Angle Deviation (milliseconds)', y = 'Percent of Maximum Angle Deviation') + 
  theme_bw() -> maxAngle
  
timeToAction %>%
  filter(timeA >= -1000 & timeA <= 500) %>%
  group_by(ID) %>%
  mutate(percentA = abs(relv - lv) / max(abs(relv - lv))) %>%
  filter(abs(percentA) <= 1) %>%
  ungroup() %>%
  group_by(timeA) %>%
  summarise(mean = mean(percentA)) %>%
  ggplot() + geom_line(aes(timeA, 100 * mean)) + 
  labs(x = 'Time to Maximum Absoulute Acceleration (milliseconds)', y = 'Percent of Maximum Absolute Acceleration') + 
  theme_bw() -> maxAccel

gridExtra::grid.arrange(maxAccel, maxAngle, ncol = 2)


