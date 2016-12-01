library(readr)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)

#Data Manipulation

temps = read_csv("data/temp_86071.csv")
temps$Date = as.Date(temps$Date, origin = "1899-12-30")
tempsubset = filter(temps, Date > "2010-12-31" & Date < "2015-01-01")
ltemp = temps[which(temps$Date > "2010-12-31" & temps$Date < "2015-01-01")-1,]
month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
year = 2011:2014
elecdata = NULL
for(i in 1:4){
  for(j in 1:12){
    elecdata = rbind(elecdata, read_csv(paste0("data/PRICE_AND_DEMAND_", year[i], month[j], "_VIC1.csv")))
  }
}

elecdata$Date = as.Date(elecdata$SETTLEMENTDATE)
period = as.numeric(as.factor(substr(elecdata$SETTLEMENTDATE, 12, 16))) -1
period[period == 0] = 48
elecdata$Period = as.factor(period)
elecdata$Day = wday(elecdata$Date, label = TRUE)
elecdata$Temp = tempsubset$Temp
elecdata$LagTemp = ltemp$Temp
write.csv(elecdata, "fulldata.csv", row.names = FALSE)

#Temperature vs Load
tempeffect = ungroup(mutate(group_by(elecdata, Period, Day), load = TOTALDEMAND - mean(TOTALDEMAND)))
ggplot(tempeffect, aes(Temp, load)) + geom_point()

#Daily Load vs Time
daily = ungroup(summarise(group_by(elecdata, Date), load = mean(TOTALDEMAND), maxtemp = max(Temp)))
ggplot(daily, aes(Date, load)) + geom_line()

#Half-hourly Load vs Time
halfhourly = ungroup(summarise(group_by(elecdata, Period), load = mean(TOTALDEMAND)))
ggplot(halfhourly, aes(as.numeric(Period), load)) + geom_line()

#Load vs Temp
p1 = ggplot(daily, aes(Date, load)) + geom_line()
p2 = ggplot(daily, aes(Date, maxtemp)) + geom_line()
grid.arrange(p1, p2, ncol = 1)



