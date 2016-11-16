library(readr)
library(tidyr)
library(dplyr)
library(lubridate)
library(Rcpp)
library(microbenchmark)
library(ggplot2)
library(coda)
sourceCpp("AR3MCMC.cpp")

#Data

loadjul2011 = read_csv("Vic_mrim_sum_web.csv")
loadapr2014 = read_csv("VIC_20140401_20140630.csv")
loadjul2014 = read_csv("VIC_20140701_20140930.csv")
loadoct2014 = read_csv("VIC_20141001_20141231.csv")
loadjan2015 = read_csv("VIC_20150101_20150331.csv")
loadapr2015 = read_csv("VIC_20150401_20150630.csv")
loadjul2015 = read.csv("VIC_20150701_20150930.csv")
loadoct2015 = read.csv("VIC_20151001_20151231.csv")
loadjan2016 = read_csv("VIC_20160101_20160331.csv")
temps = read_csv("temp_86071.csv")

temps$Date = as.Date(temps$Date, origin = "1899-12-30")
loadjul2011$DAILYT = apply(loadjul2011[,-c(1, 2, 51)], 1, sum)

allloads = rbind(loadjul2011, loadapr2014, loadjul2014, loadoct2014, loadjan2015, loadapr2015, loadjul2015, loadoct2015, loadjan2016)
allloads$SETTD = as.Date(allloads$SETTD, format = "%d/%m/%Y")
loadlong = select(allloads, -DCTC, - DAILYT) %>% gather(hour, load, -c(PROFILEAREA, SETTD))
loadlong = mutate(loadlong, hour = as.numeric(substr(hour, 4, 5)))
loadlong = summarise(group_by(loadlong, hour, SETTD), load = sum(load))

loadsubset = filter(loadlong, SETTD < "2015-03-01")
loadsubset = plyr::arrange(loadsubset, SETTD, hour)
tempsubset = filter(temps, Date > "2011-06-30")

alldata = data.frame(tempsubset, load = loadsubset$load)
alldata = mutate(alldata, day = factor(wday(Date)), Period = factor(Period), AbsTemp = abs(Temp - 18.3))
alldata$trend = 1:dim(alldata)[1]

X = model.matrix(load ~ Period + day + AbsTemp, data = alldata)
Y = alldata$load

#MCMC Algorithm

rep = 50000
dimx = ncol(X)
T = length(Y)

#Beta - 0 intercept,  1-47 half hour effects, 48-53 day of week, 54 temp, dimx linear time trend, base - sunday 12:00
mu = rep(0, dimx)
tau = rep(10000^2, dimx)
gamma = rep(0, 3)
lambda = rep(10, 3)
alphasig = 1
betasig = 1
alphabar = alphasig + T/2 - 168

#calculations
sumxx = xx(X, T, dimx)
sumyx = yx(Y, X, T, dimx)
sumyy = yy(Y, T)

#initial states
beta = rep(0, dimx)
phi = rep(0, 3)
sigmasq = 1

draws = matrix(0, nrow = rep, ncol = dimx + 4)

for(i in 1:(rep/1000)){
  draws[((i-1)*1000+1):(i*1000) , ] = MCMC_looponly(1000, Y, X, T, dimx, sumyy, sumyx, sumxx, beta, phi, sigmasq, mu, tau, lambda, gamma, alphabar, betasig)
  beta = draws[i*1000, 1:dimx]
  phi = draws[i*1000, (dimx + 1) : (dimx + 3)]
  sigmasq = draws[i * 1000, (dimx + 4)]
}


#MCMC = MCMC_load(Y, X, rep)

#Diagnostics
effectiveSize(MCMC[(rep/5+1):rep, ])
posmean = apply(MCMC[(rep/5+1):rep, ], 2, mean)

xbetahat = colSums(t(X) * posmean[1:55])
yhat = vector(length = length(Y))
yhat[1:336] = Y[1:336]
for(i in 337:length(Y)){
  yhat[i] = xbetahat[i] + posmean[56]*(yhat[i-1] - xbetahat[i-1]) + posmean[57]*(yhat[i-48] - xbetahat[i-48]) +
    posmean[58]*(yhat[i-336] - xbetahat[i-336]) 
}
ypred = vector(length = length(Y))
ypred[1:336] = Y[1:336]
for(i in 337:length(Y)){
  ypred[i] = xbetahat[i] + posmean[57]*(Y[i-1] - xbetahat[i-1]) + posmean[58]*(Y[i-48] - xbetahat[i-48]) +
    posmean[59]*(Y[i-336] - xbetahat[i-336]) 
}

plot(alldata$Date, Y-yhat)


