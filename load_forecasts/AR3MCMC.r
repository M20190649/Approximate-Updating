library(readr)
library(tidyr)
library(dplyr)
library(lubridate)
library(Rcpp)
library(microbenchmark)
library(ggplot2)
library(coda)
sourceCpp("AR3MCMC.cpp")
sourceCpp("AR4MCMC.cpp")

#Data

loaddata = select(read_csv("fulldata.csv"), SETTLEMENTDATE, TOTALDEMAND, Date, Period, Day, Temp)
loaddata$SETTLEMENTDATE = ymd_hms(loaddata$SETTLEMENTDATE, tz = "AEST")

GFriday = c("2011-04-22", "2012-04-06", "2013-03-29", "2014-04-18")
EasterS = c("2011-04-24", "2012-04-08", "2013-03-31", "2014-04-20")
EasterM = c("2011-04-26", "2012-04-09", "2013-04-01", "2014-04-21")
Ausday = c("2011-01-26", "2012-01-26", "2013-01-28", "2013-01-27")
Labour = c("2011-03-14", "2012-03-12", "2013-03-11", "2014-03-10")
Anzac = c("2011-04-25", "2012-04-25", "2013-04-25", "2014-04-25")
Queens = c("2011-06-13", "2012-06-11", "2013-06-10", "2014-06-09")
MelCup = c("2011-11-01", "2012-11-06", "2013-11-05", "2014-11-04")
xmasshutdown = c("2011-12-24", "2011-12-27", "2011-12-28", "2011-12-29", 
                 "2011-12-30", "2012-01-02", "2012-01-03", "2012-12-24", 
                 "2012-12-27", "2012-12-28", "2012-12-29", "2012-12-30",
                 "2013-01-02", "2013-01-03", "2013-12-23", "2013-12-24", 
                 "2013-12-27", "2013-12-28", "2013-12-29", "2013-12-30",
                 "2014-01-02", "2014-01-03", "2014-12-22", "2014-12-23", 
                 "2014-12-24", "2014-12-27", "2014-12-28", "2014-12-29", 
                 "2014-12-30")

loaddata = mutate(loaddata, Christmas = ifelse(substr(Date, 6, 10) == "12-25", 1, 0), Boxingday = ifelse(substr(Date, 6, 10) == "12-26", 1, 0),
                  NYE = ifelse(substr(Date, 6, 10) == "12-31", 1, 0), NYD = ifelse(substr(Date, 6, 10) == "01-01", 1, 0),
                  Easter = ifelse(substr(Date, 1, 10) %in% EasterS, 1, 0), GoodFriday = ifelse(substr(Date, 1, 10)%in% GFriday, 1, 0), 
                  OtherHoliday = ifelse(substr(Date, 1, 10) %in% c(EasterM, Ausday, Labour, Anzac, Queens, MelCup, xmasshutdown), 1, 0))

loaddata = mutate(loaddata, AbsTemp = abs(Temp - 20))


X = model.matrix(TOTALDEMAND ~ Day + Christmas + Boxingday + NYE + NYD + Easter + GoodFriday + OtherHoliday + AbsTemp, data = loaddata)
Y = loaddata$TOTALDEMAND

#MCMC Algorithm

MCMC = MCMC_load(Y, X, 25000)
#MCMC = MCMC_load4(Y, X, 5000)

#Diagnostics
posterior.statistics = function(x){
  l95 = quantile(x, probs = 0.025)
  mean = mean(x)
  median = median(x)
  u95 = quantile(x, probs = 0.975)
  return(c(l95, mean, u95))
}

keep = seq(5001, 25000, 40)
MCMCdf = as.data.frame(MCMC[keep,])
colnames(MCMCdf) = c(colnames(X), "Phi1", "Phi48", "Phi336", "Sigmasq")
MCMCdf$iter = keep
MCMC.l = gather(MCMCdf, parameter, draw, -iter)
ggplot(MCMC.l, aes(iter, draw)) + geom_line() + facet_wrap(~parameter, scales = "free") + theme(axis.ticks.x = element_blank(), 
                                                                                                axis.text.x = element_blank(),
                                                                                                axis.ticks.y = element_blank(),
                                                                                                axis.text.y = element_blank())

effectiveSize(MCMCdf)
posterior = apply(MCMCdf[,1:19], 2, posterior.statistics)
row.names(posterior) = c("Lower 95", "Mean", "Median", "Upper 95")
posterior
posmean = posterior[2,]
xbetahat = colSums(t(X) * posmean[1:15])
demeaned = Y - xbetahat
qplot(loaddata$SETTLEMENTDATE, demeaned, geom = "line")

yhat = vector(length = length(Y))
yhat[1:336] = Y[1:336]
for(i in 337:length(Y)){
  yhat[i] = xbetahat[i] + posmean[16]*(yhat[i-1] - xbetahat[i-1]) + posmean[17]*(yhat[i-48] - xbetahat[i-48]) + posmean[18]*(yhat[i-336] - xbetahat[i-336]) 
}
qplot(loaddata$SETTLEMENTDATE, yhat, geom = "line")


ypred = vector(length = length(Y))
ypred[1:336] = Y[1:336]
for(i in 337:length(Y)){
  ypred[i] = xbetahat[i] + posmean[16]*(Y[i-1] - xbetahat[i-1]) + posmean[17]*(Y[i-48] - xbetahat[i-48]) + posmean[18]*(Y[i-336] - xbetahat[i-336]) 
}
qplot(loaddata$SETTLEMENTDATE, Y-ypred, geom = "line")
mean(abs(Y-ypred))

