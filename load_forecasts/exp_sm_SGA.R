library(VineCopula)
library(readr)
library(pscl)

##SGA Algorithm
#Inverse CDF function is easier in R
invert_cdf = function(u, lambda, i){
  if(i == 1 | i == 2){
    support = seq(0, 1, 0.001)
    dens = lambda[5] * dnorm(support, lambda[1], lambda[2]) + (1-lambda[5]) * dnorm(support, lambda[3], lambda[4])
    CDF = cumsum(dens)
    u = u * max(CDF)
    draws = rep(0, length(u))
    for(j in 1:length(u)){
      draws[i] = support[min(which(CDF > u))]
    }
  } else if(i == 3) {
    draws = qigamma(u, lambda[1], lambda[2])
  } else {
    draws = qnorm(u, lambda[1], lambda[2])
  }
  return(draws)
}

#Setup:
#y, x, Trans, lags, VB Initialisation, Vine object loaded
s = 100 #samples per iteration
k = 5 #minimum number of iterations inside loop
K = 2 #minimum number of iterations outside loop
m = 1000 #maximum number of iterations inside loop
M = 50 #maximum number of iterations outside loop
threshold = 0.001 #convergence threshold
steptuning = 0.1 #for adagrad algorithm, can be vectorised

y = read_csv("simy.csv")
lags = c(1, 7)
m = length(lags)
dim = sum(lags) - m +1
x = c(1, rep(-1, dim-1))
Trans = matrix(0, dim, dim)
Trans[1,1] = 1
newstates = c(1, 2+cumsum(lags-1)[-m])
for(i in 2:m){
  if(i < m){
    Trans[newstates[i], newstates[i]:(newstates[i+1]-1)] = -1
  } else {
    Trans[newstates[m], newstates[m]:dim] = -1
  }
}
for(i in (1:dim)[-newstates]){
  Trans[i, i-1] = 1
}
Vine = readRDS("Vine.rds")

lambda = list()
eta = list()
repso = 0
LBo = ELBO(lambda, eta)
diffo = 1
#Main Loop
while(diffo > threshold | repso < K) {
  repso = reps + 1
  ##Step One: Simulate theta - R only
  #Require RVM Object created from Vine Structure (prespecified) and dependence parameters (eta)
  #Matrix given by Vine object
  #Family given by Vine object
  #params given by eta
  RVM = RVineMatrix(Matrix, family, par, par2)
  simvine = RVineSim(s, RVM)
  #Require Inverse CDF transformation for all marginals from mean field parameters (lambda)
  sim = matrix(0, s, 10)
  for(i in 1:10){
    sim[,i] = invert_cdf(simvine[,s], lambda[[s]], i)
  }
  
##Step Two - subloop: Maximise over lambda - can probably be evaluated in cpp (unsure on copulas)
  diffi = 1
  LBi = LBo
  repsi = 0
  while(diffi > threshold | repsi < k) {
    repsi = repsi + 1
    lambdanew = SGA_lambda()
    LBnew = ELBO(lambdanew, eta)
    diffi = abs(LBnew - LBi)
    LBi = LBnew
    lambda = lambdanew
  }
#Evalute log joint density and log q density from last lambda, eta
#Evaluate derivative of log q density wrt lambda
#Calculate change in lambda
#Check ELBO Convergence

##Step Three: Maximise over eta
  diffi = 1
  repsi = 0
  while(diffi > threshold | repsi < k) {
    repsi = repsi + 1
    etanew = SGA_eta(lambda, etanew)
    LBnew = ELBO(lambdanew)
    diffi = abs(LBnew - LBi)
    LBi = LBnew
    eta = etanew
  }
#Essentially the same as Step Two

#Step Four:
#Check ELBO Convergence, decide to repeat loop or not
  diff = abs(LBi - LBo)
  LBo = LBi
}