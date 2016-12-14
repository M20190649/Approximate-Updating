library(VineCopula)
library(readr)
library(pscl)

##SGA Algorithm

#R only functions
invert_cdf = function(u, lambda, i){
  if(i == 1 | i == 2){
    support = seq(0, 1, 0.001)
    dens = lambda[5*i] * dnorm(support, lambda[5*(i-1) + 1], lambda[5*(i-1) + 2]) + 
      (1-lambda[5*i]) * dnorm(support, lambda[5*(i-1) + 3], lambda[5*(i-1) + 4])
    CDF = cumsum(dens)
    u = u * max(CDF)
    draws = rep(0, length(u))
    for(j in 1:length(u)){
      draws[i] = support[min(which(CDF > u))]
    }
  } else if(i == 3) {
    draws = qigamma(u, lambda[11], lambda[12])
  } else {
    draws = qnorm(u, lambda[5+2*i], lambda[6+2*i])
  }
  return(draws)
}

simtheta= function(s, RVM, lambda){
  #Simulate theta - R only
  #Require RVM Object created from Vine Structure (prespecified) and dependence parameters (eta)
  simvine = RVineSim(s, RVM)
  #Require Inverse CDF transformation for all marginals from mean field parameters (lambda)
  sim = matrix(0, s, 10)
  for(i in 1:10){
    sim[,i] = invert_cdf(simvine[,s], lambda[[s]], i)
  }
  return(sim)
}

#Setup:
#SGA parameters, reading data, initialisation of lambda, eta
s = 100 #samples per iteration
k = 5 #minimum number of iterations inside loop
K = 2 #minimum number of iterations outside loop
m = 1000 #maximum number of iterations inside loop
M = 50 #maximum number of iterations outside loop
threshold = 0.001 #convergence threshold
steptuning = 0.1 #for adagrad algorithm, can be vectorised

y = read_csv("simy.csv")
Vine = readRDS("Vine.rds")

lambdalist = list(alpha = c(),
              delta = c(),
              sigmasq = c(),
              il = c(),
              is1 = c(),
              is2 = c(),
              is3 = c(),
              is4 = c(),
              is5 = c(),
              is5 = c()) #use list for easy initialisation
lambda = unlist(lambda)
eta = list(eta1 = Vine$par,
           eta2 = Vine$par2) #keep in same format as the RVM object
cop_structure = matrix(0, 45, 2)
repso = 0
RVM = RVineMatrix(Matrix, family, eta[[1]], eta[[2]])
theta = simtheta(10, RVM, lambda)
LBo = ELBO(y, lambda, eta[[1]], eta[[2]], theta, x, Trans, dim, T, lags, 10)
diffo = 1
#Main Loop
while(diffo > threshold | repso < K) {
  repso = reps + 1

  
  ##Step One - subloop: Maximise over lambda - can probably be evaluated in cpp (unsure on copulas)
  #Simulate Theta
  #Evalute log joint density and log q density from last lambda, eta
  #Evaluate derivative of log q density wrt lambda
  #Calculate change in lambda
  #Check ELBO Convergence
  diffi = 1
  RVM = RVineMatrix(Matrix, family, eta[[1]], eta[[2]])
  repsi = 0
  LBi = LB
  Gt = rep(0, 26)
  while(diffi > threshold | repsi < k){
    if(repsi > M){
      break
    }
    repsi = repsi + 1
    theta = simtheta(s, RVM, lambda)
    update = SGA_lambda(y, theta, lambda, eta[[1]], eta[[2]], steptuning, Gt, cop_structure, Vine$type, s)
    lambda = update[1,]
    Gt = update[2,]
    LBnew = ELBO(y, lambda, eta[[1]], eta[[2]], theta, x, Trans, dim, T, lags, 10)
    diffi = abs(LBnew - LBi)
    Lbi = LBnew
  }
  
  ##Step Two: Maximise over eta
  #Essentially the same as Step One, need to rebuild RVM object after each change in eta
  diffi = 1
  repsi = 0
  Gt = array(0, dim = c(10, 10, 2))
  while(diffi > threshold | repsi < k) {
    if(repsi > M){
      break
    }
    repsi = repsi + 1
    RVM = RVineMatrix(Matrix, family, eta[[1]], eta[[2]])
    theta = simtheta(s, RVM, lambda)
    update = SGA_eta1(y, theta, lambda, eta[[1]], eta[[2]], steptuning, Gt[,,1], cop_structure, Vine$Type, s)
    eta[[1]] = update[,,1]
    Gt[,,1] = update[,,2]
    update = SGA_eta2(y, theta, lambda, eta[[1]], eta[[2]], steptuning, Gt[,,2], cop_structure, Vine$Type, s)
    eta[[2]] = update[,,1]
    Gt[,,2] = update[,,2]
    LBnew = ELBO(y, lambda, eta[[1]], eta[[2]], theta, x, Trans, dim, T, lags, 10)
    diffi = abs(LBnew - LBi)
    LBi = LBnew
  }

  ##Step Three: Check ELBO Convergence, decide to repeat loop or not
  diff = abs(LBi - LBo)
  LBo = LBi
}