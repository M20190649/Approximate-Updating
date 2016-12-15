library(VineCopula)
library(readr)
library(pscl)
source("SGA_R_functions.R")
sourceCPP("emp_sm_SGA.cpp")
##SGA Algorithm

#Setup:
#SGA parameters, reading data, initialisation of lambda, eta
s = 100 #samples per iteration
k = 5 #minimum number of iterations inside loop
K = 2 #minimum number of iterations outside loop
m = 1000 #maximum number of iterations inside loop
M = 50 #maximum number of iterations outside loop
threshold = 0.001 #convergence threshold
steptuning = 0.1 #for adagrad algorithm, can be vectorised

y = read_csv("simy.csv")$x
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
cop_structure1 = NULL
cop_structure2 = NULL
for(i in 1:10){
  for(j in 1:10){
    if(Vine$par[i, j] != 0){
      cop_structure1 = rbind(cop_structure1, c(i, j))
    }
    if(Vine$par2[i, j] != 0){
      cop_structure2 = rbind(cop_structure2, c(i, j))
    }
  }
}
cop_structure = list(eta1 = cop_structure1, eta2 = cop_structure2)
r = 0
RVM = RVineMatrix(Matrix, family, eta[[1]], eta[[2]])
theta = simtheta(10, RVM, lambda)
LB = ELBO(y, lambda, eta[[1]], eta[[2]], theta, x, Trans, dim, T, lags, 10)
diff = 1

#Main Loop
while(diff > threshold | reps < K) {
  if(r> M){
    break
  }
  r = r + 1

  ##Step One - subloop: Maximise over lambda - can probably be evaluated in cpp (unsure on copulas)
  #Simulate Theta
  #Evalute log joint density and log q density from last lambda, eta
  #Evaluate derivative of log q density wrt lambda
  #Calculate change in lambda
  #Check ELBO Convergence
  lambda = SGA_lambda(y, theta, lamdba, eta, s, LB, Vine, k, M)
  LB = lambda[27]
  ##Step Two: Maximise over eta
  #Essentially the same as Step One, need to rebuild RVM object after each change in eta
  eta = SGA_eta(y, theta, lambda, eta, s, LB, Vine, k, M, cop_structure)
  LBnew = eta[[3]]
  ##Step Three: Check ELBO Convergence, decide to repeat loop or not
  diff = abs(LB - LBnew)
  LB = LBnew
}