library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)
sourceCpp("mcmc.cpp")

#Generate Data
set.seed(753115)
psi <- 5
tmu <- 5
tsig <- 1
y <- tmu + tsig*sqrt((psi-2)/psi)*rt(50, psi)
y <- as.data.frame(y)
ggplot(data=y, aes(x=y, y=..density..)) + geom_histogram(binwidth=0.5)

logpkernel.f = function(y_in, mu_in, sig_in){
  logkernel = -(50+1)*log(sig_in)
  terms = -0.5*(5+1)*log(1+(((y_in - mu_in)/sig_in)^2)/(5-2))
  out = logkernel+sum(terms)
  return(out)
}

#Alg 1 stuff
IG.f = function(sigin, vin, sigsqhatin){
  out = (2/gamma(vin/2))*((vin*sigsqhatin/2)^(vin/2))*exp(-vin*sigsqhatin/(2*sigin^2))/(sigin^(vin+1))
  return(out)
}

mu_MH.f = function(mu_last, sig_in, y_in){
  ybar = mean(y_in)
  sy = sd(y_in) 
  mu_c = rnorm(1,mean=ybar,sd=(sy/sqrt(50))) 
  logalpha = logpkernel.f(y_in, mu_c, sig_in) 
  logalpha = logalpha - logpkernel.f(y_in, mu_last, sig_in)
  logalpha = logalpha - log(dnorm(mu_c, mean = ybar, sd = (sy/sqrt(50)))) 
  logalpha = logalpha + log(dnorm(mu_last, mean = ybar, sd = (sy/sqrt(50))))
  alpha = min(1, exp(logalpha))
  if(alpha==1){
    out = mu_c
  } else if (runif(1,0,1)<=alpha){
    out = mu_c
  } else {
    out = mu_last
  }
  return(out)
}

sig_MH.f = function(mu_in, sig_last, y_in){
  sighatsq = sum((y_in - mu_in)^2)/50
  sig_shape = 25
  sig_scale = 50*sighatsq/2
  sig_c = 1/sqrt(rgamma(1,shape=sig_shape,rate = sig_scale)) 
  logalpha = logpkernel.f(y_in, mu_in, sig_c) 
  logalpha = logalpha - logpkernel.f(y_in, mu_in, sig_last)
  logalpha = logalpha - log(IG.f(sig_c, 50, sighatsq))
  logalpha = logalpha + log(IG.f(sig_last, 50, sighatsq))
  alpha = min(1, exp(logalpha))
  if(alpha==1){
    out = sig_c
  } else if (runif(1,0,1)<=alpha){
    out = sig_c
  } else {
    out = sig_last
  }
  return(out)
}

Alg1 <- function(B, NG, y) {
  MHGibbs = matrix(c(0),nrow=(B+NG),ncol=2) 
  MHGibbs[1,] = c(mean(y),sd(y))
  for(iter in 2:(B+NG)){
    MHGibbs[iter,1] = mu_MH.f(MHGibbs[(iter-1),1], MHGibbs[(iter-1),2], y)
    MHGibbs[iter,2] = sig_MH.f(MHGibbs[iter,1], MHGibbs[(iter-1),2], y)
  }
  return(MHGibbs)
}

Alg1C <- function(B, NG, y) {
  MHGibbs <- MHchain(B, NG, y)
  return(MHGibbs)
}

MHGibbs <- data.frame(Mu = MHGibbs[,1], Sigma = MHGibbs[,2], rep = 1:(B+NG))
p1 <- ggplot(data=MHGibbs, aes(y=Mu, x=rep)) + geom_line(colour="red") +
  geom_vline(xintercept=B, linetype="dotted") + 
  labs(title = "Algorithm 1 Mu estimates", x = "repetition") + ylim(4.4, 5.6)
p2 <- ggplot(data=MHGibbs, aes(y=Sigma, x=rep)) + geom_line(colour="blue") +
  geom_vline(xintercept=B, linetype="dotted") + 
  labs(title = "Algorithm 1 Sigma estimates", x = "repetition") + ylim(0.5, 1.5)
grid.arrange(p1, p2, ncol=1)

