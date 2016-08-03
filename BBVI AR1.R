##Black Box Stochastic Gradients
N = 100
M = 1000

#Initialisation
psi = vector(length=N)
lambda = vector(length=N)
sh = vector(length=N)
b = vector(length=N)
yhat = vector(length=N)
gamma = vector(length=N)

#Starting values
psi[1] = 0.5
lambda[1] = 0.3
sh[1] = 100
b[1] = 100
yhat[1] = 2
gamma[1] = 0.3

p = 0.5
r = 0.99

tpsi = 0.8
tsig2 = 0.5
T = 10

y = arima.sim(list(ar=tpsi), T, sd=sqrt(tsig2))
typred = tpsi*y[T] + rnorm(1, 0, sqrt(tsig2))

log_joint = function(y, psi, sigma2, ypred){
  out = -(T+2)*log(sqrt(sigma2)) - 1/(2*sigma2)*(sum((y[2:T]-psi*y[1:(T-1)])^2) + (ypred-psi*y[T])^2)
  return(out)
}

deriv_norm_mean = function(x, mean, var){
  a = log(dnorm(x, mean, sqrt(var)))
  h = 0.000001
  b = log(dnorm(x, (mean+h), sqrt(var)))
  return((b-a)/h)
}

deriv_norm_var = function(x, mean, var){
  a = log(dnorm(x, mean, sqrt(var)))
  h = 0.000001
  b = log(dnorm(x, mean, sqrt(var+h)))
  return((b-a)/h)
}

deriv_IG_shape = function(x, shape, scale) {
  a = log(dgamma(1/x, shape, rate=scale))
  h = 0.000001
  b = log(1/dgamma(x, (shape+h), rate=scale))
  return((b-a)/h)
}

deriv_IG_scale = function(x, shape, scale) {
  a = log(dgamma(1/x, shape, rate=scale))
  h = 0.000001
  b = log(dgamma(1/x, shape, rate=(scale+h)))
  return((b-a)/h)
}

f = vector(length=M)
h = vector(length=M)

for(i in 2:N){
  sample = cbind(rnorm(M, psi[(i-1)], sqrt(lambda[(i-1)])),
                 1/rgamma(M, sh[(i-1)], rate=b[(i-1)]),
                 rnorm(M, yhat[(i-1)], sqrt(gamma[(i-1)])))
  
  #Psi mean
  for(j in 1:M){
    h[j] = deriv_norm_mean(sample[j,1], psi[(i-1)], sqrt(lambda[(i-1)]))
    f[j] = h[j]*(log_joint(y, sample[j,1], sample[j,2], sample[j,3])-dnorm(sample[j,1], psi[(i-1)], sqrt(lambda[(i-1)])))
  }
  change = 1/M*sum(f)
  psi[i] = psi[(i-1)]+p*change
  #Psi variance
  for(j in 1:M){
    h[j] = deriv_norm_var(sample[j,1], psi[(i-1)], sqrt(lambda[(i-1)]))
    f[j] = h[j]*(log_joint(y, sample[j,1], sample[j,2], sample[j,3])-dnorm(sample[j,1], psi[(i-1)], sqrt(lambda[(i-1)])))
  }
  change = 1/M*sum(f)
  if(lambda[(i-1)]+p*change < 0){
    lambda[i] = lambda[(i-1)]
  } else {
  lambda[i] = lambda[(i-1)]+p*change
  }
  #Sigma2 Shape
  for(j in 1:M){
    h[j] = deriv_IG_shape(sample[j,2], sh[(i-1)], b[(i-1)])
    f[j] = h[j]*(log_joint(y, sample[j,1], sample[j,2], sample[j,3])-dgamma(1/sample[j,2], shape = sh[(i-1)], rate = b[(i-1)]))
  }
  change = 1/M*sum(f)
  if(-p*change > sh[(i-1)]){
    sh[i] = sh[(i-1)]
  } else {
  sh[i] = sh[(i-1)]+p*change
  }
  #Sigma2 Scale
  for(j in 1:M){
    h[j] = deriv_IG_scale(sample[j,2], sh[(i-1)], b[(i-1)])
    f[j] = h[j]*(log_joint(y, sample[j,1], sample[j,2], sample[j,3])-dgamma(1/sample[j,2], shape = sh[(i-1)], rate = b[(i-1)]))
  }
  change = 1/M*sum(f)
  if(-p*change > b[(i-1)]){
    b[i] = b[(i-1)]
  } else {
    b[i] = b[(i-1)]+p*change
  }
  #Yhat
  for(j in 1:M){
    h[j] = deriv_norm_mean(sample[j,3], yhat[(i-1)], sqrt(gamma[(i-1)]))
    f[j] = h[j]*(log_joint(y, sample[j,1], sample[j,2], sample[j,3])-dnorm(sample[j,3], yhat[(i-1)], sqrt(gamma[(i-1)])))
  }
  change = 1/M*sum(f)
  yhat[i] = yhat[(i-1)]+p*change
  #Gamma
  for(j in 1:M){
    h[j] = deriv_norm_var(sample[j,3], yhat[(i-1)], sqrt(gamma[(i-1)]))
    f[j] = h[j]*(log_joint(y, sample[j,1], sample[j,2], sample[j,3])-dnorm(sample[j,3], yhat[(i-1)], sqrt(gamma[(i-1)])))
  }
  change = 1/M*sum(f)
  if(gamma[(i-1)]+p*change < 0){
    gamma[i] = gamma[(i-1)]
  } else {
    gamma[i] = gamma[(i-1)]+p*change
  }
  p = p*r
}
