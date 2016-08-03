log.norm.dens <- function(n, y, mu, sigma2) {
  out = -(n/2)*log(sigma2) - 1/(2*sigma2)*sum((y-mu)^2)
  return(out)
}

n = 10
y = rnorm(n, 3, 2)
muhat = mean(y)
sig2hat = var(y)*(n-1)/n

mu.deriv.dens = function(n, y, mu, sigma2) {
  a = log.norm.dens(n, y, mu, sigma2)
  h = 0.000001
  b = log.norm.dens(n, y, (mu+h), sigma2)
  return((b-a)/h)
}

sig2.deriv.dens = function(n, y, mu, sigma2) {
  a = log.norm.dens(n, y, mu, sigma2)
  h = 0.000001
  b = log.norm.dens(n, y, mu, (sigma2+h))
  return((b-a)/h)
}

muvec = vector(length=100)
muvec[1] = 3
sig2vec = vector(length=100)
sig2vec[1] = 2
p = 0.5

a = Sys.time()
for(i in 2:100){
  muderiv = mu.deriv.dens(n, y, muvec[(i-1)], sig2vec[(i-1)])
  muvec[i] = muvec[(i-1)] + (1/i)^p*muderiv
  sig2deriv = sig2.deriv.dens(n, y, muvec[(i-1)], sig2vec[(i-1)])
  sig2vec[i] = sig2vec[(i-1)] +  (1/i)^p*sig2deriv
}
b = Sys.time()

