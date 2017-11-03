library(tensorflow)
library(greta)
install_tensorflow()
install_tensorflow(version="gpu")

load('standata.rds')

ar1_cov <- function (sigSq, phi, times1, times2 = times1) {
  time_diff <- abs(outer(times1, times2,`-`))
  (sigSq * phi ^ time_diff) / (1 - phi ^ 2)
}

T <- data$T
N <- data$N
x <- as_data(data$x[2:T,])
y <- as_data(data$y[2:T,])

thetaHat <- multivariate_normal(c(-8, -8, 0, 0), diag(1, 4))
thetaVar <- inverse_gamma(1, 0.001, dim = 4)

theta <- zeros(4, N)
for(i in 1:4){
  theta[i,] <- normal(thetaHat[i], sqrt(thetaVar[i]), dim = N)
}
mError <- inverse_gamma(1, 0.01, dim = 2)

a <- zeros(N, T-1)
d <- zeros(N, T-1)
v <- zeros(T-1, N)
xA <- zeros(T-1, N)
yA <- zeros(T-1, N)

SigmaA <- zeros(T-1, (T-1)*N)
SigmaD <- zeros(T-1, (T-1)*N)

for(i in 1:N){
  assign(paste0('SigmaA', i), ar1_cov(exp(theta[i, 1]), theta[i, 3], seq_len(T-1)))
  assign(paste0('SigmaD', i), ar1_cov(exp(theta[i, 2]), theta[i, 4], seq_len(T-1)))
  a[i,] <- multivariate_normal(rep(0, T-1), get(paste0('SigmaA', i)))
  d[i,] <- multivariate_normal(rep(0, T-1), get(paste0('SigmaD', i)))
  v[,i] <- data$initV[o] + cumsum(t(a[i,]))
  xA[,i] <- data$x[1, i] + cumsum(v[,i] * cos(t(d[i,]) + pi/2))
  yA[,i] <- data$x[1, i] + cumsum(v[,i] * cos(t(d[i,]) + pi/2))
}
distribution(x) <- normal(xA, sqrt(mError[1]))
distribution(y) <- normal(yA, sqrt(mError[2]))

gretaMod <- model(mError, thetaHat, thetaVar, theta, a, d)
plot(gretaMod)
draws <- mcmc(gretaMod, n_samples = 10000, warmup = 5000)


## Single

T <- data$T
x1 <- as_data(data$x[2:T,1])
y1 <- as_data(data$y[2:T,1])

sigSq <- inverse_gamma(1, 1, dim = 4)
phiV <- beta(1, 1)
phiD <- beta(1, 1)

SigmaA <- ar1_cov(sigSq[1], phiV, seq_len(T-1))
SigmaD <- ar1_cov(sigSq[2], phiD, seq_len(T-1))
a <- multivariate_normal(rep(0, T-1), SigmaA)
d <- multivariate_normal(rep(0, T-1), SigmaD)
v <- data$initV[1] + cumsum(t(a))
xA <- data$x[1, 1] + cumsum(v * cos(t(d) + pi/2))
yA <- data$y[1, 1] + cumsum(v * sin(t(d) + pi/2))

distribution(x1) <- normal(xA, sqrt(sigSq[3]))
distribution(y1) <- normal(yA, sqrt(sigSq[4]))

gretaMod <- model(sigSq, phiV, phiD, a, d, v, xA, yA)
#plot(gretaMod)
draws <- mcmc(gretaMod, n_samples = 100000, warmup = 25000)

