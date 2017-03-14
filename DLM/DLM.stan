data {
  int<lower=1> T;
  vector[T] y;
}
parameters {
  vector[T+1] x;
  real mu;
  real<lower=-1, upper=1> phi;
  real<lower=0> sigmaSqX;
  real<lower=0> sigmaSqY;
}
model {
  mu ~ normal(0, 10);
  phi ~ uniform(-1, 1);
  sigmaSqX ~ inv_gamma(1, 1);
  sigmaSqY ~ inv_gamma(1, 1);
  
  x[1] ~ normal(0, sqrt(sigmaSqX));
  for(t in 1:T)
    x[t+1] ~ normal(phi*x[t], sqrt(sigmaSqX));
  for(t in 1:T)
    y[t] ~ normal(mu + x[t+1], sqrt(sigmaSqY));
}
