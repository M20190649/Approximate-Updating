data {
	int<lower=1> N;
	int<lower=1> T;
	real d[T,N];
	real v[T,N];
}

parameters {
	matrix[6, N] theta;
	vector[6] thetaHat1[6];
	vector[6] thetaHat2[6];
	real<lower=0, upper=1> pi[N];
	
	corr_matrix[6] cor1;
	vector<lower=0>[6] sigma1;
	corr_matrix[6] cor2;
	vector<lower=0>[6] sigma2;
}

transformed parameters {
  cov_matrix[6] Sig1;
  cov_matrix[6] Sig2;
  real<lower=0> sigSqD[N];
  real<lower=0> sigSqV[N];
  
  Sig1 = quad_form_diag(cor1, sigma1);
  Sig2 = quad_form_diag(cor2, sigma2);
  
  for(i in 1:N){
    sigSqV[i] = exp(theta[1, i]);
    sigSqD[i] = exp(theta[2, i]);
  }
}
		
model {
  for(i in 1:6){
    sigma1[i] ~ cauchy(0, 5);
    sigma2[i] ~ cauchy(0, 5);
  }
  
  for(i in 1:N){
    pi[i] ~ uniform(0, 1);
    target +=(log_mix(pi[i],
      multi_normal_lpdf(theta[,i] | thetaHat1, Sig1),
      multi_normal_lpdf(theta[,i] | thetaHat2, Sig2)));
    
    for(t in 3:T){
      v[t, i] ~ normal(v[t-1, i] * theta[3, i] + v[t-2, i] * theta[4, i], sqrt(sigSqV[i]));
      d[t, i] ~ normal(d[t-1, i] * theta[5, i] + d[t-2, i] * theta[6, i], sqrt(sigSqD[i]));
    }
  }
}
