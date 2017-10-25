data {
	int<lower=0> N;
	int<lower=0> T;
	real d[T,N];
	real v[T,N];
	real hyperMean[8];
	real hyperVar[8];
}

parameters {
	real arD[N];
	real arV[N];
	real Lsigma2D[N];
	real Lsigma2V[N];
	real thetaHat[4];
	real<lower=0> thetaVar[4];
}

transformed parameters {
	real<lower=0> sigma2D[N];
	real<lower=0> sigma2V[N];
	for(n in 1:N){
		sigma2D[n] = exp(Lsigma2D[n]);
		sigma2V[n] = exp(Lsigma2V[n]);
	}
}
		

model {
	for(i in 1:4){
		thetaHat[i] ~ normal(hyperMean[2*i - 1], sqrt(hyperMean[2*i]));
		thetaVar[i] ~ gamma(hyperVar[2*i - 1], hyperVar[2*i]);
	}
		
	
	for(n in 1:N){
		arD[n] ~ normal(thetaHat[1], sqrt(thetaVar[1]));
		arV[n] ~ normal(thetaHat[2], sqrt(thetaVar[2]));
		Lsigma2D[n] ~ normal(thetaHat[3], sqrt(thetaVar[3]));		
		Lsigma2V[n] ~ normal(thetaHat[4], sqrt(thetaVar[4]));
		
		for(t in 2:T){
			d[t, n] ~ normal(arD[n] * d[t-1, n], sqrt(sigma2D));
			v[t, n] ~ normal(arV[n] * v[t-1, n], sqrt(sigma2V));
		}
	}
}



