data {
	int<lower=0> N;
	int<lower=0> T;
	real xM[T,N];
	real yM[T,N];
	real initV[N];
	real length[N];
	real hyperMean[12];
	real hyperVar[12];
}

parameters {
	real ar1D[N];
	real ar1V[N];
	real ar2D[N];
	real ar2V[N];
	real Lsigma2D[N];
	real Lsigma2V[N];
	real thetaHat[6];
	real<lower=0> thetaVar[6];
	real a[T, N];
	real d[T, N];
	real<lower=0> sigma2X;
	real<lower=0> sigma2Y;
}

transformed parameters {
	real<lower=0> sigma2D[N];
	real<lower=0> sigma2V[N];
	real v[T, N];
	real xA[T, N];
	real yA[T, N];
	for(n in 1:N){
		sigma2D[n] = exp(Lsigma2D[n]);
		sigma2V[n] = exp(Lsigma2V[n]);
		v[1, n] = initV[n];
		xA[1, n] = xM[1, n];
		yA[1, n] = yM[1, n];
		for(t in 2:length[n]){
		  v[t, n] = v[t-1, n] + a[t, n];
		  xA[t, n] = xA[t-1, n] + v[t, n] * cos(d[t, n] + 1.5706);
		  yA[t, n] = yA[t-1, n] + v[t, n] * sin(d[t, n] + 1.5706);
		}
	}
	
}
		

model {
  sigma2X ~ gamma(1, 1);
  sigma2Y ~ gamma(1, 1);
  
	for(i in 1:6){
		thetaHat[i] ~ normal(hyperMean[2*i - 1], sqrt(hyperMean[2*i]));
		thetaVar[i] ~ gamma(hyperVar[2*i - 1], hyperVar[2*i]);
	}
		
	
	for(n in 1:N){
		Lsigma2D[n] ~ normal(thetaHat[1], sqrt(thetaVar[1]));		
		Lsigma2V[n] ~ normal(thetaHat[2], sqrt(thetaVar[2]));
		ar1D[n] ~ normal(thetaHat[3], sqrt(thetaVar[3]));
		ar1V[n] ~ normal(thetaHat[4], sqrt(thetaVar[4]));
		ar2D[n] ~ normal(thetaHat[5], sqrt(thetaVar[5]));
		ar2V[n] ~ normal(thetaHat[6], sqrt(thetaVar[6]));
		
		a[1, n] ~ normal(0, sqrt(sigma2V));
		d[1, n] ~ normal(0, sqrt(sigma2D));
		a[2, n] ~ normal(ar1V[n] * a[1, n], sqrt(sigma2V));
		d[2, n] ~ normal(ar1D[n] * d[1, n], sqrt(sigma2D));
		
		for(t in 3:length[n]){
			d[t, n] ~ normal(ar1D[n] * d[t-1, n] + ar2D[n] * d[t-2, n], sqrt(sigma2D));
			a[t, n] ~ normal(ar1V[n] * a[t-1, n] + ar2V[n] * a[t-2, n], sqrt(sigma2V));
			xM[t, n] ~ normal(xA[t, n], sqrt(sigma2X));
			yM[t, n] ~ normal(yA[t, n], sqrt(sigma2Y));
		}
	}
}




