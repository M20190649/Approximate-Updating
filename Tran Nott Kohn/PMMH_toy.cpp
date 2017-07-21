// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
vec particleFilterToy (vec y, rowvec theta, int N, int T){
  // Particles and indices are arranged into T * N matrix
  // Each column is a particle, and each row is a time
  mat X (T, N);
  // B maps the indices used in the resampling
  // If X_t^2 = X_t-1^1, then B(t, 2) = 1.
  // To get a path of a particle, sample column i of B, and then X(t, B(t, i)) tracks particle i at time t. 
  mat B (T, N);
  // At time 0, each particle index points to the original particle
  for(int k = 0; k < N; ++k){
    B(0, k) = k;
  }
  vec pi (N); pi.fill(1.0/N); //normalised weights
  vec piSum(N); // cumulative sum of pi
  vec omega (N); // unnormalised weights
  double sigSqV = theta(0);
  double sigSqW = theta(1);
  double yDens = 0; // actually the log density
  
  // x0 from the stationary distribution
  X.row(0) = sqrt(5) * randn<rowvec>(N);
  
  // Particle Filter loop
  for(int t = 1; t < T; ++t){
    // Create CDF to use for resampling
    piSum(0) = pi(0);
    for(int k = 1; k < N; ++k){
      piSum(k) = pi(k) + piSum(k-1);
    }
    // Resample xt via inverse CDF
    for(int k = 0; k < N; ++k){
      double u = randu<vec>(1)[0];
      for(int i = 0; i < N; ++i){
        if(u < piSum(i)){
          B(t, k) = i;
          break;
        }
      }
    }
    // Calculate weights
    double omegaSum = 0;
    for(int k = 0; k < N; ++k){
      double Xold = X(t-1, B(t, k));
      X(t, k) = Xold/2  +  25 * Xold / (1 + pow(Xold, 2))  +  8 * cos(1.2*(t+1))  +  sqrt(sigSqV) * randn<vec>(1)[0]; // step ahead transition
      omega(k) = 1.0 / sqrt(2*3.14159*sigSqW) * exp(-pow(y(t-1) - pow(X(t, k), 2)/20, 2) / (2*sigSqW)); // y ~ N(0, exp(Xt))
      omegaSum += omega(k); // sum of weights for normalisation
    }
    // Normalise weights
    for(int k = 0; k < N; ++k){
      pi(k) = omega(k) / omegaSum;
    }
    // log p(y_1:T | theta)) = sum_t log p(y_t | theta)
    yDens += log(omegaSum / N);
  } // end main loop
  
  // Randomly sample a column of B
  piSum(0) = pi(0);
  for(int k = 1; k < N; ++k){
    piSum(k) = pi(k) + piSum(k-1);
  }
  double u = randu<vec>(1)[0];
  double index;
  for(int k = 0; k < N; ++k){
    if(u < piSum(k)){
      index = k;
    }
  }
  // Output is { p(y|theta), x_0, x_1, ..., x_T } where {x_t} is the path of the particle given by the sampled index.
  vec output(T+1);
  output(0) = yDens;
  for(int i = 0; i < T; ++i){
    output(i+1) = X(i, B(i, index));
  }
  return output;
}

double priorDensity (rowvec theta, vec hyperParams){
  //theta[0] = sigSqV ~ IG(hp[0], hp[1])
  //theta[1] = sigSqW ~ IG(hp[2], hp[3])
  double density = -(hyperParams(0) + 1) * log(theta(0))- hyperParams(1) / theta(0) -
    (hyperParams(2) + 1) * log(theta(1))- hyperParams(3) / theta(1);
  return density;
}

// [[Rcpp::export]]
Rcpp::List PMMHToy(vec y, int reps, int burn, int N, vec hyperParams, rowvec stepSize){
  int T = y.n_elem;
  mat theta(reps, 2);
  mat x(reps, T);
  rowvec initial = {10, 10};
  theta.row(0) = initial;
  double prevPrior = priorDensity(initial, hyperParams);
  vec prevPF = particleFilterToy(y, initial, N, T);
  x.row(0) = prevPF.tail(T).t();
  
  double accept = 0;
  for(int iter = 1; iter < reps + burn; ++iter){
    
    rowvec canTheta = theta.row(iter-1) + stepSize % randn<rowvec>(2);
    double canPrior = priorDensity(canTheta, hyperParams);
    vec canPF = particleFilterToy(y, canTheta, N, T);
    double ratio = exp(canPrior + canPF[0] - prevPrior - prevPF[0]);
    double u = randu<vec>(1)[0];
    if(u < ratio){
      prevPrior = canPrior;
      prevPF = canPF;
      theta.row(iter) = canTheta;
      x.row(iter) = canPF.tail(T).t();
      if(iter > burn){
        accept += 1;  
      }
    } else {
      theta.row(iter) = theta.row(iter-1);
      x.row(iter) = x.row(iter-1);
    }
  }
  Rcpp::Rcout << accept/reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta.tail_rows(reps), Rcpp::Named("x") = x.tail_rows(reps));
}