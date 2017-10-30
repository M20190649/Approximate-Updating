
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat particleFilter (mat data, rowvec theta, int N, int T, double initV){
  // Particles and indices are arranged into T+1 * N matrix
  // Each column is a particle, and each row is a time
  mat a (T, N), d (T, N), v (T, N), x (T, N), y (T, N);
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
  vec lWeights(N); // auxillary weights and sum.
  vec lSum(N);
  vec kDraws(N); //APF indices
  double sigSqV = theta(0);
  double sigSqD = theta(1);
  double sigSqX = theta(2);
  double sigSqY = theta(3);
  double arV = theta(4);
  double arD = theta(5);
  
  double Dens = 0; // actually the log density
  double Pi = 3.141593;
  
  // x0 from the stationary distribution
  a.row(0) = sqrt(sigSqV / (1 - pow(arV, 2))) * randn<rowvec>(N);
  d.row(0) = sqrt(sigSqD / (1 - pow(arD, 2))) * randn<rowvec>(N);
  v.row(0).fill(initV);
  x.row(0).fill(data(0, 0));
  y.row(0).fill(data(0, 1));
  
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

    
    // Calculate Step Ahead Means and k density
    for(int k = 0; k < N; ++k){
      double aAhead = arV * a(t-1, B(t, k));
      double dAhead = arD * d(t-1, B(t, k));
      double vAhead = v(t-1, B(t, k))  +  aAhead;
      double xAhead = x(t-1, B(t, k))  +  vAhead * cos(dAhead + Pi/2);
      double yAhead = y(t-1, B(t, k))  +  vAhead * sin(dAhead + Pi/2);
      lWeights(k) = -0.5 * log(2 * Pi * sigSqX)  -  pow(data(t, 0) - xAhead, 2) / (2*sigSqX)  -
        0.5 * log(2 * Pi * sigSqY)  -  pow(data(t, 1) - yAhead, 2) / (2*sigSqY);
      if(k == 0){
        lSum(k) = exp(lWeights(k));
      } else {
        lSum(k) = exp(lWeights(k)) + lSum(k-1);
      }
    }
    // Normalise
    for(int k = 0; k < N; ++k){
      lSum(k) = lSum(k) / lSum(N-1);
    }
    // Draw index k using lambda as weights
    for(int k = 0; k < N; ++k){
      double u = randu<vec>(1)[0];
      for(int i = 0; i < N; ++i){
        if(u < lSum(i)){
          kDraws(k) = i;
        }
      }
    }
    // Calculate weights
    double omegaSum = 0;
    for(int k = 0; k < N; ++k){
      // Step Ahead transitions
      a(t, k) = arV * a(t-1, B(t, kDraws(k)))  +  sqrt(sigSqV) * randn<vec>(1)[0];
      d(t, k) = arD * d(t-1, B(t, kDraws(k)))  +  sqrt(sigSqD) * randn<vec>(1)[0];
      v(t, k) = v(t-1, B(t, kDraws(k)))  +  a(t, k);
      x(t, k) = x(t-1, B(t, kDraws(k)))  +  v(t, k) * cos(d(t, k) + Pi/2);
      y(t, k) = y(t-1, B(t, kDraws(k)))  +  v(t, k) * sin(d(t, k) + Pi/2);
      // Measurement density
      omega(k) =  - 0.5 * log(2 * Pi * sigSqX)  -  pow(data(t, 0) - x(t, k), 2) / (2*sigSqX)  -
        0.5 * log(2 * Pi * sigSqY)  -  pow(data(t, 1) - y(t, k), 2) / (2*sigSqY);
      omega(k) -=  lWeights(kDraws(k));
      // sum of weights for normalisation
      omegaSum += exp(omega(k));
    }
    // Normalise weights
    for(int k = 0; k < N; ++k){
      pi(k) = exp(omega(k)) / omegaSum;
    }
    // log p(y_1:T | theta)) = sum_t log p(y_t | theta)
    Dens += log(omegaSum / N);
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
  // Output is { p(y|theta), x_0, x_1, ..., x_T } where {x_t} is a vector giving the path of the particle given by the sampled index.
  mat output(T+1, 2);
  output.row(0).fill(Dens);
  for(int i = 0; i < T; ++i){
    output(i+1, 0) = a(i, B(i, index));
    output(i+1, 1) = d(i, B(i, index));
  }
  return output;
}

double priorDensity (rowvec theta, vec hyperParams){
  // Evaluate log(p(theta))
  double prior = -(hyperParams(0) + 1) * log(theta(0))  -  hyperParams(1) / theta(0)  -     //sigSqA ~ IG 
    (hyperParams(2) + 1) * log(theta(1))  -  hyperParams(3) / theta(1)  +                   //sigSqD ~ IG
    (hyperParams(4) + 1) * log(theta(2))  -  hyperParams(5) / theta(2) -                    //sigSqX ~ IG
    (hyperParams(6) + 1) * log(theta(3))  -  hyperParams(7) / theta(3);                     //sigSqY ~ IG
    (hyperParams(8) - 1) * log(theta(4)+1)  +  (hyperParams(9)-1) * log(1-theta(4))  +      //arV ~ St.B
    (hyperParams(10) - 1) * log(theta(5)+1)  +  (hyperParams(11)-1) * log(1-theta(5));      //arD ~ St.B
  return prior;
}

// [[Rcpp::export]]
Rcpp::List PMMH(mat data, int reps, int N, vec hyperParams, rowvec stepSize, rowvec initial, double initV){
  int T = data.n_rows;
  mat theta(reps, 6);
  mat a(reps, T), d(reps, T);
  theta.row(0) = initial;
  double prevPrior = priorDensity(initial, hyperParams);
  mat prevPF = particleFilter(data, initial, N, T, initV);
  double prevLoglik = prevPF(0, 0);
  a.row(0) = prevPF.col(0).tail(T).t();
  d.row(0) = prevPF.col(1).tail(T).t();
  double accept = 0;
  for(int iter = 1; iter < reps; ++iter){
    rowvec canTheta = theta.row(iter-1) + stepSize % randn<rowvec>(6);
    double canPrior = priorDensity(canTheta, hyperParams);
    mat canPF = particleFilter(data, canTheta, N, T, initV);
    double canLoglik = canPF(0, 0);
    double ratio = exp(canPrior + canLoglik - prevPrior - prevLoglik);
    double u = randu<vec>(1)[0];
    if(u < ratio){
      prevPrior = canPrior;
      prevLoglik = canLoglik;
      theta.row(iter) = canTheta;
      a.row(iter) = canPF.col(0).tail(T).t();
      d.row(iter) = canPF.col(1).tail(T).t();
      accept += 1;
    } else {
      theta.row(iter) = theta.row(iter-1);
      a.row(iter) = a.row(iter-1);
      d.row(iter) = d.row(iter-1);
    }
    if(iter % 250 == 0){
      Rcpp::Rcout << "Iteration " << iter << std::endl;
    }
 
  }
  Rcpp::Rcout << accept/reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta,
                            Rcpp::Named("a") = a,
                            Rcpp::Named("d") = d);
}