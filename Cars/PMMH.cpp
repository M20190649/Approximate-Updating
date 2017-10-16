// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat particleFilter (mat pos, rowvec theta, int N, int T){
  // Particles and indices are arranged into T+1 * N matrix
  // Each column is a particle, and each row is a time
  mat a (T+1, N), d (T+1, N), s (T+1, N), v (T+1, N), x (T+1, N), y (T+1, N);
  // B maps the indices used in the resampling
  // If X_t^2 = X_t-1^1, then B(t, 2) = 1.
  // To get a path of a particle, sample column i of B, and then X(t, B(t, i)) tracks particle i at time t. 
  mat B (T+1, N);
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
  double sigSqA = theta(0);
  double sigSqD = theta(1);
  double alpha = theta(2);
  double beta = theta(3);
  double sigSqX = theta(4);
  double sigSqY = theta(5);
  double lambda = theta(6);
  double phi = theta(7);
  double gamma = theta(8);
  double Dens = 0; // actually the log density
  double Pi = 3.14159;
  
  
  // x0 from the stationary distribution
  a.row(0) =  sqrt(sigSqA / (1 - pow(phi, 2))) * randn<rowvec>(N);
  s.row(0) = 0;
  d.row(0) = Pi/2 +  sqrt(sigSqD) * randn<rowvec>(N);
  v.row(0) = 5.0 + randn<rowvec>(N);
  x.row(0) = 5.0 + randn<rowvec>(N);
  y.row(0) = randn<rowvec>(N);
  
  // Particle Filter loop
  for(int t = 1; t <= T; ++t){
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
    double aAhead = phi * a(t-1, B(t, k));;
    double sAhead = s(t-1, B(t, k));
    double dAhead = Pi/2  +  (1 - sAhead) * lambda * (d(t-1, B(t, k)) - Pi/2)  +
      sAhead * gamma * (x(t-1, B(t, k)) - 5);
    double vAhead = v(t-1, B(t, k))  +  aAhead;
    double xAhead = x(t-1, B(t, k))  +  vAhead * cos(dAhead);
    double yAhead = y(t-1, B(t, k))  +  yAhead * sin(dAhead);
    lWeights(k) = 1.0 / sqrt(2 * Pi * sigSqX) * exp(-pow(pos(t, 0) - xAhead, 2) / (2*sigSqX)) *
      1.0 / sqrt(2 * Pi * sigSqY) * exp(-pow(pos(t, 1) - yAhead, 2) / (2*sigSqY));
    if(k == 0){
      lSum(k) = lWeights(k);
    } else {
      lSum(k) = lWeights(k) + lSum(k-1);
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
      a(t, k) = phi * a(t-1, B(t, kDraws(k)))  +  sqrt(sigSqA) * randn<vec>(1)[0];
      double u = randu<vec>(1)[0];
      if((s(t-1, B(t, kDraws(k)))  == 0 & u < alpha * abs(x(t-1, B(t, kDraws(k)))  - 5)) |
         (s(t-1, B(t, kDraws(k)))  == 1 & u < beta * abs(x(t-1, B(t, kDraws(k)))  - 5))){
        s(t, k) = 1;
      } else {
        s(t, k) = 0;
      }
      d(t, k) =  Pi/2  +  (1 - s(t, k)) * lambda * (d(t-1, B(t, kDraws(k)))  - Pi/2)  +
        s(t, k) * gamma * (x(t-1, B(t, kDraws(k)))- 5)  +  sqrt(sigSqD) * randn<vec>(1)[0];
      v(t, k) = v(t-1, B(t, kDraws(k)))  +  a(t, k);
      x(t, k) = x(t-1, B(t, kDraws(k)))  +  v(t, k) * cos(d(t, k));
      y(t, k) = y(t-1, B(t, kDraws(k)))  +  v(t, k) * sin(d(t, k));
      // Measurement density
      omega(k) =  1.0 / sqrt(2 * Pi * sigSqX) * exp(-pow(pos(t, 0) - x(t, k), 2) / (2*sigSqX)) *
        1.0 / sqrt(2 * Pi * sigSqY) * exp(-pow(pos(t, 1) - y(t, k), 2) / (2*sigSqY));
      // sum of weights for normalisation
      omegaSum += omega(k);
    }
    // Normalise weights
    for(int k = 0; k < N; ++k){
      pi(k) = omega(k) / omegaSum;
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
  mat output(T+2, 6);
  output.row(0) = Dens;
  for(int i = 0; i <=T; ++i){
    output(i+1, 0) = a(i, B(i, index));
    output(i+1, 1) = s(i, B(i, index));
    output(i+1, 2) = d(i, B(i, index));
    output(i+1, 3) = v(i, B(i, index));
    output(i+1, 4) = x(i, B(i, index));
    output(i+1, 5) = y(i, B(i, index));
  }
  return output;
}

double priorDensity (rowvec theta, vec hyperParams){
  // Evaluate log(p(theta))
  double prior = -(hyperParams(0) + 1) * log(theta(0))  -  hyperParams(1) / theta(0)  -     //sigSqA ~ IG 
    (hyperParams(2) + 1) * log(theta(1))  -  hyperParams(3) / theta(1)  +                   //sigSqD ~ IG
    (hyperParams(4) - 1) * log(theta(2))  +  (hyperParams(5)-1) * log(1-theta(2))  +        //lambda ~ B
    (hyperParams(6) - 1) * log(theta(3))  +  (hyperParams(7)-1) * log(1-theta(3))  +        //phi ~ B
    (hyperParams(8) - 1) * log(theta(4))  -  hyperParams(9) * theta(4)  +                   //alpha ~ G
    (hyperParams(10) - 1) * log(theta(5))  -  hyperParams(11) * theta(5)  -                 //beta ~ G
    pow(theta(6) - hyperParams(12), 2) / (2*hyperParams(13))  -                             //gamma ~ N
    (hyperParams(14) + 1) * log(theta(7))  -  hyperParams(15) / theta(7) -                  //sigSqX ~ IG
    (hyperParams(16) + 1) * log(theta(8))  -  hyperParams(17) / theta(8);                   //sigSqY ~ IG
  return prior;
}

// [[Rcpp::export]]
Rcpp::List PMMH(mat pos, int reps, int N, vec hyperParams, rowvec stepSize){
  int T = pos.n_rows;
  mat theta(reps, 9);
  mat a(reps, T+1), s(reps, T+1), d(reps, T+1), v(reps, T+1), x(reps, T+1), y(reps, T+1);
  rowvec initial = {0.2, 0.2, 0.9, 0.9, 0.5, 0.5, 0.001, 0.1, 0.1};
  theta.row(0) = initial;
  double prevPrior = priorDensity(initial, hyperParams);
  vec prevPF = particleFilter(pos, initial, N, T);
  a.row(0) = prevPF.col(0).tail(T+1).t();
  s.row(0) = prevPF.col(1).tail(T+1).t();
  d.row(0) = prevPF.col(2).tail(T+1).t();
  v.row(0) = prevPF.col(3).tail(T+1).t();
  x.row(0) = prevPF.col(4).tail(T+1).t();
  y.row(0) = prevPF.col(5).tail(T+1).t();
  
  double accept = 0;
  for(int iter = 1; iter < reps; ++iter){
    
    rowvec canTheta = theta.row(iter-1) + stepSize % randn<rowvec>(9);
    double canPrior = priorDensity(canTheta, hyperParams);
    vec canPF = particleFilter(y, canTheta, N, T);
    double ratio = exp(canPrior + canPF[0] - prevPrior - prevPF[0]);
    double u = randu<vec>(1)[0];
    if(u < ratio){
      prevPrior = canPrior;
      prevPF = canPF;
      theta.row(iter) = canTheta;
      a.row(iter) = canPF.col(0).tail(T+1).t();
      s.row(iter) = canPF.col(1).tail(T+1).t();
      d.row(iter) = canPF.col(2).tail(T+1).t();
      v.row(iter) = canPF.col(3).tail(T+1).t();
      x.row(iter) = canPF.col(4).tail(T+1).t();
      y.row(iter) = canPF.col(5).tail(T+1).t();
      accept += 1;
    } else {
      theta.row(iter) = theta.row(iter-1);
      a.row(iter) = a.row(iter-1);
      s.row(iter) = s.row(iter-1);
      d.row(iter) = d.row(iter-1);
      v.row(iter) = v.row(iter-1);
      x.row(iter) = x.row(iter-1);
      y.row(iter) = y.row(iter-1);
    }
    if(iter % 250 == 0){
      Rcpp::Rcout << "Iteration " << iter << std::endl;
    }
  }
  Rcpp::Rcout << accept/reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta.tail_rows(reps),
                            Rcpp::Named("a") = a.tail_rows(reps),
                            Rcpp::Named("s") = s.tail_rows(reps),
                            Rcpp::Named("d") = d.tail_rows(reps),
                            Rcpp::Named("v") = v.tail_rows(reps),
                            Rcpp::Named("x") = x.tail_rows(reps),
                            Rcpp::Named("y") = y.tail_rows(reps));
}