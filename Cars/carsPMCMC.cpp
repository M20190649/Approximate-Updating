// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
vec particleFilter (vec x, rowvec theta, int N, int T, double velocMean, double velocSd, double df, bool Normal){
  // Particles and indices are arranged into T * N matrix
  // Each column is a particle, and each row is a time
  mat vx (T, N);
  // B maps the indices used in the resampling
  // If VX_t^2 = VX_t-1^1, then B(t, 2) = 1.
  // To get a path of a particle, sample column i of B, and then X(t, B(t, i)) tracks particle i at time t. 
  mat B (T, N);
  // At time 0, each particle index points to the original particle
  for(int k = 0; k < N; ++k){
    B(0, k) = k;
  }
  vec pi (N); pi.fill(1.0/N); //normalised weights
  vec piSum(N); // cumulative sum of pi
  vec omega (N); // unnormalised weights
  double sigSqX = theta(0);
  double sigSqE = theta(1);
  double nu = theta(2);
  double delta = theta(3);
  double yDens = 0; // actually the log density
  // x0 from the stationary distribution
  if(Normal){
    vx.row(0) = velocMean  +  velocSd * randn<rowvec>(N);
  } else {
    boost::math::students_t_distribution<> tdist(df);
    vec u = randu<vec>(N);
    for(int i = 0; i < N; ++i){
      vx(0, i) = quantile(tdist, u(i));
    }
  }
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
      // Sample from the skew-t distribution
      double w = randg<vec>(1, distr_param(nu/2, 2.0/nu))[0];
      bool positive = false;
      double z;
      while(!positive){
        z = sqrt(1.0 / w) * randn<vec>(1)[0];
        if(z > 0){
          positive = true;
        }
      }
      double skewT = delta * z  +  sqrt(sigSqE / w) * randn<vec>(1)[0];
      vx(t, k) =  vx(t-1, B(t, k))  +  skewT; // step ahead transition
      omega(k) = 1.0 / sqrt(2*3.14159*sigSqX) * exp(-pow(x(t) - x(t-1) - vx(t, k), 2) / (2*sigSqX)); // x_t ~ N(x_t-1 + 0.01 * vx_t, sigSq)
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
    output(i+1) = vx(i, B(i, index));
  }
  return output;
}

double priorDensity (rowvec theta, vec hyperParams){
  //theta[0] = sigSqX ~ IG(hp[0], hp[1])
  //theta[1] = sigSqV ~ IG(hp[2], hp[3])
  //theta[2] = nu ~ G(hp[4], hp[5])
  //theta[3] = delta ~ N((hp[6], hp[7]))
  double density = - (hyperParams(0) + 1) * log(theta(0))  -   hyperParams(1) / theta(0)  - 
    (hyperParams(2) + 1) * log(theta(1))  -  hyperParams(3) / theta(1)  +  
    (hyperParams(4) - 1) * log(theta(2))  -  hyperParams(5) * theta(2)  -
    pow(theta(3) - hyperParams(6), 2) / (2*hyperParams(7));
  return density;
}

// [[Rcpp::export]]
Rcpp::List PMCMC(vec x, rowvec initialTheta, int reps, int N, vec hyperParams, rowvec stepSize, 
                 double velocMean, double velocSd, double df = 1, bool Normal = true){
  int T = x.n_elem;
  mat theta(reps, 4);
  mat vx(reps, T);
  theta.row(0) = initialTheta;
  double prevPrior = priorDensity(theta.row(0), hyperParams);
  vec prevPF = particleFilter(x, theta.row(0), N, T, velocMean, velocSd, df, Normal);
  vx.row(0) = prevPF.tail(T).t();
  double sigmaSqShape = 1 + (T-1)/2;
  
  
  double accept = 0;
  for(int iter = 1; iter < reps; ++iter){
    //Sigma X Squared 
    //double sigmaSqScale = 1;
    //for(int t = 1; t < T; ++t){
    //    sigmaSqScale = sigmaSqScale + pow(x(t) - x(t-1) - 0.01*vx(iter-1, t), 2) / 2;
    //}
    //rowvec canTheta(4);
    //canTheta(0) = 1.0 / randg<vec>(1, distr_param(sigmaSqShape, 1.0/sigmaSqScale))[0];
    rowvec canTheta = theta.row(iter-1) + stepSize % randn<rowvec>(4);
    if(canTheta(1) < 0.0001){
      canTheta(1) = 0.0001;
    }
    if(canTheta(0) < 0.0001){
      canTheta(0) = 0.0001;
    }
    if(canTheta(2) < 1){
      canTheta(2) = 1;
    }
    double canPrior = priorDensity(canTheta, hyperParams);
    vec canPF = particleFilter(x, canTheta, N, T, velocMean, velocSd, df, Normal);
    double ratio = exp(canPrior + canPF[0] - prevPrior - prevPF[0]);
    double u = randu<vec>(1)[0];
    if(u < ratio){
      prevPrior = canPrior;
      prevPF = canPF;
      theta.row(iter) = canTheta;
      vx.row(iter) = canPF.tail(T).t();
      accept += 1;  
    } else {
      theta.row(iter) = theta.row(iter-1);
      vx.row(iter) = vx.row(iter-1);
    }
    if(iter % 250 == 0){
      Rcpp::Rcout << "Iteration " << iter << std::endl;
    }
  }
  Rcpp::Rcout << accept/(reps-1) << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta,
                            Rcpp::Named("vx") = vx);
}