// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <vector>
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <boost/math/distributions.hpp> 
#include <cstdlib> 
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;
using namespace arma;

// evaluate and automatially differentiate log(p(y_1:T, x_T+1 theta)) wrt theta and x_T+1
struct logP {
  const vec& x; // data
  const int& obs; // T is taken for the variable type, so obs = number of elements in x
  const int& N; // Number of particles
  const vec& hyperParams;
  logP(const vec& xIn, const int& obsIn, const int& nIn, const vec& hyperIn) : // constructor
    x(xIn), obs(obsIn), N(nIn), hyperParams(hyperIn) {}
  template <typename T> // T will become stan's var that can be automatically differentiated - any variable differentiated in the chain rule for d logp / d theta must be type T
  T operator ()(const Matrix<T, Dynamic, 1>& theta) // derivative is with respect to theta (with xT+1 included), so theta matrix is the operator () input.
    const{
    using std::log; using std::pow; using std::exp; using std::sqrt; using std::tgamma; using std::abs;
    // explicitly declare mathematical operators used so stan knows to differentiate these.
    
    // Get elements of theta
    T sigSqX = theta(0);
    T sigSqE = theta(1);
    T nu = theta(2);
    T phi = theta(3);
    T VTplus1 = theta(4);
    // bottom T*N rows of theta are w, a function of nu, w ~ Gamma(nu/2, nu/2)
    // but w = f(nu) and dw/dnu are not analytically available. 
    // We will get the direct effect of nu in the prior, and the indirect effect of nu through w in the likelihood.
    // Combine these two outside of the auto-diff to get the full d log P / d nu

    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqX)  -  hyperParams(1) / sigSqX  - 
      (hyperParams(2) + 1) * log(sigSqE)  -  hyperParams(3) / sigSqE  +  
      (hyperParams(4) - 1) * log(nu)  -  hyperParams(5) * nu  -
      pow(phi - hyperParams(6), 2) / (2 * hyperParams(7));

    // Set up for particle filtering
    Matrix<T, Dynamic, 1> vOld(N), vNew(N), vResample(N), pi(N), piSum(N), omega(N);
    pi.fill(1.0/N); // initial weights
    T dxAllDens, dxtDens, omegaSum, VTplus1Dens;
    
    // Sample v0 from stationary distribution - switch to students t for real data
    for(int k = 0; k < N; ++k){
      vOld(k) = sqrt(sigSqE / (1 - pow(phi, 2))) * randn<vec>(1)[0];
    }
    dxAllDens = 0;
    // Particle Filter loop
    for(int t = 1; t < obs; ++t){
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
            vResample(k) = vOld(i);
            break;
          }
        }
      }
      // Calculate weights
      omegaSum = 0;
      for(int k = 0; k < N; ++k){
        vNew(k) = phi * vResample(k)  +  sqrt(sigSqE / theta(5+t*N+k)) * randn<vec>(1)[0];; // Step Ahead transition
        omega(k) = 1.0 / sqrt(2*3.14159*sigSqX) * exp(-pow(x(t) - x(t-1) - vNew(k), 2) / (2*sigSqX)); // Measurement density
        omegaSum += omega(k); // sum of weights for normalisation
      }
      // Normalise weights
      for(int k = 0; k < N; ++k){
        pi(k) = omega(k) / omegaSum;
      }
      // log(p(x_1:T | theta)) = sum_t log(p(x_t | theta))
      dxtDens = log(omegaSum / N);
      dxAllDens += dxtDens;
      // reset vOld as we step from t to t+1
      for(int k = 0; k < N; ++k){
        vOld(k) = vNew(k);
      }
    } // end of particle filter loop
    VTplus1Dens = 0; 
    // p(V_T+1 | theta, x_1:T) approx sum_k pi_k p(V_T+1 | V_T,k, theta)
    for(int k = 0; k < N; ++k){
      VTplus1Dens += pi(k) / sqrt(sigSqE) * tgamma((nu+1)/2) / (tgamma(nu/2) * sqrt(nu * 3.14159)) * 
        pow(1 + pow(VTplus1 - phi*vNew(k), 2)/(sigSqE*nu), -(nu+1)/2);
    }
    // interested in log density
    VTplus1Dens = log(VTplus1Dens);
    return dxAllDens + VTplus1Dens + prior;
  }
};


// evaluate and automatially differentiate logQ
struct logQ {
  const MatrixXd& epsilon;
  logQ(const MatrixXd& epsIn) : // constructor
    epsilon(epsIn){}
  template <typename T> // 
  T operator ()(const Matrix<T, Dynamic, 1>& lambda) // derivative is with respect to lambda
    const{
    using std::log; using std::fabs; using std::pow;

    T logEpsDens = 0;
    T logDetJ = 0;
    for(int i = 0; i < 3; ++i){
      logDetJ += lambda(i);
      for(int j = 0; j <= i; ++j){
        logDetJ += lambda(5*(i+1)+j) * epsilon(j);
      }
    }
    for(int i = 0; i < 5; ++i){
      logEpsDens -= 0.5 * pow(epsilon(i), 2);
      logDetJ +=  log(fabs(lambda((i+1)*5+i)));
    }
    return logEpsDens - logDetJ;
  }
};

// d logq / d lambda = d log(det(J)) / d lambda
Matrix<double, Dynamic, 1> QDeriv(Matrix<double, Dynamic, 1> epsilon, Matrix<double, Dynamic, 1>  lambda){
  Matrix<double, Dynamic, 1> output(30); 
  output.fill(0);
  for(int i = 0; i < 3; ++i){
    output(i) = -1;
    for(int j = 0; j <=i; ++j){
      output((i+1)*5+j) = - epsilon(j);
    }
  }
  for(int i = 0; i < 5; ++i){
    output(6*i + 5) -= 1.0 / lambda(6*i + 5);
  }
  return output;
}

// d log theta/xTplus1 / d lambda - Really need to make this better, getting a bit out of hand now.
Matrix<double, Dynamic, 1> FDeriv (Matrix<double, Dynamic, 1> epsilon, Matrix<double, Dynamic, 1>  lambda){
  Matrix<double, Dynamic, 1> output(30); output.fill(0);
  // exponential link parameters
  for(int i = 0; i < 3; ++i){
    output(i) = lambda(i);
    for(int j = 0; j <= i; ++j){
      output(i) += lambda(5*(i+1) + j) * epsilon(j); 
    }
    output(i) = exp(output(i));
    for(int j = 0; j <= i; ++j){
      output(5*(i+1) + j) = output(i) * epsilon(j);
    }
  }
  // identity link parameters
  for(int i = 3; i < 5; ++i){
    output(i) = 1;
    for(int j = 0; j <= i; ++j){
      output(5*(i+1) + j) = epsilon(j);
    }
  }
 
  return output;
}

// d w / d nu
vec WDeriv (Matrix<double, Dynamic, 1> w, vec u, double nu){
  int N = u.n_elem;
  double h = 0.000001;
  vec output(N);
  boost::math::gamma_distribution<> wDist((nu+h)/2, (nu+h)/2);
  for(int i = 0; i < N; ++i){
    double wh = quantile(wDist, u(i));
    output(i) = (wh - w(i))/h;
  }
  return output;
}

// generate sobol points
mat sobol_points(int N, int D) {
  ifstream infile("new-joe-kuo-6.21201" ,ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 
  
  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [N];
  C[0] = 1;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  mat POINTS(N, D, fill::zeros);
  
  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1
  
  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [N];
  X[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS(i, 0) = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  
  // Clean up
  delete [] V;
  delete [] X;
  
  
  // ----- Compute the remaining dimensions -----
  for (unsigned j=1;j<=D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
        V[i] = V[i-s] ^ (V[i-s] >> s); 
        for (unsigned k=1;k<=s-1;k++) 
          V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
      POINTS(i, j) = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
    }
    
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
  
  return POINTS;
}

// randomly shuffle sobol points
mat shuffle(mat sobol){
  int N = sobol.n_rows;
  int D = sobol.n_cols;
  mat output(N, D, fill::zeros);
  // draw a random rule of: switch 1 and 0  /  do not switch for each binary digit.
  vec rule = randu<vec>(16);
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < D; ++j){
      // grab element of the sobol sequence
      double x = sobol(i, j);
      // convert to a binary representation
      uvec binary(16, fill::zeros);
      for(int k = 1; k < 17; ++k){
        if(x > pow(2, -k)){
          binary(k-1) = 1;
          x -= pow(2, -k);
        }
      }
      // apply the transform of tilde(x_k) = x_k + a_k mod 2, where a_k = 1 if rule_k > 0.5, 0 otherwise
      for(int k = 0; k < 16; ++k){
        if(rule(k) > 0.5){
          binary(k) = (binary(k) + 1) % 2;
        }
      }
      // reconstruct base 10 number from binary representation
      for(int k = 0; k < 16; ++k){
        if(binary(k) == 1){
          output(i, j) += pow(2, -(k+1));
        }
      }
    }
  }
  return output;
}

// Main VB algorithm
// [[Rcpp::export]]
Rcpp::List VB_Cars (vec x, Rcpp::NumericMatrix lambdaIn, vec hyperParams, int S = 1, int N = 100, int maxIter = 5000,
                    double alpha = 0.1, double threshold = 0.01, double thresholdIS = 0){
  
  // convert to Eigen format
  Map<MatrixXd> lambdaMat(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  // Reshape matrix as a vector
  Map<MatrixXd> lambda(lambdaMat.data(), 30, 1);
  // Initialise various components for the algorithm
  int T = x.n_elem;
  Matrix<double, Dynamic, 1> gradientQ(30), gradientF(30), meanGradient(30), Mt(30), Vt(30), LB(maxIter+1), logPeval(S), gradP(5 + N*T);
  MatrixXd theta(5, S), gradientP(30, S), epsilon(5, S);
  LB.fill(0); Mt.fill(0); Vt.fill(0);
  double logQeval;
  // Initial SOBOL QMC numbers. first few numbers in each sequence are too similar so we want to work from 101'th to 100+S'th number in each sequence
  mat sobol = sobol_points(S+100, 5);
  // store quasi random "uniform" [0, 1] numbers that result from shuffling the bits of the sobol numbers
  mat unif;
  // Standard normal density to transform [0,1] numbers to standard normal
  boost::math::normal_distribution<> epsDist(0, 1);
  // Loop control
  int iter = 0;
  double diff = threshold + 1;
  double meanLB = 0;
  double omeanLB;
  
  while(diff > threshold){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // Reset derivatives to zero
    meanGradient.fill(0);
    Vt.fill(0);

    // Importance Sample check
    double checkIS = randu<vec>(1)[0];
    // do not IS on first iteration! 
    if(iter == 1 | checkIS > thresholdIS){
      // Did not importance sample: Draw new epsilon/theta, calculate new dlogp / dtheta and other gradients
      
      // Randomly shuffle the sobol numbers for RQMC
      unif = shuffle(sobol);
      // Create estimates of E(dLB/dlam) and E(dLB/dlam^2)
      for(int s = 0; s < S; ++s){
        // transform uniform numbers to standard normal
        for(int i = 0; i < 5; ++i){
          // Quantile function is unstable for extreme values of unif
          if(unif(s+100, i) > 0.9995){
            epsilon(i, s) = 3.290527;
          } else if (unif(s+100, i) < 0.0005) {
            epsilon(i, s) = -3.290527;
          } else {
            epsilon(i, s) = quantile(epsDist, unif(s+100, i));
          }
        }
        // create thetas, stored for importance sampler
        for(int i = 0; i < 5; ++i){
          theta(i, s) = lambda(i);
          for(int j = 0; j <= i; ++j){
            theta(i, s) += lambda((i+1)*5+j) * epsilon(j, s);
          }
        }
        theta(0, s) = exp(theta(0, s)); // log(sigSq) -> sigSq
        theta(1, s) = exp(theta(1, s));
        theta(2, s) = exp(theta(2, s)); // log(nu) -> nu
        if(theta(3, s) >= 1){ // keep phi stationary
          theta(3, s) = 0.99;
        }
        // generate w
        vec u = randu<vec>(T*N);
        boost::math::gamma_distribution<> wDist(theta(2, s)/2, theta(2, s)/2);
        Matrix<double, Dynamic, 1>  w(T*N);
        for(int i = 0; i < T*N; ++i){
          if(u(i) > 0.9995){
            u(i) = 0.9995;
          } else if(u(i) < 0.0005) {
            u(i) = 0.0005;
          }
          w(i) = quantile(wDist, u(i));
        }
        
        // attach w to theta.
        Matrix<double, Dynamic, 1>  parameters(5 + N*T);
        parameters.topRows(5) = theta.col(s);
        parameters.bottomRows(N*T) = w;

        // Autodiff
        logP p(x, T, N, hyperParams);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(p, parameters, logPeval(s), gradP);
        
        Matrix<double, Dynamic, 1> eps = epsilon.col(s);
        logQ q(eps);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(q, lambda, logQeval, gradientQ);

        
        //  deriv nu = dL/dnu + sum dw/dnu * dw/dL
        vec GradientW = WDeriv(w, u, theta(2, s));
        for(int i = 0; i < N*T; ++i){
          gradP(2) += gradP(5 + i) * GradientW(i);
        }

        // take derivative of theta
        gradientF = FDeriv(epsilon.col(s), lambda);

        // update estimates of deriv(logp - logq) and (deriv(logp - logq))^2
        for(int i = 0; i < 30; ++i){
          // grab relevant element of dlogp / dtheta
          if(i < 5){
            gradientP(i, s) = gradP(i);
          } else {
            gradientP(i, s) = gradP(floor((i-5)/5));
          }
          //update gradients
          meanGradient(i) += (gradientP(i, s) * gradientF(i) - gradientQ(i)) / S;
          //meanGradientSq(i) += pow(S * meanGradient(i), 2) / S;
        }
        LB(iter-1) += (logPeval(s) - logQeval)/S;

      }
    } else {
      // Importance Sampling: reuse derivatives and simulated values from most recent non IS iteration
      
      for(int s = 0; s < S; ++s){
        // Transform theta back to the epsilon that would have generated them with current lambda values
        Matrix<double, Dynamic, 1>  impliedEpsilon(5);
        impliedEpsilon(0) = (log(theta(0, s)) - lambda(0)) / lambda(5); 
        impliedEpsilon(1) = (log(theta(1, s)) - lambda(1) - lambda(10) * impliedEpsilon(0)) / lambda(11);
        impliedEpsilon(2) = (log(theta(2, s)) - lambda(2) - lambda(15) * impliedEpsilon(0) - lambda(16) * impliedEpsilon(1)) / lambda(17);
        for(int i = 3; i < 5; ++i){
          impliedEpsilon(i) = theta(i, s) - lambda(i);
          for(int j = 0; j < i; ++j){
            impliedEpsilon(i) -= lambda((i+1)*5 + j) * impliedEpsilon(j);
          }
          impliedEpsilon(i) = impliedEpsilon(i) / lambda((i+1)*5 + i);
        }
        // Take derivatives of Q and F with these new values
        logQ q(impliedEpsilon);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(q, lambda, logQeval, gradientQ);
        gradientF = FDeriv(impliedEpsilon, lambda);
        // Calculate importance sample weight as q(impEps)/q(eps), q ~ MVN(0, I)
        double weight = 0;
        for(int i = 0; i < 5; ++i){
          weight += pow(epsilon(i, s), 2) - pow(impliedEpsilon(i), 2);
        }
        weight = exp(weight/2);
        // Calculate derivative estimates
        for(int i = 0; i < 30; ++i){
          meanGradient(i) += (gradientP(i, s) * gradientF(i) - gradientQ(i)) * weight / S;
        }
        // Update ELBO
        LB(iter-1) += (logPeval(s) - logQeval) / S;
      }
    }
    // adagrad updating rule for lambda values
    for(int i = 0; i < 30; ++i){
      Mt(i) += pow(meanGradient(i), 2);
      if(Mt(i) != 0){
        Vt(i) = pow(Mt(i), -0.5);
      }
      if(iter > 1){
        lambda(i) += alpha * Vt(i) * meanGradient(i);
      }
    }
    // check convergence by looking at the difference of the average of past five LB values and the average of the five before that
    if(iter % 5 == 0){
      omeanLB = meanLB;
      meanLB = 0.2 * (LB(iter-1) + LB(iter-2) + LB(iter-3) + LB(iter-4) + LB(iter-5));
      diff = std::fabs(meanLB - omeanLB);
    } 
    // report progress from time to time
    if(iter % 25 == 0){
      Rcpp::Rcout << "Iteration: " << iter << ", ELBO: " << meanLB << std::endl;
    }
  } // End while loop
  if(iter <= maxIter & iter > 2){ // did converge in time
    // iter goes up by one before checking the maxIter condition, so need to use LB(iter-2)
    Rcpp::Rcout << "Converged after " << iter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  } else if(iter < maxIter){ // reached maxIter, so did not converge
    Rcpp::Rcout << "Warning, failed to converge after " << maxIter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  }
  // Grab elements of lambda as the more meaningful approximation mean vector and upper triangular cholesky factor of the variance matrix
  MatrixXd Mu = lambda.topRows(5);
  Map<MatrixXd> U(lambda.bottomRows(25).data(), 5, 5);
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("U") = U,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}

