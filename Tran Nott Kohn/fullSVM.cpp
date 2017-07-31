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

// evaluate and automatially differentiate log(p(y_1:T, x_T+1, theta)) wrt theta and x_T+1
struct logP {
  const MatrixXd& y; // data
  const int& obs; // T is taken for the variable type, so obs = number of elements in y
  logP(const MatrixXd& yIn, const int& obsIn) : // constructor
    y(yIn), obs(obsIn) {}
  template <typename T> // T will become stan's var that can be automatically differentiated - any variable differentiated in the chain rule for d logp / d theta must be type T
  T operator ()(const Matrix<T, Dynamic, 1>& theta) // derivative is with respect to theta (with xT+1 included), so theta matrix is the operator () input.
    const{
    using std::log; using std::pow; using std::exp; using std::sqrt; // explicitly declare mathematical operators used so stan knows to differentiate these.
    
    // Get elements of theta
    T sigSq = theta(0);
    T mu = theta(1);
    T tau = theta(2);
    T phi = 2*tau - 1;
    Matrix<T, Dynamic, 1> x = theta.bottomRows(obs+2);
    // Evaluate log(p(theta)) + log(p(x_0 | theta))
    T prior = -pow(mu, 2) / 20  +  19 * log(tau)  + 0.5 * log(1-tau) - 3.5 * log(sigSq)  -  0.025 / sigSq;
    T states = - 0.5 * log(sigSq)  +  0.5 * log(1 - pow(phi, 2))  -  pow(theta(3) - mu, 2) * (1 - pow(phi, 2)) / (2*sigSq);
    T data = 0;
    // Run through sum_t log(p(y_t | x_t, theta) * p(x_t | x_t-1, theta))
    for(int t = 1; t <= obs; ++t){
      states +=  - 0.5 * log(sigSq)  -   pow(theta(3+t) - mu - phi*(theta(2+t)-mu), 2) / (2*sigSq);
      data += - 0.5 * theta(3+t)  -  pow(y(t-1), 2) / (2*exp(theta(3+t)));
    }
    // log(p(x_T+1 | x_T, theta))
    states += - 0.5 * log(sigSq)  -  pow(theta(obs+4) - mu - phi * (theta(obs+3)-mu), 2) / (2*sigSq);
    
    return prior + states + data;
  }
};

// evaluate log(q(theta, xT+1)) = log(q(eps)) - log(det(jacobian))
double evalLogQ (Matrix<double, Dynamic, 1> epsilon, MatrixXd lambda, int T){
  double logEpsDens = 0;
  double logDetJ = lambda(0) + lambda(4)*epsilon(0);
  for(int i = 0; i < T+5; ++i){
    logEpsDens -= 0.5 * pow(epsilon(i), 2);
    logDetJ +=  log(fabs(lambda(T+5 + i)));
  }
  return logEpsDens - logDetJ;
}

// d logq / d lambda = d log(det(J)) / d lambda
Matrix<double, Dynamic, 1> QDeriv(Matrix<double, Dynamic, 1> epsilon, MatrixXd lambda, int T){
  Matrix<double, Dynamic, 1> output(2*T+10); output.fill(0); 
  output.fill(0);
  output(0) = - 1; //mu_1
  output(T+5) = - (epsilon(0) + 1.0/lambda(T+5)); //U_11
  for(int t = T+6; t < 2*T+10; ++t){
    output(t) = - 1.0 / lambda(t); //U_ii
  }
  return output;
}

// d log theta/xT / d lambda
Matrix<double, Dynamic, 1> FDeriv (Matrix<double, Dynamic, 1> epsilon, MatrixXd lambda, int T){
  Matrix<double, Dynamic, 1> output(2*T+10); output.fill(1);
  output(0) = exp(lambda(0) + lambda(T+5) * epsilon(0));
  output(T+5) = epsilon(0) * output(0);
  for(int t = T+6; t < 2*T + 10; ++t){
    output(t) = epsilon(t-T+5);
  }
  return output;
}

// generate sobol points
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
Rcpp::List ADVI (Rcpp::NumericMatrix yIn, Rcpp::NumericMatrix lambdaIn, int S, int maxIter = 5000,
                    double alpha = 0.1, double beta1 = 0.9, double beta2 = 0.999, double threshold = 0.01, double thresholdIS = 0){
  
  // convert to Eigen format
  Map<MatrixXd> y(Rcpp::as<Map<MatrixXd> >(yIn));
  Map<MatrixXd> lambdaMat(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  int T = y.rows();
  // Reshape matrix as a vector
  Map<MatrixXd> lambda(lambdaMat.data(), 2*T+10, 1);
  // Initialise various components for the algorithm
  Matrix<double, Dynamic, 1> gradientQ(2*T+10), gradientF(2*T+10), meanGradient(2*T+10), meanGradientSq(2*T+10), 
  Mt(2*T+10), Vt(2*T+10), LB(maxIter+1), logPeval(S), gradP(T+5);
  MatrixXd theta(T+5, S), gradientP(2*T+10, S), epsilon(T+5, S);
  LB.fill(0); Mt.fill(0); Vt.fill(0);
  double logQeval;
  
  
  // Initial SOBOL QMC numbers. first few numbers in each sequence are too similar so we want to work from 101'th to 100+S'th number in each sequence
  mat sobol = sobol_points(S+100, T+5);
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
    meanGradientSq.fill(0);
    
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
        for(int i = 0; i < T+5; ++i){
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
        for(int i = 0; i < T+5; ++i){
          theta(i, s) = lambda(i) + lambda(T+5+i) * epsilon(i, s);
        }
        theta(0, s) = exp(theta(0, s)); // log(sigSq) -> sigSq
        // a bad way to force stationarity
        if(theta(2, s) > 0.99) {
          theta(2, s) = 0.99;
        } else if(theta(2, s) < 0.01){
          theta(2, s) = 0.01;
        }
        if(theta(0, s) > 0.3){
          theta(0, s) = 0.3;
        }
        // take derivative of log joint
        logP p(y, T);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(p, theta.col(s), logPeval(s), gradP);
        // take derivative of q
        logQeval = evalLogQ(epsilon.col(s), lambda, T);
        gradientQ = QDeriv(epsilon.col(s), lambda, T);
        // take derivative of theta
        gradientF = FDeriv(epsilon.col(s), lambda, T);
        // update estimates of deriv(logp - logq) and (deriv(logp - logq))^2
        for(int i = 0; i < 2*T+5; ++i){
          // grab relevant element of dlogp / dtheta
          if(i < T+5){
            gradientP(i, s) = gradP(i);
          } else {
            gradientP(i, s) = gradP(i - T+5);
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
        Matrix<double, Dynamic, 1> impliedEpsilon(T+5);
        for(int i = 0; i < T+5; ++i){
          if(i == 0){
            impliedEpsilon(i) = (log(theta(i, s)) - lambda(i)) / lambda(T+5 + i);
          } else {
            impliedEpsilon(i) = (theta(i, s) - lambda(i)) / lambda(T+5 + i);
          }
        }
        // Take derivatives of Q and F with these new values
        logQeval = evalLogQ(impliedEpsilon, lambda, T);
        gradientQ = QDeriv(impliedEpsilon, lambda, T);
        gradientF = FDeriv(impliedEpsilon, lambda, T);
        // Calculate importance sample weight as q(impEps)/q(eps), q ~ MVN(0, I)
        double weight = 0;
        for(int i = 0; i < T+5; ++i){
          weight += pow(epsilon(i, s), 2) - pow(impliedEpsilon(i), 2);
        }
        weight = exp(weight/2);
        // Calculate derivative estimates
        for(int i = 0; i < 2*T+10; ++i){
          meanGradient(i) += (gradientP(i, s) * gradientF(i) - gradientQ(i)) * weight / S;
          // meanGradientSq(i) += pow(S * meanGradient(i), 2) / S;
        }
        // Update ELBO
        LB(iter-1) += (logPeval(s) - logQeval) / S;
      }
    }
    // adam updating rule for lambda values
    for(int i = 0; i < 2*T+10; ++i){
      //Mt(i) = Mt(i) * beta1  +  (1 - beta1) * meanGradient(i);
      //Vt(i) = Vt(i) * beta2  +  (1 - beta2) * meanGradientSq(i);
      //double MtHat = Mt(i) / (1 - pow(beta1, iter));
      //double VtHat = Vt(i) / (1 - pow(beta2, iter));
      //lambda(i) += alpha * MtHat / (sqrt(VtHat) + pow(1, -8));
      //lambda(i) += alpha / (1 + iter) * meanGradient(i);
      Mt(i) += pow(meanGradient(i), 2);
      if(Mt(i) != 0){
        Vt(i) = pow(Mt(i), -0.5);
      } 
      if(iter > 1){
        lambda(i) += alpha * Vt(i) * meanGradient(i);
      }
    }
    if(lambda(2) > 0.99) {
      lambda(2) = 0.99;
    } else if(lambda(2) < 0.01) {
      lambda(2) = 0.01;
    }
    if(lambda(0) > 0) {
      lambda(0) = 0;
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
  if(iter <= maxIter){ // did converge in time
    // iter goes up by one before checking the maxIter condition, so need to use LB(iter-2)
    Rcpp::Rcout << "Converged after " << iter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  } else { // reached maxIter, so did not converge
    Rcpp::Rcout << "Warning, failed to converge after " << maxIter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  }
  // Grab elements of lambda as the more meaningful approximation mean vector and upper triangular cholesky factor of the variance matrix
  MatrixXd Mu = lambda.topRows(T+5);
  MatrixXd Sd = lambda.bottomRows(T+5);
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("Sd") = Sd,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}


