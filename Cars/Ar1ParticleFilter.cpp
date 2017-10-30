// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <boost/math/distributions.hpp> 

using namespace arma;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// [[Rcpp::export]]
mat sobol_points(int N, int D) {
  std::ifstream infile("new-joe-kuo-6.21201" ,ios::in);
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

struct arPF {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  const int lags;
  const int P;
  const double initialV;
  arPF(const mat& dataIn, const vec& epsIn, const vec& hypIn, const int& lagsIn, const int& Pin, const double& initVIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn), lags(lagsIn), P(Pin), initialV(initVIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs; using std::sin; using std::cos;
    int N = data.n_rows;
    double Pi = 3.141593;
    // Create theta as Mu + U %*% Eps
    int dim = 4 + 2 * lags;
    Matrix<T, Dynamic, 1> theta(dim);
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(dim*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1)), sigSqX = exp(theta(2)), sigSqY = exp(theta(3));
    // Constrained (-1, 1)
    T arD = 2.0 / (1 + exp(-theta(4))) - 1, arV = 2.0 / (1 + exp(-theta(5))) - 1;
    // Evaluate log(p(theta))
    T prior =  -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  - 
      (hyperParams(4) + 1) * log(sigSqX)  -  hyperParams(5) / sigSqX  - 
      (hyperParams(6) + 1) * log(sigSqY)  -  hyperParams(7) / sigSqY;
    
    // Evaluate Log |J|  
    T logdetJ = theta(0)  +  theta(1)  +  theta(2)  +   theta(3)  -  theta(4)  -  theta(5) + 2 * log(1 + arV) + 2 * log(1 + arD);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda((dim+1)*i+dim)));
    }
    // Evaluate log likelihood
    T logLik = 0;
    Matrix<T, Dynamic, Dynamic> states(5, P), resampled(5, P);
    Matrix<T, Dynamic, 1>  pi(P), piSum(P), omega(P);
    T omegaSum;
    pi.fill(1.0 / P);
   
    // Sample initial particles
    for(int k = 0; k < P; ++k){
      states(0, k) = sqrt(sigSqV / (1 - pow(arV, 2))) * randn<vec>(1)[0];
      states(1, k) = sqrt(sigSqD / (1 - pow(arD, 2))) * randn<vec>(1)[0];
    }  
    states.row(2).fill(initialV);
    states.row(3).fill(data(0, 0));
    states.row(4).fill(data(0, 1));
    for(int t = 1; t < N; ++t){
      // Create CDF to use for resampling
      piSum(0) = pi(0);
      for(int k = 1; k < P; ++k){
        piSum(k) = pi(k) + piSum(k-1);
      }
      // Resample states via inverse CDF
      for(int k = 0; k < P; ++k){
        double u = randu<vec>(1)[0];
        for(int i = 0; i < P; ++i){
          if(u < piSum(i)){
            for(int j = 0; j < 5; ++j){
              resampled(j, k) = states(j, i);
            }
            break;
          }
        }
      }
      omegaSum = 0;
      // Sample next set of particles and calculate weights
      for(int k = 0; k < P; ++k){
        states(0, k) = arV * resampled(0, k) + sqrt(sigSqV) * randn<vec>(1)[0];           // Accelleration
        states(1, k) = arD * resampled(1, k) + sqrt(sigSqD) * randn<vec>(1)[0];           // Steering Angle
        states(2, k) = resampled(2, k) + states(0, k);                                    // Velocity
        states(3, k) = resampled(3, k) + states(2, k) * cos(Pi/2 + states(1, k));         // Position in Relative X space
        states(4, k) = resampled(4, k) + states(2, k) * sin(Pi/2 + states(1, k));         // Position in Relative Y space
        omega(k) = - log(2*Pi)  -  0.5 * (theta(2) + theta(3))  -                         // Log-likelihood
          pow(states(3, k) - data(t, 0), 2) / (2 * sigSqX)  - 
          pow(states(4, k) - data(t, 1), 2) / (2 * sigSqY);                             
        omegaSum += exp(omega(k));                                                        
      }
      // Normalise weights
      for(int k = 0; k < P; ++k){
        pi(k) = exp(omega(k)) / omegaSum;
      }
      // log(p(x_1:T | theta)) = sum_t log(p(x_t | theta))
      logLik += log(omegaSum / P);
    }
    return prior + logLik + logdetJ;
    }
};

// [[Rcpp::export]]
Rcpp::List arPFDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyper, int lags, int P, double initV){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1>  grad(42);
  int N = data.n_rows;
  // Autodiff
  arPF p(data, epsilon, hyper, lags, P, initV);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

struct arSSM {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  const int lags;
  const double initialV;
  arSSM(const mat& dataIn, const vec& epsIn, const vec& hypIn, const int& lagsIn, const double& initVIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn), lags(lagsIn), initialV(initVIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs; using std::sin; using std::cos;
    int N = data.n_rows;
    // Create theta as Mu + U %*% Eps
    int dim = 4 + 2 * lags;
    Matrix<T, Dynamic, 1> theta(dim);
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(dim*(i+1) + j) * epsilon(j);
      }
    }
    T logdetJ = 0;
    Matrix<T, Dynamic, Dynamic> states(N, 5);
    for(int i = 0; i < N; ++i){
      states(i, 0) = lambda(dim*(dim+1) + 2*i) + lambda(dim*(dim+1) + 2*i + 1) * epsilon(dim + i);
      logdetJ += log(fabs(lambda(dim*(dim+1) + 2*i + 1)));
      states(i, 1) = lambda(dim*(dim+1) + 2*N + 2*i) + lambda(dim*(dim+1) + 2*N + 2*i + 1) * epsilon(N + dim + i);
      logdetJ += log(fabs(lambda(dim*(dim+1) + 2*N + 2*i)));
    }
    states(0, 2) = initialV;
    states(0, 3) = data(0, 0);
    states(0, 4) = data(0, 1);
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1)), sigSqX = exp(theta(2)), sigSqY = exp(theta(3));
    // Constrained (-1, 1)
    T arD = 2.0 / (1 + exp(-theta(4))) - 1, arV = 2.0 / (1 + exp(-theta(5))) - 1;
    // Evaluate log(p(theta))
    T prior =  -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD - 
      (hyperParams(4) + 1) * log(sigSqX)  -  hyperParams(5) / sigSqX  - 
      (hyperParams(6) + 1) * log(sigSqY)  -  hyperParams(7) / sigSqY;
    
    // Evaluate Log |J|  
    logdetJ += theta(0)  +  theta(1)  + theta(2) + theta(3) -  theta(4)  -  theta(5) + 2 * log(1 + arV) + 2 * log(1 + arD);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda((dim+1)*i+dim)));
    }
    // Evaluate log(p(y | x, theta)) + log(p(x | theta))
    T loglik = 0;
    T logstates = 0;
    for(int t = 1; t < N; ++t){
      logstates += - 0.5 * theta(0)  -  pow(states(t, 0) - arV * states(t-1, 0), 2) / (2 * sigSqV) -
        0.5 * theta(1)  -  pow(states(t, 1) - arD * states(t-1, 1), 2) / (2 * sigSqD);
      states(t, 2) = states(t-1, 2) + states(t, 0);
      states(t, 3) = states(t-1, 3) + states(t, 2) * cos(states(t, 1) + 1.570796);
      states(t, 4) = states(t-1, 4) + states(t, 2) * sin(states(t, 1) + 1.570796);
      loglik += - 0.5 * theta(2)  -  pow(states(t, 3) - data(t, 0), 2) / (2 * sigSqX)  -  
        0.5 * theta(3)  -  pow(states(t, 4) - data(t, 1), 2) / (2 * sigSqY);
    }
    return prior + logdetJ + loglik + logstates;
  }
};

// [[Rcpp::export]]
Rcpp::List SSM(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyper, int lags, double initV){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  int N = data.n_rows;
  Matrix<double, Dynamic, 1>  grad(42 + 4 * N); // lags = 1
  // Autodiff
  arSSM p(data, epsilon, hyper, lags, initV);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

