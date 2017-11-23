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

struct arUpdate {
  const mat data;
  const vec epsilon;
  const vec mean;
  const mat Linv;
  const int lags;
  arUpdate(const mat& dataIn, const vec& epsIn, const vec& meanIn, const mat& LinvIn, const int& lagsIn) :
    data(dataIn), epsilon(epsIn), mean(meanIn), Linv(LinvIn), lags(lagsIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = data.n_rows;
    int dim = 2 + 2 * lags;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(dim);
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(dim*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Evaluate log(p(theta))
    T prior = 0;
    Matrix<T, Dynamic, 1> kernel(dim);
    kernel.fill(0);
    for(int i = 0; i < dim; ++i){
      for(int j = 0; j <= i; ++j){
        kernel(i) += (theta(j) - mean(j)) * Linv(i, j);
      }
      prior += - 0.5 * pow(kernel(i), 2);
    }
    
    // Evaluate Log Det J
    T logdetJ = theta(0)  +  theta(1);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda((dim+1)*i+dim)));
    }
    // Evaluate log likelihood
    T logLik = 0;
    Matrix<T, Dynamic, Dynamic> loglikKernel(N, 2);
    for(int t = lags; t < N; ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lags; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(1 + i);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(1 + lags + i);
      }
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List arUpdater(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec mean, mat Linv, int lags){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1>  grad(6 + 10 * lags + 4 * lags * lags);
  // Autodiff
  arUpdate p(data, epsilon, mean, Linv, lags);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// [[Rcpp::export]]
cube evalFcDens (mat data, mat means, mat L, int M, int S, vec asup, vec dsup, Rcpp::List MCMC){
  int N = asup.n_elem;
  int K = means.n_cols;
  cube densities (N * S, 4, K, fill::zeros);
  boost::math::normal_distribution<> Zdist(0, 1);
  for(int k = 0; k < K; ++k){
    mat MCMCdraws = MCMC(k);
    for(int m = 0; m < M; ++m){
      // VB forecast densities
      double afc1 = data(1, 0), afc2 = data(0, 0), dfc1 = data(1, 1), dfc2 = data(0, 1), afc, dfc;
      vec draws = means.col(k) + L.rows(k*6, (k+1)*6-1) * randn<vec>(6);
      for(int h = 0; h < S; ++h){
        afc = afc1 * draws(2) + afc2 * draws(3);
        dfc = dfc1 * draws(4) + dfc2 * draws(5);
        for(int n = 0; n < N; ++n){
          densities(N*h + n, 0, k) += pdf(Zdist, (asup(n) - afc) / sqrt(std::exp(draws(0)))) / M;
          densities(N*h + n, 1, k) += pdf(Zdist, (dsup(n) - dfc) / sqrt(std::exp(draws(1)))) / M;
        }
        afc2 = afc1;
        afc1 = afc;
        dfc2 = dfc1;
        dfc1 = dfc;
      }
      // MCMC forecast densities
      afc1 = data(1, 0); afc2 = data(0, 0); dfc1 = data(1, 1); dfc2 = data(0, 1);
      draws = MCMCdraws.row(1003 + 4 * m).t();
      for(int h = 0; h < S; ++h){
        afc = afc1 * draws(2) + afc2 * draws(3);
        dfc = dfc1 * draws(4) + dfc2 * draws(5);
        for(int n = 0; n < N; ++n){
          densities(N*h + n, 2, k) += pdf(Zdist, (asup(n) - afc) / sqrt(std::exp(draws(0)))) / M;
          densities(N*h + n, 3, k) += pdf(Zdist, (dsup(n) - dfc) / sqrt(std::exp(draws(1)))) / M;
        }
        afc2 = afc1;
        afc1 = afc;
        dfc2 = dfc1;
        dfc1 = dfc;
      }
    }
  }
  return densities;
}

// [[Rcpp::export]]
double nlogDensity (mat data, vec theta, vec mu, mat varInv){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  for(int t = 2; t < T; ++t){
    dens +=  - 0.5 * (theta(0) + theta(1)) - 
      pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (2 * exp(theta(0))) -
      pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (2 * exp(theta(1)));
  }
  return dens;
}

// [[Rcpp::export]]
double tlogDensity (mat data, vec theta, vec mu, mat varInv){
  int T = data.n_rows;
  double sigSqV = exp(theta(0)), sigSqD = exp(theta(1)), nuV = exp(theta(6)), nuD = exp(theta(7)),
    dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  for(int t = 2; t < T; ++t){
    dens += lgamma((1+nuV)/2)  +  lgamma((1+nuD)/2)  -  lgamma(nuV/2)  -  lgamma(nuD/2) -
      0.5 * (theta(0) + theta(1) + theta(6) + theta(7))  - 
      (nuV+1)/2 * log(1 + pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (nuV * sigSqV))  -
      (nuD+1)/2 * log(1 + pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (nuD * sigSqD));
  }
  return dens;
}