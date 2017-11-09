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

struct ar {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  const int lags;
  ar(const mat& dataIn, const vec& epsIn, const vec& hypIn, const int& lagsIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn), lags(lagsIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = data.n_rows;
    // Create theta as Mu + U %*% Eps
    int dim = 2 + 2 * lags;
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
    T prior =  -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD;
    for(int i = 0; i < 2 * lags; ++i){
      prior -= pow(theta(i) - hyperParams(4 + 2 * i), 2) / (2 * hyperParams(5 + 2 * i));
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
        loglikKernel(t, 0) -= data(t-i, 0) * theta(2*i);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(2*i+1);
      }
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List arDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperParams, int lags){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1> grad(6 + 10 * lags + 4 * lags * lags);
  // Autodiff
  ar p(data, epsilon, hyperParams, lags);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

struct heirAr {
  const mat data;
  const vec epsilon;
  const vec hyperMean;
  const mat hyperLinv;
  const mat Linv;
  const vec obsSum;
  const int lags;
  const bool diag;
  heirAr(const mat& dataIn, const vec& epsIn, const vec& hypMIn, const mat& hypLiIn, const mat& LinvIn, const vec& obsIn, const int& lagsIn, const bool& diagIn) :
    data(dataIn), epsilon(epsIn), hyperMean(hypMIn), hyperLinv(hypLiIn), Linv(LinvIn), obsSum(obsIn), lags(lagsIn), diag(diagIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = obsSum.n_elem;
    // Create thetaHat as Mu + U %*% Eps
    int dim = 2 + 2 * lags;
    Matrix<T, Dynamic, 1> thetaHat(dim);
    T logdetJ = 0;
    int lCounter = dim*(N+1);
    for(int i = 0; i < dim; ++i){
      thetaHat(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        if(i == j){
          thetaHat(i) += lambda(lCounter) * epsilon(i);
          logdetJ += log(fabs(lambda(lCounter)));
        }  else if(!diag){
          thetaHat(i) += lambda(lCounter) * epsilon(j);
        }
        lCounter += 1;
      }
    }
    Matrix<T, Dynamic, Dynamic> theta(dim, N);
    for(int k = 0; k < N; ++k){
      for(int i = 0; i < dim; ++i){
        theta(i, k) = lambda(dim*(k+1) + i);
        for(int j = 0; j <= i; ++j){
          if(i == j){
            theta(i, k) += lambda(lCounter) * epsilon(dim*(k+1) + i);
            logdetJ += log(fabs(lambda(lCounter)));
          } else if(!diag){
            theta(i, k) += lambda(lCounter) * epsilon(dim*(k+1) + j);
          }
          lCounter += 1;
        }
      }
    }
    Matrix<T, Dynamic, 1> sigSqV(N), sigSqD(N);
    for(int k = 0; k < N; ++k){
      sigSqV(k) = exp(theta(0, k));
      sigSqD(k) = exp(theta(1, k));
    }
  
    // Finish evaluation of Log Det J
    for(int k = 0; k < N; ++k){
      logdetJ += theta(0, k) + theta(1, k);
    }
    // Evaluate log(p(theta hat))
    T hyperprior = 0;
    Matrix<T, Dynamic, 1> hypKernel(dim);
    hypKernel.fill(0);
    for(int i = 0; i < dim; ++i){
      for(int j = i; j < dim; ++j){
        hypKernel(i) += (thetaHat(j) - hyperMean(j)) * hyperLinv(j, i);
      }
      hyperprior += - 0.5 * pow(hypKernel(i), 2);
    }
    
    // Evaluate log(p(theta | theta hat))
    T prior = 0;
    Matrix<T, Dynamic, Dynamic> kernel(dim, N);
    kernel.fill(0);
    for(int k = 0; k < N; ++k){
      for(int i = 0; i < dim; ++i){
        for(int j = i; j < dim; ++j){
          kernel(i, k) += (theta(j , k) - thetaHat(j)) * Linv(j, i);
        }
        prior += - 0.5 * pow(kernel(i, k), 2);
      }
    }

    // Evaluate log likelihood
    T logLik = 0;
    Matrix<T, Dynamic, Dynamic> loglikKernel(data.n_rows, 2);
    for(int t = lags; t < obsSum(0); ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lags; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(2*i, 0);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(2*i+1, 0);
      }
      logLik += - 0.5 * (theta(0, 0) + theta(1, 0)) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV(0)) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD(0));
    }
    for(int k = 1; k < N; ++k){
      for(int t = obsSum(k-1) + lags; t < obsSum(k); ++t){
        loglikKernel(t, 0) = data(t, 0);
        loglikKernel(t, 1) = data(t, 1);
        for(int i = 1; i <= lags; ++i){
          loglikKernel(t, 0) -= data(t-i, 0) * theta(2*i, k);
          loglikKernel(t, 1) -= data(t-i, 1) * theta(2*i+1, k);
        }
        logLik += - 0.5 * (theta(0, k) + theta(1, k)) - 
          pow(loglikKernel(t, 0), 2) / (2 * sigSqV(k)) -
          pow(loglikKernel(t, 1), 2) / (2 * sigSqD(k));
      }
    }
    return hyperprior + prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List heirArDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperMean, mat hyperLinv, mat Linv, vec obsSum, int lags, bool diag = false){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  int N = obsSum.n_elem;
  int dimT = 2 + 2 * lags;
  int dimU = 0;
  for(int i = 0; i < dimT; ++i){
    dimU += i;
  }
  Matrix<double, Dynamic, 1> grad((dimT + dimU) * (N + 1));
  // Autodiff
  heirAr p(data, epsilon, hyperMean, hyperLinv, Linv, obsSum, lags, diag);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
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
      for(int j = i; j < dim; ++j){
        kernel(i) += (theta(j) - mean(j)) * Linv(j, i);
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
        loglikKernel(t, 0) -= data(t-i, 0) * theta(2*i);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(2*i+1);
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


struct arUpdateMix {
  const mat data;
  const vec epsilon;
  const vec mean;
  const mat Linv;
  const int lags;
  const vec priorComp;
  arUpdateMix(const mat& dataIn, const vec& epsIn, const vec& meanIn, const mat& LinvIn, const int& lagsIn, const vec& priorCompIn) :
    data(dataIn), epsilon(epsIn), mean(meanIn), Linv(LinvIn), lags(lagsIn), priorComp(priorCompIn)  {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = data.n_rows;
    int dim = 2 + 2 * lags;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta1(dim);
    for(int i = 0; i < dim; ++i){
      theta1(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta1(i) += lambda(dim*(i+2) + j) * epsilon(j);
      }
    }
    Matrix<T, Dynamic, 1> theta2(dim);
    for(int i = 0; i < dim; ++i){
      theta2(i) = lambda(dim + i);
      for(int j = 0; j <= i; ++j){
        theta2(i) += lambda(dim*(dim + 2 + i) + j) * epsilon(dim + j);
      }
    }
    T w = 1.0 / (1 - exp(-lambda(2*dim*(dim+1))));
    Matrix<T, Dynamic, 1> theta = w * theta1  +  (1-w) * theta2;  
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Evaluate log(p(theta)), start by evaluative the quadrating in the MVN exponents
    Matrix<T, Dynamic, Dynamic> kernel(dim, 2);
    kernel.fill(0);
    T exponent1 = 0, exponent2 = 0;
    for(int i = 0; i < dim; ++i){
      for(int j = i; j < dim; ++j){
        kernel(i, 0) += (theta(j) - mean(j)) * Linv(j, i);
        kernel(i, 1) += (theta(j) - mean(j+dim)) * Linv(j+dim, i);
      }
      exponent1 += pow(kernel(i, 0), 2);
      exponent2 += pow(kernel(i, 1), 2);
    }
    // Put it together with weights and precalculated det(2*pi*Sig)^(-0.5) as priorComp
    T prior = log(priorComp(0) * priorComp(1) * exp(-0.5 * exponent1)  +
      (1 - priorComp(0)) * priorComp(2) * exp(-0.5 * exponent2));
    
    
    // Evaluate Log Det J
    T logdetJ = theta(0)  +  theta(1)  +  log(w*lambda(dim) + (1-w)*lambda(dim*(dim+2)))  +
      log(w*lambda(2*dim+1) + (1-w)*lambda(dim*(dim+3)+1)) +  (2*lags) * (log(w) + log(1-w));
    for(int i = 2; i < dim; ++i){
      logdetJ += log(fabs(lambda(dim*(i+2) + i)))  +  log(fabs(lambda(dim*(dim + 2 + i) + i)));
    }
    // Evaluate log likelihood
    T logLik = 0;
    Matrix<T, Dynamic, Dynamic> loglikKernel(N, 2);
    for(int t = lags; t < N; ++t){
      loglikKernel(t, 0) = data(t, 0);
      loglikKernel(t, 1) = data(t, 1);
      for(int i = 1; i <= lags; ++i){
        loglikKernel(t, 0) -= data(t-i, 0) * theta(2*i);
        loglikKernel(t, 1) -= data(t-i, 1) * theta(2*i+1);
      }
      logLik += - 0.5 * theta(0) - 0.5 * theta(1) - 
        pow(loglikKernel(t, 0), 2) / (2 * sigSqV) -
        pow(loglikKernel(t, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List arUpdaterMix(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec mean, mat Linv, int lags, vec priorComp){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1>  grad(13 + 20 * lags + 8 * lags * lags);
  // Autodiff
  arUpdateMix p(data, epsilon, mean, Linv, lags, priorComp);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}


