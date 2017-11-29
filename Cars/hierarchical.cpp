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

struct hierAr {
  const mat data;
  const vec epsilon;
  const vec hyperMean;
  const mat hyperLinv;
  const mat Linv;
  const vec obsSum;
  const int lags;
  const bool diag;
  hierAr(const mat& dataIn, const vec& epsIn, const vec& hypMIn, const mat& hypLiIn, const mat& LinvIn, const vec& obsIn, const int& lagsIn, const bool& diagIn) :
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
      for(int j = 0; j <= i; ++j){
        hypKernel(i) += (thetaHat(j) - hyperMean(j)) * hyperLinv(i, j);
      }
      hyperprior += - 0.5 * pow(hypKernel(i), 2);
    }
    
    // Evaluate log(p(theta | theta hat))
    T prior = 0;
    Matrix<T, Dynamic, Dynamic> kernel(dim, N);
    kernel.fill(0);
    for(int k = 0; k < N; ++k){
      for(int i = 0; i < dim; ++i){
        for(int j = 0; j <= i; ++j){
          kernel(i, k) += (theta(j , k) - thetaHat(j)) * Linv(i, j);
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
Rcpp::List hierArDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperMean, mat hyperLinv, mat Linv, vec obsSum, int lags, bool diag = false){
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
  hierAr p(data, epsilon, hyperMean, hyperLinv, Linv, obsSum, lags, diag);
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
  const vec mean;
  const mat Linv;
  const vec logdets;
  const vec weights;
  arUpdateMix(const mat& dataIn, const vec& meanIn, const mat& LinvIn, const vec& logdetsIn, const vec& weightsIn) :
    data(dataIn), mean(meanIn), Linv(LinvIn), logdets(logdetsIn), weights(weightsIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& theta)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = data.n_rows;
    int K = weights.n_elem;
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
    Matrix<T, Dynamic, Dynamic> kernel(6, K);
    kernel.fill(0);
    Matrix<T, Dynamic, 1> exponents(K);
    exponents.fill(0);
    T prior = 0;
    for(int k = 0; k < K; ++k){
      for(int i = 0; i < 6; ++i){
        for(int j = 0; j <= i; ++j){
          kernel(i, k) += (theta(j) - mean(k*6 + j)) * Linv(k*6 + i, j);
        }
        exponents(k) += pow(kernel(i, k), 2);
      }
      prior += weights(k) * logdets(k) * exp(-0.5 * exponents(k));
    }
    // Put it together with weights and precalculated det(2*pi*Sig)^(-0.5) as priorComp
    prior = log(prior);
    
    // Evaluate log likelihood
    T logLik = 0;
    for(int t = 2; t < N; ++t){
      logLik += - 0.5 * (theta(0) + theta(1)) - pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (2 * sigSqV) -
        pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (2 * sigSqD);
    }
    return prior + logLik;
  }
};

mat dThetadLambda(mat lambda, vec eps){
  int K = lambda.n_cols;
  mat out(12, K, fill::ones);
  for(int k = 0; k < 6; ++k){
    out.row(k+6) = eps(k);
  }
  return out;
}

mat dJdLambda(mat lambda, int k){
  int K = lambda.n_cols;
  mat out(12, K, fill::zeros);
  for(int i = 6; i < 12; ++i){
    out(i, k) = 1.0 / std::abs(lambda(i, k));
  }
  return out;
}

mat dThetadPi(Matrix<double, Dynamic, 1> theta, vec pi, mat lambda){
  int K = pi.n_elem;
  mat out(6, K);
  double f;
  for(int i = 0; i < 6; ++i){
    f = 0;
    for(int k = 0; k < K; ++k){
      boost::math::normal_distribution<> dist(lambda(k), lambda(6+k));
      out(i, k) = - cdf(dist, theta(i));
      f += pi(k) * pdf(dist, theta(i));
    }
    out.row(i) /= f;
  }
  return out;
}

mat dPidZ(vec z){
  int K = z.n_elem;
  mat out(K, K);
  vec ez = exp(z);
  double sumSq = std::pow(sum(ez), 2);
 
  for(int i = 0; i < K; ++i){
    for(int j = i; j < K; ++j){
      if(i == j){
        double other = sumSq - ez(i);
        out(i, j) = ez(i) * other / sumSq;
      } else {
        out(i, j) = out(j, i) = ez(i) * ez(j) / sumSq;
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List arUpdaterMix(mat data, mat lambda, vec pi, vec eps, vec z, Rcpp::NumericMatrix thetaIn, int k, vec mean, mat Linv, vec logdets, vec weights){
  Map<MatrixXd> theta(Rcpp::as<Map<MatrixXd> >(thetaIn));
  int K = lambda.n_cols;
  double eval;
  Matrix<double, Dynamic, 1>  thetaGrad(6);
  
  // Autodiff log p wrt theta
  arUpdateMix p(data, mean, Linv, logdets, weights);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, theta, eval, thetaGrad);

  // diff theta wrt lambda
  mat lambdaGrad = dThetadLambda(lambda, eps);
  
  // diff log J wrt lambda
  mat jGrad = dJdLambda(lambda, k);
  
  // combine lambda gradients
  for(int i = 0; i < 6; ++i){
    lambdaGrad.row(i) *= thetaGrad(i);
    lambdaGrad.row(6+i) *= thetaGrad(i);
  }
  lambdaGrad += jGrad;
  
  // diff theta wrt pi
  mat piGrad = dThetadPi(theta, pi, lambda);
  
  // diff pi wrt z
  mat zGrad = dPidZ(z);
  
  // combine z gradients
  vec zGradFull(K, fill::zeros);
  for(int i = 0; i < K; ++i){
    for(int k = 0; k < K; ++k){
      for(int j = 0; j < 6; ++j){
        zGradFull(i) += thetaGrad(j) * piGrad(j, k) * zGrad(k, i);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("Lgrad") = lambdaGrad,
                            Rcpp::Named("Zgrad") = zGradFull,
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
      draws = MCMCdraws.row(MCMCdraws.n_rows - 1000 +m).t();
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