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

using namespace std;
using namespace arma;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// evaluate and automatially differentiate log(p(x_1:T, theta)) + log det J  wrt lambda
struct logPJ {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  logPJ(const mat& dataIn, const vec& epsIn, const vec& hypIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt;
    
    int N = data.n_rows;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(4);
    for(int i = 0; i < 4; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(4*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Constrained to (0, 1)
    T phi = 1.0 / (1 + exp(-theta(2))), gamma = 1.0 / (1 + exp(-theta(3)));
    double pi = 3.14159;
    
    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  +
      (hyperParams(4) - 1) * log(phi)  +  (hyperParams(5)-1) * log(1-phi)  +
      (hyperParams(6) - 1) * log(gamma)  +  (hyperParams(7)-1) * log(1-gamma);
    
    // Evaluate Log Det J
    T logdetJ = log(lambda(4))  +  log(lambda(9))  +  log(lambda(14))  +  log(lambda(19))  +  lambda(0)  +  theta(0)  +  theta(1)  +
      theta(2)  +  theta(3)  - 2 * log(phi + 1)  - 2 * log(gamma + 1);
    
    // Evaluate log likelihood
    T logLik = 0;
    for(int t = 2; t < N; ++t){
      logLik += - 0.5 * theta(0) - pow((data(t, 0) - data(t-1, 0)) - phi * (data(t-1, 0) - data(t-2, 0)), 2) / (2*sigSqV);
      logLik += - 0.5 * theta(1) - pow(data(t, 1) - pi/2 - gamma * (data(t-1, 1) - pi/2), 2) / (2*sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

struct hamiltonFilter {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  const double M;
  hamiltonFilter (const mat& dataIn, const vec& epsIn, const vec& hypIn, const double& MIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn), M(MIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    
    int N = data.n_rows;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(7);
    for(int i = 0; i < 7; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(7*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Constrained to (0, 1)
    //T phi = 1.0 / (1 + exp(-theta(2)));
    // Unconstrained
    T zeta = theta(2), a0 = theta(3), b0 = theta(4), a1 = theta(5), b1 = theta(6);
    double pi = 3.14159;
    
    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  -
      // (hyperParams(4) - 1) * log(phi)  +  (hyperParams(5)-1) * log(1-phi)  -
      pow(zeta - hyperParams(4), 2) / (2 * hyperParams(5)) - 
      pow(a0 - hyperParams(6), 2) / (2 * hyperParams(7)) - 
      pow(b0 - hyperParams(8), 2) / (2 * hyperParams(9)) - 
      pow(a1 - hyperParams(10), 2) / (2 * hyperParams(11)) - 
      pow(b1 - hyperParams(12), 2) / (2 * hyperParams(13));
    
    // Evaluate Log Det J
    T logdetJ = theta(0)  +  theta(1);// -  theta(2)  +  2 * log(phi);
    for(int i = 0; i < 7; ++i){
      logdetJ += log(fabs(lambda(8*i + 7)));
    }
    
    // Initialise the filter
    Matrix<T, Dynamic, 1> xi(N), rho01(N), rho11(N), pv(N), pd0(N), pd1(N), likelihood(N);
    xi(0) = 0.95;
    T loglik = 0;
    for(int t = 1; t < N; ++t){
      rho01(t) = 1.0 / (1 + exp(-(a0 + b0 * pow(data(t, 0) - M, 2))));
      rho11(t) = 1.0 / (1 + exp(-(a1 + b1 * pow(data(t, 0) - M, 2))));
      pv(t) = - 0.5 * log(2 * pi * sigSqV)  -  pow(data(t, 1) - data(t-1, 1), 2) / (2*sigSqV);// - phi * (data(t-1, 1) - data(t-2, 1)), 2)/(2*sigSqV);
      pd0(t) = -0.5 * log(2 * pi * sigSqD)  -  pow(data(t, 2) - data(t-1, 2), 2) / (2*sigSqD);
      pd1(t) = -0.5 * log(2 * pi * sigSqD)  -  pow(data(t, 2) - pi/2 - zeta * (data(t, 0) - M), 2) / (2*sigSqD);
      likelihood(t) =  xi(t-1) * (1 - rho01(t)) * exp(pv(t) + pd0(t))  +   //i = 0, j = 0
        xi(t-1) * rho01(t) * exp(pv(t) + pd1(t))  +                        //i = 0, j = 1
        (1 - xi(t-1)) * (1 - rho11(t)) * exp(pv(t) + pd0(t))  +            //i = 1, j = 0
        (1 - xi(t-1)) * rho11(t) * exp(pv(t) + pd1(t));                    //i = 1, j = 1
      xi(t) = ((1 - rho01(t)) * xi(t-1) * exp(pv(t) + pd0(t))  +  (1 - rho11(t)) * (1 - xi(t-1)) * exp(pv(t) + pd0(t))) / likelihood(t);
      //Rcpp::Rcout << pd0(t) << " " << pd1(t) << std::endl;
      loglik += log(likelihood(t));

    }
    
    // Evaluate log(q(eps))
    T logQ  = 0;
    for(int i = 0; i < 7; ++i){
      logQ += -pow(epsilon(i), 2) / 2;
    }
    return prior + loglik + logdetJ - logQ;
  }
};

struct arimaVD {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  arimaVD(const mat& dataIn, const vec& epsIn, const vec& hypIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt;
    int N = data.n_rows;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(6);
    for(int i = 0; i < 6; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(6*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Constrained to (0, 1)
    T arV = 1.0 / (1 + exp(-theta(2))), maV = 1.0 / (1 + exp(-theta(3))),
      arD = 1.0 / (1 + exp(-theta(4))), maD = 1.0 / (1 + exp(-theta(5)));
    double pi = 3.14159;
    // Evaluate log(p(theta))
    T prior =  -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  +
      (hyperParams(4) - 1) * log(arD)  +  (hyperParams(5)-1) * log(1-arD)  +
      (hyperParams(6) - 1) * log(arV)  +  (hyperParams(7)-1) * log(1-arV)  +
      (hyperParams(8) - 1) * log(maD)  +  (hyperParams(9)-1) * log(1-maD)  +
      (hyperParams(10) - 1) * log(maV)  +  (hyperParams(11)-1) * log(1-maV);
    
    // Evaluate Log Det J
     T logdetJ = theta(0)  +  theta(1)  -  theta(2)  -  theta(3)  -  theta(4)  -  theta(5)  +  
       2 * (log(arD)  +  log(arV)  + log(maD)  + log(maV));
    for(int i = 0; i < 6; ++i){
      logdetJ += theta(i) + log(fabs(lambda(7*i+6)));
    }
    // Filter errors
    Matrix<T, Dynamic, 1> epsilonV(N), epsilonD(N);
    epsilonV(0) = 0;
    epsilonD(0) = 0;
    for(int t = 1; t < N; ++t){
      epsilonV(t) = data(t, 0)  -  arV * data(t-1, 0)  -  maV * epsilonV(t-1);
      epsilonD(t) = data(t, 1)  -  pi/2   -  arD * (data(t-1, 1) - pi/2)  -  maD * epsilonD(t-1);
    }
    // Evaluate log likelihood
    T logLik = 0;
    for(int t = 2; t < N; ++t){
      logLik += - 0.5 * theta(0) - pow(data(t, 0) - arV * data(t-1, 0) - maV * epsilonV(t-1), 2) / (2*sigSqV);
      logLik += - 0.5 * theta(1) - pow(data(t, 1) - pi/2 - arD * (data(t-1, 1) - pi/2) - maD * epsilonD(t-1), 2) / (2*sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};

struct ar1VD {
  const mat data;
  const vec epsilon;
  const vec hyperParams;
  ar1VD(const mat& dataIn, const vec& epsIn, const vec& hypIn) :
    data(dataIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt;
    int N = data.n_rows;
    // Create theta as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> theta(4);
    for(int i = 0; i < 4; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(4*(i+1) + j) * epsilon(j);
      }
    }
    // Constrained Positive
    T sigSqV = exp(theta(0)), sigSqD = exp(theta(1));
    // Constrained to (-1, 1)
    T arV = 2.0 / (1 + exp(-theta(2))) - 1, arD = 2.0 / (1 + exp(-theta(3))) - 1;
    double pi = 3.14159;
    // Evaluate log(p(theta))
    T prior =  -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  +
      (hyperParams(4) - 1) * log(1+arD)  +  (hyperParams(5)-1) * log(1-arD)  +
      (hyperParams(6) - 1) * log(1+arV)  +  (hyperParams(7)-1) * log(1-arV);
    
    // Evaluate Log Det J
    T logdetJ = theta(0)  +  theta(1)  -  theta(2)  -  theta(3)  +  2 * (log(arD)  +  log(arV));
    for(int i = 0; i < 4; ++i){
      logdetJ += log(fabs(lambda(5*i+4)));
    }
    // Evaluate log likelihood
    T logLik = 0;
    for(int t = 1; t < N; ++t){
      logLik += - 0.5 * theta(0) - pow(data(t, 0) - arV * data(t-1, 0), 2) / (2*sigSqV);
      logLik += - 0.5 * theta(1) - pow(data(t, 1) - pi/2 - arD * (data(t-1, 1) - pi/2), 2) / (2*sigSqD);
    }
    return prior + logLik + logdetJ;
  }
};
  
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

// [[Rcpp::export]]
Rcpp::List PJDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperParams){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1> grad(20);
  // Autodiff
  logPJ p(data, epsilon, hyperParams);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// [[Rcpp::export]]
Rcpp::List HFDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperParams, double M = 5){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1> grad(56);
  // Autodiff
  hamiltonFilter hf(data, epsilon, hyperParams, M);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(hf, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// [[Rcpp::export]]
Rcpp::List arimaDeriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperParams){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1> grad(42);
  // Autodiff
  arimaVD p(data, epsilon, hyperParams);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}


// [[Rcpp::export]]
Rcpp::List ar1Deriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperParams){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1> grad(20);
  // Autodiff
  ar1VD p(data, epsilon, hyperParams);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}



double MCMCfilter (mat data, rowvec theta, vec hyperParams, double M){
  int N = data.n_rows;
  double xi, rho01, rho11, pv, pd0, pd1, likelihood;
  xi = 0.95;
  double loglik = 0, sigSqV = theta(0), sigSqD = theta(1), zeta = theta(2), a0 = theta(3), b0 = theta(4), a1 = theta(5), b1 = theta(6), pi = 3.14159;
  
  
  for(int t = 1; t < N; ++t){
    rho01 = 1.0 / (1 + exp(-(a0 + b0 * pow(data(t, 0) - M, 2))));
    rho11 = 1.0 / (1 + exp(-(a1 + b1 * pow(data(t, 0) - M, 2))));
    pv = - 0.5 * log(2 * pi * sigSqV)  -  pow(data(t, 1) - data(t-1, 1), 2) / (2*sigSqV);// - phi * (data(t-1, 1) - data(t-2, 1)), 2)/(2*sigSqV);
    pd0 = -0.5 * log(2 * pi * sigSqD)  -  pow(data(t, 2) - data(t-1, 2), 2) / (2*sigSqD);
    pd1 = -0.5 * log(2 * pi * sigSqD)  -  pow(data(t, 2) - pi/2 - zeta * (data(t, 0) - M), 2) / (2*sigSqD);
    likelihood =  xi * (1 - rho01) * exp(pv + pd0)  +   //i = 0, j = 0
      xi * rho01 * exp(pv + pd1)  +                        //i = 0, j = 1
      (1 - xi) * (1 - rho11) * exp(pv + pd0)  +            //i = 1, j = 0
      (1 - xi) * rho11 * exp(pv + pd1);                    //i = 1, j = 1
    xi = ((1 - rho01) * xi * exp(pv + pd0)  +  (1 - rho11) * (1 - xi) * exp(pv + pd0)) / likelihood;
    loglik += log(likelihood);
  }
  
  double prior = -(hyperParams(0) + 1) * log(sigSqV)  -  hyperParams(1) / sigSqV  - 
    (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  -
    // (hyperParams(4) - 1) * log(phi)  +  (hyperParams(5)-1) * log(1-phi)  -
    pow(zeta - hyperParams(4), 2) / (2 * hyperParams(5)) - 
    pow(a0 - hyperParams(6), 2) / (2 * hyperParams(7)) - 
    pow(b0 - hyperParams(8), 2) / (2 * hyperParams(9)) - 
    pow(a1 - hyperParams(10), 2) / (2 * hyperParams(11)) - 
    pow(b1 - hyperParams(12), 2) / (2 * hyperParams(13));
  
  return loglik + prior;
}

// [[Rcpp::export]]
mat HFMCMC (mat data, int reps, rowvec stepsize, vec hyperParams, double Mid){
  int N = data.n_rows;
  mat theta(reps, 7);
  rowvec initial = {0.01, 0.00001, 0.05, -5, 1, 3, 12};
  theta.row(0) = initial;
  int accept = 0;
  double loglik = MCMCfilter(data, initial, hyperParams, Mid);
  for(int iter = 1; iter < reps; ++iter){
    rowvec candidate = theta.row(iter-1) + stepsize % randn<rowvec>(7);
    double canLoglik = MCMCfilter(data, candidate, hyperParams, Mid);
    double ratio = exp(canLoglik - loglik);
    double u = randu<vec>(1)[0];
    if(u < ratio){
      accept += 1;
      theta.row(iter) = candidate;
      loglik = canLoglik;
    } else {
      theta.row(iter) = theta.row(iter-1);
    }
    if(iter % 10000 == 0){
      Rcpp::Rcout << "Iteration: " << iter << std::endl;
    }
  }
  Rcpp::Rcout << "Iteration: " << reps << std::endl;
  Rcpp::Rcout << accept * 1.0 / reps << std::endl;
  return theta; 
}

