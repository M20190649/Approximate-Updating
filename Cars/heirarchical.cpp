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

struct heirAr1 {
  const mat data;
  const vec epsilon;
  const vec hyperMean;
  const mat hyperLinv;
  const vec obsSum;
  const mat Linv;
  heirAr1(const mat& dataIn, const vec& epsIn, const vec& hypMIn, const mat& hypLiIn, const mat& LiIn, const vec& obsIn) :
    data(dataIn), epsilon(epsIn), hyperMean(hypMIn), hyperLinv(hypLiIn), Linv(LiIn), obsSum(obsIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::fabs;
    int N = obsSum.n_elem;
    // Create thetaHat as Mu + U %*% Eps
    Matrix<T, Dynamic, 1> thetaHat(4);
    T logdetJ = 0;
    int lCounter = 4*(N+1);
    for(int i = 0; i < 4; ++i){
      thetaHat(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        thetaHat(i) += lambda(lCounter) * epsilon(j);
        if(i == j){
          logdetJ += log(fabs(lambda(lCounter)));
        }
        lCounter += 1;
      }
    }
    //Matrix<T, Dynamic, Dynamic> L(4, 4);
    //int vCounter = 0;
    //for(int i = 0; i < 4; ++i){
    //  for(int j = 0; j <= i; ++j)
    //    L(i, j) = lambda(4 + vCounter) + lambda(lCounter) * epsilon(4 + vCounter);
    //    lCounter += 1;
    //    vCounter += 1;
    //}
    
    Matrix<T, Dynamic, Dynamic> theta(4, N);
    for(int k = 0; k < N; ++k){
      for(int i = 0; i < 4; ++i){
        theta(i, k) = lambda(4 + 4*k + i);
        for(int j = 0; j <= i; ++j){
          theta(i, k) += lambda(lCounter) * epsilon(4 + 4*k + j);
          if(i == j){
            logdetJ += log(fabs(lambda(lCounter)));
          }
          lCounter += 1;
        }
      }
    }
    Matrix<T, Dynamic, 1> sigSqV(N), sigSqD(N), arV(N), arD(N);
    for(int k = 0; k < N; ++k){
      sigSqV(k) = exp(theta(0, k));
      sigSqD(k) = exp(theta(1, k));
      arV(k) = 2.0 / (1 + exp(-theta(2, k))) - 1;
      arD(k) = 2.0 / (1 + exp(-theta(3, k))) - 1;
    }
    double pi = 3.14159;
    // Evaluate Log Det J
    for(int k = 0; k < N; ++k){
      logdetJ += theta(0, k) + theta(1, k) - theta(2, k) - theta(3, k) + 2 * log(1+arV(k)) + 2 * log(1+arD(k));
    }
    // Evaluate log(p(theta hat))
    T hyperprior =  -log(fabs(hyperLinv(0, 0))) - log(fabs(hyperLinv(1, 1)))  -log(fabs(hyperLinv(2, 2))) - log(fabs(hyperLinv(3, 3))) -
           0.5 * (pow((thetaHat(0) - hyperMean(0))*hyperLinv(0, 0) + (thetaHat(1) - hyperMean(1))*hyperLinv(1, 0) + 
                      (thetaHat(2) - hyperMean(2))*hyperLinv(2, 0) + (thetaHat(3) - hyperMean(3))*hyperLinv(3, 0), 2) + 
                  pow((thetaHat(1) - hyperMean(1))*hyperLinv(1, 1) + (thetaHat(2) - hyperMean(2))*hyperLinv(2, 1) + 
                      (thetaHat(3) - hyperMean(3))*hyperLinv(3, 1), 2) + 
                  pow((thetaHat(2) - hyperMean(2))*hyperLinv(2, 2) + (thetaHat(3) - hyperMean(3))*hyperLinv(3, 2), 2) + 
                  pow((thetaHat(3) - hyperMean(3))*hyperLinv(3, 3), 2));

    // L * t(L) ~ IW(hyperWish)

    // Invert the L matrix
    //Matrix<T, Dynamic, Dynamic> Linv(4, 4);
    //for(int i = 0; i < 4; ++i){
    //  for(int j = i; j >= 0; --j){
    //    if(i == j){
    //      Linv(i, j) = 1.0 / L(i, j);
    //    } else if(i == j + 1) {
    //      Linv(i, j) = - L(i, j) * Linv(j, j) / L(i, i);
    //    } else if(i == j + 2) {
    //      Linv(i, j) = - (L(i, j) * Linv(j, j) + L(i, j+1) * Linv(j+1, j)) / L(i, i);
    //    } else if(i == j + 3) {
    //      Linv(i, j) = - (L(i, j) * Linv(j, j) + L(i, j+1) * Linv(j+1, j) + L(i, j+2) * Linv(j+2, j)) / L(i, i);
    //    }
    //  }
    //}
    // Evaluate log(p(theta | theta hat))
    T prior = 0;
    for(int k = 0; k < N; ++k){
      prior += -log(fabs(Linv(0, 0))) - log(fabs(Linv(1, 1))) - log(fabs(Linv(2, 2))) - log(fabs(Linv(3, 3)))  -
        0.5 * (pow((theta(0, k) - thetaHat(0))*Linv(0, 0) + (theta(1, k) - thetaHat(1))*Linv(1, 0) + 
        (theta(2, k) - thetaHat(2))*Linv(2, 0) + (theta(3, k) - thetaHat(3))*Linv(3, 0), 2) + 
        pow((theta(1, k) - thetaHat(1))*Linv(1, 1) + (theta(2, k) - thetaHat(2))*Linv(2, 1) + (theta(3, k) - thetaHat(3))*Linv(3, 1), 2) + 
        pow((theta(2, k) - thetaHat(2))*Linv(2, 2) + (theta(3, k) - thetaHat(3))*Linv(3, 2), 2) + 
        pow((theta(3, k) - thetaHat(3))*Linv(3, 3), 2));
    }
    
    // Evaluate log likelihood
    T logLik = 0;
    for(int t = 1; t < obsSum(0); ++t){
      logLik += - 0.5 * theta(0, 0) - pow(data(t, 0) - arV(0) * data(t-1, 0), 2) / (2*sigSqV(0));
      logLik += - 0.5 * theta(1, 0) - pow(data(t, 1) - pi/2 - arD(0) * (data(t-1, 1) - pi/2), 2) / (2*sigSqD(0));
    }
    for(int k = 1; k < N; ++k){
      for(int t = obsSum(k-1) + 1; t < obsSum(k); ++t){
        logLik += - 0.5 * theta(0, k) - pow(data(t, 0) - arV(k) * data(t-1, 0), 2) / (2*sigSqV(k));
        logLik += - 0.5 * theta(1, k) - pow(data(t, 1) - pi/2 - arD(k) * (data(t-1, 1) - pi/2), 2) / (2*sigSqD(k));
      }
    }
    return hyperprior + prior + logLik + logdetJ;
  }
};

// [[Rcpp::export]]
Rcpp::List heirAr1Deriv(mat data, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperMean, mat hyperLinv, mat Linv, vec obsPer){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  int N = obsPer.n_elem;
  Matrix<double, Dynamic, 1> grad(14 * (N + 1));
  // Autodiff
  heirAr1 p(data, epsilon, hyperMean, hyperLinv, Linv, obsPer);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}


