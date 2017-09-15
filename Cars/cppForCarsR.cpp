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

struct logPZ {
  const mat x;
  const int obs;
  const int N;
  const vec epsilon;
  const vec hyperParams;
  logPZ(const mat& xIn, const int& obsIn, const int& nIn, const vec& epsIn, const vec& hypIn) :
    x(xIn), obs(obsIn), N(nIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambdaV)
    const{
    using std::log; using std::exp; using std::pow; using std::cos; using std::sin;
    
    // Create theta
    Matrix<T, Dynamic, 1> theta(14);
    for(int i = 0; i < 14; ++i){
      theta(i) = lambdaV(i) + sqrt(pow(lambdaV(14+i), 2)) * epsilon(i);
    }
    // Constrained Positive
    T sigSqA = exp(theta(0)), sigSqD = exp(theta(1)), sigSqX = exp(theta(2)), sigSqY = exp(theta(3));
    // Constrained to (0, 1)
    T lambda = 1.0 / (1 + exp(-theta(4))), phi = 1.0 / (1 + exp(-theta(5))), sT1 = 1.0 / (1 + exp(-theta(12)));
    // Unconstrained
    T gamma = theta(6), alpha0 = theta(7), beta0 = theta(8), alpha1 = theta(9), beta1 = theta(10), aT1 = theta(11), dT1 = theta(13);
    double pi = 3.14159;
    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqA)  -  hyperParams(1) / sigSqA  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  -
      (hyperParams(4) + 1) * log(sigSqX)  -  hyperParams(5) / sigSqX -
      (hyperParams(6) + 1) * log(sigSqY)  -  hyperParams(7) / sigSqY + 
      (hyperParams(8) - 1) * log(lambda)  +  (hyperParams(9)-1) * log(1-lambda)  +
      (hyperParams(10) - 1) * log(phi)  +  (hyperParams(11)-1) * log(1-phi)  -  
      pow(gamma - hyperParams(12), 2) / (2*hyperParams(13))  -
      pow(alpha0 - hyperParams(14), 2) / (2*hyperParams(15))  -
      pow(beta0 - hyperParams(16), 2) / (2*hyperParams(17))  -
      pow(alpha1 - hyperParams(18), 2) / (2*hyperParams(19))  -
      pow(beta1 - hyperParams(20), 2) / (2*hyperParams(21)); 
    // Set up for particle filtering
    T xAllDens, xtDens, omegaSum, ZT1Dens;
    Matrix<T, Dynamic, Dynamic> States(N, 6), Resample(N, 6), Ahead(N, 6);
    Matrix<T, Dynamic, 1> pweights(N), piSum(N), omega(N), lweights(N), lambdaSum(N);
    Matrix<int, Dynamic, 1> kDraws(N);
    pweights.fill(1.0/N); // initial weights
    
    // Sample z0 from stationary distribution
    for(int k = 0; k < N; ++k){
      States(k, 0) = sqrt(sigSqA / (1 - pow(phi, 2))) * randn<vec>(1)[0];
      States(k, 1) = 0;
      States(k, 2) = pi / 2  +  sqrt(sigSqD) * randn<vec>(1)[0];
      States(k, 3) = 5.0 + 0.01 * randn<vec>(1)[0];
      States(k, 4) = 5.0 + 0.01 * randn<vec>(1)[0];
      States(k, 5) = 0.1 * randn<vec>(1)[0];
    }
    xAllDens = 0;
    
    // Particle Filter loop
    for(int t = 0; t < obs; ++t){
      // Create CDF to use for resampling
      piSum(0) = pweights(0);
      for(int k = 1; k < N; ++k){
        piSum(k) = pweights(k) + piSum(k-1);
      }
      // Resample xt via inverse CDF
      for(int k = 0; k < N; ++k){
        double u = randu<vec>(1)[0];
        for(int i = 0; i < N; ++i){
          if(u < piSum(i)){
            for(int j = 0; j < 6; ++j){
              Resample(k, j) = States(i, j);
            }
            //Resample.block(1, 6, k, 6*t) = Old.block(1, 6, i, 6*t);
            break;
          }
        }
      }
      
      // Calculate Step Ahead Means and k density
      for(int k = 0; k < N; ++k){
        Ahead(k, 0) = phi * Resample(k, 0);
        Ahead(k, 1) = Resample(k, 1);
        Ahead(k, 2) = pi/2  +  (1 - Ahead(k, 1)) * lambda * (Resample(k, 2) - pi/2)  +
          Ahead(k, 1) * gamma * (Resample(k, 4) - 5);
        Ahead(k, 3) = Resample(k, 3)  +  Ahead(k, 0);
        Ahead(k, 4) = Resample(k, 4)  +   Ahead(k, 3) * cos(Ahead(k, 2));
        Ahead(k, 5) = Resample(k, 5)  +   Ahead(k, 3) * sin(Ahead(k, 2));
        lweights(k) = 1.0 / sqrt(2 * pi * sigSqX) * exp(-pow(x(t, 0) - Ahead(k, 4), 2) / (2*sigSqX)) *
          1.0 / sqrt(2 * pi * sigSqY) * exp(-pow(x(t, 1) - Ahead(k, 5), 2) / (2*sigSqY));
        if(k == 0){
          lambdaSum(k) = lweights(k);
        } else {
          lambdaSum(k) = lweights(k) + lambdaSum(k-1);
        }
      }
      // Normalise
      for(int k = 0; k < N; ++k){
        lambdaSum(k) = lambdaSum(k) / lambdaSum(N-1);
      }
      
      
      // Draw index k using lambda as weights
      for(int k = 0; k < N; ++k){
        double u = randu<vec>(1)[0];
        for(int i = 0; i < N; ++i){
          if(u < lambdaSum(i)){
            kDraws(k) = i;
          }
        }
      }
      
      // Calculate weights
      omegaSum = 0;
      for(int k = 0; k < N; ++k){
        // Step Ahead transitions
        States(k, 0) = phi * Resample(kDraws(k), 0)  +  sqrt(sigSqA) * randn<vec>(1)[0];
        double u = randu<vec>(1)[0];
        if((Resample(kDraws(k), 1) == 0 & u < 1.0 / (1 + exp(- (alpha0 + beta0 * abs(Resample(kDraws(k), 4) - 5))))) |
           (Resample(kDraws(k), 1) == 1 & u < 1.0 / (1 + exp(- (alpha1 + beta1 * abs(Resample(kDraws(k), 4) - 5)))))){
          States(k, 1) = 1;
        } else {
          States(k, 1) = 0;
        }
        States(k, 2) =  pi/2  +  (1 - States(k, 1)) * lambda * (Resample(kDraws(k), 2) - pi/2)  +
          States(k, 1) * gamma * (Resample(kDraws(k), 4) - 5)  +  sqrt(sigSqD) * randn<vec>(1)[0];
        States(k, 3) = Resample(kDraws(k), 3)  +  States(k, 0);
        States(k, 4) = Resample(kDraws(k), 4)  +  States(k, 3) * cos(States(k, 2));
        States(k, 5) = Resample(kDraws(k), 5)  +  States(k, 3) * sin(States(k, 2));
        // Measurement density
        omega(k) =  1.0 / sqrt(2 * pi * sigSqX) * exp(-pow(x(t, 0) - States(k, 4), 2) / (2*sigSqX)) *
          1.0 / sqrt(2 * pi * sigSqY) * exp(-pow(x(t, 1) - States(k, 5), 2) / (2*sigSqY));
        // sum of weights for normalisation
        omegaSum += omega(k);
      }
      
      // Normalise weights
      for(int k = 0; k < N; ++k){
        pweights(k) = omega(k) / omegaSum;
      }
      // log(p(x_1:T | theta)) = sum_t log(p(x_t | theta))
      xtDens = log(omegaSum / N);
      xAllDens += xtDens;
    }
    // end of particle filter loop
    // step ahead latent states density
    ZT1Dens = 0; 
    Matrix<T, Dynamic, 1> atDens(N), dtDens(N), stDens(N);
    // p(V_T+1 | theta, x_1:T) approx sum_k pi_k p(V_T+1 | V_T,k, theta)
    for(int k = 0; k < N; ++k){
      atDens(k) = 1.0 / sqrt(2*pi*sigSqA) * exp(-pow(aT1 - phi * States(k, 0), 2) / (2*sigSqA));
      dtDens(k) = 1.0 / sqrt(2*pi*sigSqD) * exp(-pow(dT1- pi/2 - (1-sT1)*lambda*(States(k, 2) -pi/2) + sT1*gamma*(States(k, 4) - 5), 2) / (2*sigSqD));
      stDens(k) = (1-States(k, 1)) * pow(1 / (1 + exp(-alpha0 - beta0 * abs(States(k, 4)-5))), sT1) * pow(1 - 1/(1 + exp(-alpha0 - beta0 * abs(States(k, 4)-5))), 1-sT1) + 
        States(k, 1) * pow(1 / (1 + exp(-alpha1 - beta1 * abs(States(k, 4)-5))), sT1) * pow(1 - 1/(1 + exp(-alpha1 - beta1 * abs(States(k, 4)-5))), 1-sT1);
      ZT1Dens += pweights(k) * atDens(k) * dtDens(k) * stDens(k);
    }
    // interested in log density
    ZT1Dens = log(ZT1Dens);
    return xAllDens + ZT1Dens + prior;
  }
};

struct logP{
  const mat x;
  const int obs;
  const int N;
  const vec epsilon;
  const vec hyperParams;
  logP(const mat& xIn, const int& obsIn, const int& nIn, const vec& epsIn, const vec& hypIn) :
    x(xIn), obs(obsIn), N(nIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambdaV)
    const{
    using std::log; using std::exp; using std::pow; using std::cos; using std::sin;
    
    // Create theta
    Matrix<T, Dynamic, 1> theta(11);
    for(int i = 0; i < 11; ++i){
      theta(i) = lambdaV(i) + sqrt(pow(lambdaV(11+i), 2)) * epsilon(i);
    }
    // Constrained Positive
    T sigSqA = exp(theta(0)), sigSqD = exp(theta(1)), sigSqX = exp(theta(2)), sigSqY = exp(theta(3));
    // Constrained to (0, 1)
    T lambda = 1.0 / (1 + exp(-theta(4))), phi = 1.0 / (1 + exp(-theta(5)));
    // Unconstrained
    T gamma = theta(6), alpha0 = theta(7), beta0 = theta(8), alpha1 = theta(9), beta1 = theta(10);
    double pi = 3.14159;
    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqA)  -  hyperParams(1) / sigSqA  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  -
      (hyperParams(4) + 1) * log(sigSqX)  -  hyperParams(5) / sigSqX -
      (hyperParams(6) + 1) * log(sigSqY)  -  hyperParams(7) / sigSqY + 
      (hyperParams(8) - 1) * log(lambda)  +  (hyperParams(9)-1) * log(1-lambda)  +
      (hyperParams(10) - 1) * log(phi)  +  (hyperParams(11)-1) * log(1-phi)  -  
      pow(gamma - hyperParams(12), 2) / (2*hyperParams(13))  -
      pow(alpha0 - hyperParams(14), 2) / (2*hyperParams(15))  -
      pow(beta0 - hyperParams(16), 2) / (2*hyperParams(17))  -
      pow(alpha1 - hyperParams(18), 2) / (2*hyperParams(19))  -
      pow(beta1 - hyperParams(20), 2) / (2*hyperParams(21)); 
    // Set up for particle filtering
    T xAllDens, xtDens, omegaSum;
    Matrix<T, Dynamic, Dynamic> States(N, 6), Resample(N, 6);//, Ahead(N, 6);
    Matrix<T, Dynamic, 1> pweights(N), piSum(N), omega(N); //lweights(N), lambdaSum(N);
    //Matrix<int, Dynamic, 1> kDraws(N);
    pweights.fill(1.0/N); // initial weights
    
    // Sample z0 from stationary distribution
    for(int k = 0; k < N; ++k){
      States(k, 0) = sqrt(sigSqA / (1 - pow(phi, 2))) * randn<vec>(1)[0];
      States(k, 1) = 0;
      States(k, 2) = pi / 2  +  sqrt(sigSqD) * randn<vec>(1)[0];
      States(k, 3) = 5.0 + 0.01 * randn<vec>(1)[0];
      States(k, 4) = 5.0 + 0.01 * randn<vec>(1)[0];
      States(k, 5) = 0.1 * randn<vec>(1)[0];
    }
    xAllDens = 0;
    
    // Particle Filter loop
    for(int t = 0; t < obs; ++t){
      // Create CDF to use for resampling
      piSum(0) = pweights(0);
      for(int k = 1; k < N; ++k){
        piSum(k) = pweights(k) + piSum(k-1);
      }
      // Resample xt via inverse CDF
      for(int k = 0; k < N; ++k){
        double u = randu<vec>(1)[0];
        for(int i = 0; i < N; ++i){
          if(u < piSum(i)){
            for(int j = 0; j < 6; ++j){
              Resample(k, j) = States(i, j);
            }
            //Resample.block(1, 6, k, 6*t) = Old.block(1, 6, i, 6*t);
            break;
          }
        }
      }
      
      // Calculate Step Ahead Means and k density
      //for(int k = 0; k < N; ++k){
      //  Ahead(k, 0) = phi * Resample(k, 0);
      //  Ahead(k, 1) = Resample(k, 1);
      //  Ahead(k, 2) = pi/2  +  (1 - Ahead(k, 1)) * lambda * (Resample(k, 2) - pi/2)  +
      //    Ahead(k, 1) * gamma * (Resample(k, 4) - 5);
      //  Ahead(k, 3) = Resample(k, 3)  +  Ahead(k, 0);
      //  Ahead(k, 4) = Resample(k, 4)  +   Ahead(k, 3) * cos(Ahead(k, 2));
      //  Ahead(k, 5) = Resample(k, 5)  +   Ahead(k, 3) * sin(Ahead(k, 2));
      //  lweights(k) = 1.0 / sqrt(2 * pi * sigSqX) * exp(-pow(x(t, 0) - Ahead(k, 4), 2) / (2*sigSqX)) *
      //    1.0 / sqrt(2 * pi * sigSqY) * exp(-pow(x(t, 1) - Ahead(k, 5), 2) / (2*sigSqY));
      //  if(k == 0){
      //    lambdaSum(k) = lweights(k);
      //  } else {
      //    lambdaSum(k) = lweights(k) + lambdaSum(k-1);
      //  }
      //}
      // Normalise
      //for(int k = 0; k < N; ++k){
      //  lambdaSum(k) = lambdaSum(k) / lambdaSum(N-1);
      //}
      
      
      // Draw index k using lambda as weights
      //for(int k = 0; k < N; ++k){
      //  double u = randu<vec>(1)[0];
      //  for(int i = 0; i < N; ++i){
      //    if(u < lambdaSum(i)){
      //      kDraws(k) = i;
      //    }
      //  }
      //}
      
      // Calculate weights
      omegaSum = 0;
      for(int k = 0; k < N; ++k){
        // Step Ahead transitions
        States(k, 0) = phi * Resample(k, 0)  +  sqrt(sigSqA) * randn<vec>(1)[0];
        States(k, 1) = 1.0 / (1 + exp(- (1-Resample(k, 1))*(alpha0 + beta0 * abs(Resample(k, 4) - 5)) - 
          Resample(k, 1)*(alpha1 + beta1 * abs(Resample(k, 4) - 5)) + 0.1*randn<vec>(1)[0]));
        States(k, 2) =  pi/2  +  (1 - States(k, 1)) * lambda * (Resample(k, 2) - pi/2)  +
          States(k, 1) * gamma * (Resample(k, 4) - 5)  +  sqrt(sigSqD) * randn<vec>(1)[0];
        States(k, 3) = Resample(k, 3)  +  States(k, 0);
        States(k, 4) = Resample(k, 4)  +  States(k, 3) * cos(States(k, 2));
        States(k, 5) = Resample(k, 5)  +  States(k, 3) * sin(States(k, 2));
        // Measurement density
        omega(k) =  1.0 / sqrt(2 * pi * sigSqX) * exp(-pow(x(t, 0) - States(k, 4), 2) / (2*sigSqX)) *
          1.0 / sqrt(2 * pi * sigSqY) * exp(-pow(x(t, 1) - States(k, 5), 2) / (2*sigSqY));
        // sum of weights for normalisation
        omegaSum += omega(k);
      }
      
      // Normalise weights
      for(int k = 0; k < N; ++k){
        pweights(k) = omega(k) / omegaSum;
      }
      // log(p(x_1:T | theta)) = sum_t log(p(x_t | theta))
      xtDens = log(omegaSum / N);
      xAllDens += xtDens;
    }
    // end of particle filter loop
    return xAllDens + prior;
  }
};

struct logDetJ {
  const vec epsilon;
  logDetJ(const vec epsIn) : // constructor
    epsilon(epsIn){}
  template <typename T> // 
  T operator ()(const Matrix<T, Dynamic, 1>& lambda) // derivative is with respect to lambda
    const{
    using std::log; using std::exp; 
    T logdet = 0;
    for(int i = 0; i < 11; ++i){
      if(i >= 6){ // linear
        logdet += log(sqrt(pow(lambda(11+i), 2)));
      } else if(i < 4){ // exponential
        logdet += log(sqrt(pow(lambda(11+i), 2))) + lambda(i) + sqrt(pow(lambda(11+i), 2))*epsilon(i);
      } else { // sigmoid
        logdet += log(sqrt(pow(lambda(11+i), 2))) + lambda(i) + sqrt(pow(lambda(11+i), 2))*epsilon(i) - 2*log(1 + exp(lambda(i) + sqrt(pow(lambda(11+i), 2))*epsilon(i)));
      }
    }
    return logdet;
  }
};

// [[Rcpp::export]]
Rcpp::List pDeriv(mat x, Rcpp::NumericMatrix lambdaIn, vec epsilon, vec hyperParams, int N){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  int T = x.n_rows;
  Matrix<double, Dynamic, 1> grad(22);
  // Autodiff
  logP p(x, T, N, epsilon, hyperParams);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(p, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
}

// [[Rcpp::export]]
Rcpp::List jDeriv(Rcpp::NumericMatrix lambdaIn, vec epsilon){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double eval;
  Matrix<double, Dynamic, 1> grad(22);
  // Autodiff
  logDetJ j(epsilon);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(j, lambda, eval, grad);
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = eval);
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

// generate theta for importance sampling.
// [[Rcpp::export]]
vec constructTheta(Rcpp::NumericMatrix lambdaIn, vec epsilon){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  // Create theta
  vec theta(14);
  for(int i = 0; i < 14; ++i){
    theta(i) = lambda(i) + lambda(14+i) * epsilon(i);
  }
  // Constrained Positive
  theta(0) = exp(theta(0));
  theta(1) = exp(theta(1));
  theta(2) = exp(theta(2));
  theta(3) = exp(theta(3));
  // Constrained to (0, 1)
  theta(4) = 1.0 / (1 + exp(-theta(2)));
  theta(5) = 1.0 / (1 + exp(-theta(3)));
  theta(12) = 1.0 / (1 + exp(-theta(12)));
  return theta;
}
