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

struct logP {
  const MatrixXd& y;
  const int& obs;
  logP(const MatrixXd& yIn, const int& obsIn) :
   y(yIn), obs(obsIn) {}
  template <typename T>
  T operator ()(const Matrix<T, Dynamic, 1>& theta)
    const{
    using std::log; using std::pow; using std::exp; using std::sqrt;
    // Get elements of theta
    T sigSq = theta(0);
    T mu = theta(1);
    T phi = theta(2);
    T xT = theta(3);
    
    // Evaluate log(p(theta)
    T prior = -pow(mu, 2) / 200  - pow(phi - 0.5, 2) / (2 * 0.01)  -  2 * log(sigSq)  -  1.0 / sigSq;
    
    // Set up for kalman filtering
    Matrix<T, Dynamic, 1> at(obs), pt(obs), vt(obs), St(obs), Kt(obs), att(obs+1), ptt(obs+1);
    // att and ptt are offset by one as they also contain the values at t = 0, while at and pt are for t > 0 only.
    att(0) = mu; // mean t = 0 | 0
    ptt(0) = 3 * sigSq;  // var t = 0 | 0
    T ytDens; // log(p(yt | theta))
    T yAllDens = 0; // log(p(y_1:T | theta))
    T xTDens = 0; // log(p(x_T | theta, y_1:T))
    // kalman loop
    for(int t = 0; t < obs; ++t){
      at(t) = mu + phi * (att(t) - mu); // mean t | t-1
      pt(t) = pow(phi, 2) * ptt(t) + sigSq; // var t | t-1
      vt(t) = y(t-1) - at(t); // innovation residual
      St(t) = 1 + pt(t); // innovation variance
      Kt(t) = pt(t) / St(t); // kalman gain
      att(t+1) =  at(t) + Kt(t) * vt(t); // update mean t|t
      ptt(t+1) = pt(t) - Kt(t) * pt(t); // update variance t|t
      ytDens = - 0.5 * log(2*3.14159) -  0.5 * log(St(t))  -  pow(y(t-1) - pt(t), 2) / (2 * St(t)); 
      yAllDens += ytDens;
    }
    // log(p(xT | theta, y_1:T))
    xTDens = - 0.5 * log(2*3.14159)  -  0.5 * log(ptt(obs))  -  pow(xT - att(obs), 2) / (2 * ptt(obs));
  
    return yAllDens + xTDens + prior;
  }
};

struct logQ {
  const MatrixXd& epsilon;
  logQ(const MatrixXd& epsIn) :
    epsilon(epsIn) {}
  template <typename T>
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::pow; using std::exp; using stan::math::fabs;
    // Evaluate log(q(theta, xT))
    T qLogDens = - lambda(0) - lambda(4)*epsilon(0);
    for(int i = 0; i < 4; ++i){
      qLogDens -=  log(fabs(lambda((i+1)*4+i))) + 0.5 * pow(epsilon(i), 2);
    }
    return qLogDens;
  }
};

double evalLogQ (Matrix<double, Dynamic, 1> epsilon, MatrixXd lambda){
  double qLogDens = - lambda(0) - lambda(4)*epsilon(0);
  for(int i = 0; i < 4; ++i){
    qLogDens -=  log(fabs(lambda((i+1)*4+i))) + 0.5 * pow(epsilon(i), 2);
  }
  return qLogDens;
}

Matrix<double, Dynamic, 1> QDeriv(Matrix<double, Dynamic, 1> epsilon, MatrixXd lambda){
  Matrix<double, Dynamic, 1> output(20); output.fill(0);
  output(0) = -1;
  output(4) = -(epsilon(0) + 1.0/lambda(4));
  output(9) = -1.0 / lambda(9);
  output(14) = -1.0 / lambda(14);
  output(19) = -1.0 / lambda(19);
  return output;
}

Matrix<double, Dynamic, 1> FDeriv (Matrix<double, Dynamic, 1> epsilon, MatrixXd lambda){
  Matrix<double, Dynamic, 1> output(20); output.fill(0);
  output(0) = exp(lambda(0) + lambda(4)*epsilon(0));
  output(1) = output(2) = output(3) = 1;
  output(4) = epsilon(0) * output(0);
  output(8) = output(12) = output(16) = epsilon(0);
  output(9) = output(13) = output(17) = epsilon(1);
  output(14) = output(18) = epsilon(2);
  output(19) = epsilon(3);
  return output;
}

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
      // apply the transform of tilde(x_k) = x_k + a_k % 2, where a_k = 1 if rule_k > 0.5, 0 otherwise
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
Rcpp::List VBIL_KF (Rcpp::NumericMatrix yIn, Rcpp::NumericMatrix lambdaIn, int S, double thresholdIS = 0.8, 
                    double alpha = 0.1, double beta1 = 0.9, double beta2 = 0.999, double threshold=0.01, int maxIter=5000){
  // convert to Eigen format
  Map<MatrixXd> y(Rcpp::as<Map<MatrixXd> >(yIn));
  Map<MatrixXd> lambdaMat(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  Map<MatrixXd> lambda(lambdaMat.data(), 20, 1);
  // Initialise various components for the algorithm
  Matrix<double, Dynamic, 1> gradientQ(20), gradientF(20), meanGradient(20), meanGradientSq(20), 
  Mt(20), Vt(20), LB(maxIter+1), logPeval(S), gradP(4);
  MatrixXd theta(4, S), gradientP(20, S), epsilon(4, S);
  LB.fill(0); Mt.fill(0); Vt.fill(0);
  double logQeval;
  int T = y.rows();
  
  // Initial SOBOL QMC numbers
  mat sobol = sobol_points(S+100, 4);
  mat unif;
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

    // Importance Sampling 
    double checkIS = randu<vec>(1)[0];
    if(iter > 1 & checkIS < thresholdIS){
      for(int s = 0; s < S; ++s){
        Matrix<double, Dynamic, 1> impliedEpsilon(4);
        for(int i = 0; i < 4; ++i){
          if(i == 0){
            impliedEpsilon(i) = log(theta(i, s)) - lambda(i);
          } else {
            impliedEpsilon(i) = theta(i, s) - lambda(i);
          }
          for(int j = 0; j <=i; ++j){
            impliedEpsilon(i) -= lambda((i+1)*4+j) * epsilon(j);
          }
        }
        gradientF = FDeriv(impliedEpsilon, lambda);
        logQ q(epsilon.col(s));
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(q, lambda, logQeval, gradientQ);
     
        double weight = 0;
        for(int i = 0; i < 4; ++i){
          weight += pow(epsilon(i, s), 2) - pow(impliedEpsilon(i), 2);
        }
        weight = exp(weight/2);
   
        for(int i = 0; i < 20; ++i){
          meanGradient(i) += (gradientP(i, s) * gradientF(i) - gradientQ(i)) * weight / S;
          meanGradientSq(i) += pow(S * meanGradient(i), 2) / S;
        }
        LB(iter-1) += (logPeval(s) - logQeval) * weight / S;
       }
    } else {
      // Randomly shuffle the sobol numbers for RQMC
      unif = shuffle(sobol);
      // Create estimates of E(dLB/dlam) and E(dLB/dlam^2)
      for(int s = 0; s < S; ++s){
        // transform uniform numbers to standard normal
        for(int i = 0; i < 4; ++i){
          if(unif(s+100, i) > 0.999){
            epsilon(i, s) = 3;
          } else if (unif(s+100, i) < 0.001) {
            epsilon(i, s) = -3;
          } else {
            epsilon(i, s) = quantile(epsDist, unif(s+100, i));
          }
        }
        //  create thetas for use later in importance sampling
        for(int i = 0; i < 4; ++i){
          theta(i, s) = lambda(i);
          for(int j = 0; j <= i; ++j){
            theta(i, s) += lambda((i+1)*4+j) * epsilon(j);
          }
        }
        theta(0, s) = exp(theta(0, s));
        // take derivative of log joint
        logP p(y, T);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(p, theta.col(s), logPeval(s), gradP);
        // take derivative of q
        logQeval = evalLogQ(epsilon.col(s), lambda);
        gradientQ = QDeriv(epsilon.col(s), lambda);
        // take derivative of theta
        gradientF = FDeriv(epsilon.col(s), lambda);
        for(int i = 0; i < 20; ++i){
          // grab relevant element of dlogp / dtheta
          if(i < 4){
            gradientP(i, s) = gradP(i);
          } else {
            gradientP(i, s) = gradP(floor((i-4)/4));
          }
          meanGradient(i) += (gradientP(i, s) * gradientF(i) - gradientQ(i)) / S;
          meanGradientSq(i) += pow(S * meanGradient(i), 2) / S;
        }
        LB(iter-1) += (logPeval(s) - logQeval)/S;
      }
    }
    
    // adam updating rule
    for(int i = 0; i < 20; ++i){
      Mt(i) = Mt(i) * beta1  +  (1 - beta1) * meanGradient(i);
      Vt(i) = Vt(i) * beta2  +  (1 - beta2) * meanGradientSq(i);
      double MtHat = Mt(i) / (1 - pow(beta1, iter));
      double VtHat = Vt(i) / (1 - pow(beta2, iter));
      if(iter > 1){
      lambda(i) += alpha * MtHat / (sqrt(VtHat) + pow(1, -8));
      }
    }
    // check convergence
    if(iter > 5){
      omeanLB = meanLB;
      meanLB = 0.2 * (LB(iter-1) + LB(iter-2) + LB(iter-3) + LB(iter-4) + LB(iter-5));
      diff = std::fabs(meanLB - omeanLB);
    } 
    // report progress
    if(iter % 25 == 0){
      Rcpp::Rcout << "Iteration: " << iter << ", ELBO: " << LB(iter-1) << std::endl;
    }
  } // End while loop
  if(iter <= maxIter){
    // iter goes up by one before checking the maxIter condition, so need to use LB(iter-2)
    Rcpp::Rcout << "Converged after " << iter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  } else {
    Rcpp::Rcout << "Warning, failed to converge after " << maxIter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  }
  MatrixXd Mu = lambda.topRows(4);
  Map<MatrixXd> U(lambda.bottomRows(16).data(), 4, 4);
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("U") = U,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}


//MCMC stuff
// Kalman FFBS for MCMC
// [[Rcpp::export]]
rowvec FFBS(vec y, rowvec theta){
  int T = y.size();
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T+3); //0 - T -> x0 - xT
  double at = theta[2];  //E(x1 | E(x0) = gamma)
  double pt = pow(theta[1], 2) * theta[0]  +  theta[0];
  double vt = y[0]  -  at;
  att[0] = at  +  pt * vt / (pt + 1); //E(x1 | y1, x0)
  ptt[0] = pt  -  pow(pt, 2) / (pt + 1);
  for(int t = 1; t < T; ++t){
    at = theta[2]  +  theta[1] * (att[t-1] - theta[2]); 
    pt = pow(theta[1], 2) * ptt[t-1]  +  theta[0];
    vt = y[t]  -  at;
    att[t] = at  +  pt * vt / (pt + 1);
    ptt[t] = pt  -  pow(pt, 2) / (pt + 1);
  }   
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T);
  vec ptT(T);
  draws[T] = att[T-1]  +  sqrt(ptt[T-1]) * randn<vec>(1)[0];
  for(int t = T - 1; t > 0; --t){
    vstar = draws[t+1]  -  theta[2]  -  theta[1] * (att[t-1] - theta[2]);
    fstar = pow(theta[1], 2) * ptt[t-1]  +  theta[0];
    mstar = ptt[t-1] * theta[1];
    atT[t] = att[t-1]  +  mstar * vstar / fstar;
    ptT[t] = ptt[t-1]  -  pow(mstar, 2) / fstar;
    draws[t] = atT[t]  +  sqrt(ptT[t]) * randn<vec>(1)[0];
  }
  atT[0] = theta[1] * draws[1]  -  theta[1] * theta[2] + pow(theta[1], 2) * theta[2]  +  1  -  pow(theta[1], 2);
  ptT[0] = theta[0];
  draws[0] = atT[0]  +  sqrt(ptT[0]) * randn<vec>(1)[0];
  draws[T+1] = att[T-1];
  draws[T+2] = ptt[T-1];
  return draws;
}

double logPhiDens(rowvec x, double sigmaSq, double phi, double gamma){
  int T = x.n_elem-1;
  
  double prior = 19.0 * log(1 + phi)  +  0.5 * log(1 - phi);
  double term1 = 0.5 * log(1-pow(phi, 2));
  double term2 = -1.0 / (2*sigmaSq);
  double term3 = (1 - pow(phi, 2)) * pow(x[0] - gamma, 2);
  double term4 = 0;
  for(int t = 1; t < T+1; ++t){
    term4 += pow(x[t]-gamma-phi*(x[t-1]-gamma), 2);
  }
  return prior + term1 + term2 * (term3 + term4);
}

double logTruncDens (double x, double mu, double sd){
  boost::math::normal_distribution<> aux (mu, sd);
  double numer = log(pdf(aux, x));
  double denom = log(sd) + log(cdf(aux, 1) - cdf(aux, -1));
  return numer - denom;
}

// Main MCMC algorithm
// [[Rcpp::export]]
Rcpp::List DLM_MCMC(vec y, int reps){
  int T = y.n_elem;
  mat theta(reps, 3);
  mat x(reps, T+3, fill::randn);
  rowvec initial = {1, 0, 0};
  theta.row(0) = initial;
  
  double sigmaSqXShape = 0.5*T + 1.5;
  double accept = 0;
  
  for(int i = 1; i < reps; ++i){

    //sigmaSqX ~ IG(shape, scale)
    double sigmaSqXScale = 1 + (1-pow(theta(i-1, 1), 2)) * pow(x(i-1, 0) - theta(i-1, 2), 2) / 2;
    for(int t = 1; t < T+1; ++t){
      sigmaSqXScale += pow(x(i-1, t) - theta(i-1, 2) - theta(i-1, 1)*(x(i-1, t-1) - theta(i-1, 2)), 2) / 2;
    }
    theta(i, 0) = (1.0 / randg<vec>(1, distr_param(sigmaSqXShape, 1.0/sigmaSqXScale)))[0];
    
    //phi ~ MH with RW Truncnorm candidate
    bool flag = true;
    double candidate;
    while(flag){
      candidate = theta(i-1, 1) + 0.1 * randn<vec>(1)[0];
      if(candidate < 1 & candidate > -1){
        flag = false;
      }
    }
    double canQDens = logTruncDens (candidate, theta(i-1, 1), 0.1);
    double oldQDens = logTruncDens (theta(i-1, 1), candidate, 0.1);
    double canPDens = logPhiDens (x(i-1, span(0, T)), theta(i, 0), candidate, theta(i-1, 2));
    double oldPDens = logPhiDens (x(i-1, span(0, T)), theta(i, 0), theta(i-1, 2), theta(i-1, 2));
    double ratio = exp(canPDens - oldPDens - canQDens + oldQDens);
    if (randu<vec>(1)[0] < ratio){
      theta(i, 1) = candidate;
      accept += 1;
    } else {
      theta(i, 1) = theta(i-1, 1);
    }
    
    //mu ~ Normal(mean, var)
    double denom = theta(i, 0) + 100 * (T + T*pow(theta(i, 1), 2) - 2*T*theta(i, 1) + 1 - pow(theta(i, 1), 2));
    double meanNumer = (1 - pow(theta(i, 1), 2)) * x(i-1, 0);
    double varNumer = 100 * theta(i, 0);
    for(int t = 1; t < T+1; ++t){
      meanNumer += x(i-1, t) + pow(theta(i, 1), 2)*x(i-1, t-1) - theta(i, 1) * (x(i-1, t-1) + x(i-1, t)); 
    }
    meanNumer *= 100;
    theta(i, 2) = meanNumer / denom + sqrt(varNumer / denom) * randn<vec>(1)[0];
    x.row(i) = FFBS(y, theta.row(i));
  }
  Rcpp::Rcout << accept / reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("x") = x);
}

double logPhiDens_AR1 (vec x, double phi, double sigmaSq){
  double prior =  19.0 * log(1 + phi)  +  0.5 * log(1 - phi);
  double term1 = 0.5 * log(1-pow(phi, 2));
  double term2 = -1.0 / (2*sigmaSq);
  double term3 = - pow(phi*x[0], 2);
  double term4 = 0;
  for(int t = 1; t < x.n_elem; ++t){
    term3 += pow(phi * x[t-1], 2);
    term4 -= 2.0* phi  *  x[t]*x[t-1];
  }
  return prior + term1 + term2 * (term3 + term4);
}


