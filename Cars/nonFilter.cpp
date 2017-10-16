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

struct logP {
  const mat data;
  const int obs;
  const int N;
  const vec epsilon;
  const vec hyperParams;
  logP(const mat& dataIn, const int& obsIn, const int& nIn, const vec& epsIn, const vec& hypIn) :
    data(dataIn), obs(obsIn), N(nIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambdaV)
    const{
    using std::log; using std::exp; using std::pow; using std::cos; using std::sin;
    
    // Create theta
    Matrix<T, Dynamic, 1> theta(9 + 3 * (obs+1));
    for(int i = 0; i < 12 + 3*obs; ++i){
      theta(i) = lambdaV(i) + lambdaV(12+3*obs+i) * epsilon(i);
    }
    // Constrained Positive
    T sigSqA = exp(theta(0)), sigSqD = exp(theta(1)), alpha = exp(theta(4)), beta = exp(theta(5)), sigSqX = exp(theta(7)), sigSqY = exp(theta(8));
    // Constrained to (0, 1)
    T lambda = 1.0 / (1 + exp(-theta(2))), phi = 1.0 / (1 + exp(-theta(3)));
    // Unconstrained
    T gamma = theta(6);
    for(int i = 9 + (obs+1); i < 9 + 2*(obs+1); ++i){
      theta(i) = 1.0 / (1 + exp(-theta(i)));
    }
    double pi = 3.14159;
    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqA)  -  hyperParams(1) / sigSqA  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  +  
      (hyperParams(4) - 1) * log(lambda)  +  (hyperParams(5)-1) * log(1-lambda)  +
      (hyperParams(6) - 1) * log(phi)  +  (hyperParams(7)-1) * log(1-phi)  +  
      (hyperParams(8) - 1) * log(alpha)  -  hyperParams(9) * alpha  +  
      (hyperParams(10) - 1) * log(beta)  -  hyperParams(11) * beta  -
      pow(gamma - hyperParams(12), 2) / (2*hyperParams(13))  -
      (hyperParams(14) + 1) * log(sigSqX)  -  hyperParams(15) / sigSqX -
      (hyperParams(16) + 1) * log(sigSqY)  -  hyperParams(17) / sigSqY; 
    // Other variables
    Matrix<T, Dyanmic, 1> x(obs + 1), y(obs+1), v(obs+1);
    x(0) = 5;
    y(0) = 0;
    v(0) = 8;
  
    T density = 0;
    for(int t = 0; t < obs; ++t){
      // at transition
      density += -0.5 * log(2*pi*sigSqA) -pow(theta(10 + i) - phi * theta(9 + i), 2) / (2*sigSqA);
      // st transition
      density += (1 / (1 + 1/(pow(alpha, 1-theta(9+obs+i)) * pow(beta, theta(9+obs+i)) * abs(x(t) - 5)))-1) * log(theta(10+obs+i)) +
        (1 - 1 / (1 + 1/(pow(alpha, 1-theta(9+obs+i)) * pow(beta, theta(9+obs+i)) * abs(x(t) - 5))) - 1) * log(1-theta(10+obs+i));
      // dt transition
      density += -0.5 * log(2*pi*sigSqD) - pow(theta(11 + 2*obs +i) - pi/2 - theta(10+obs+i) * lamdda * (theta(10+2*obs+i) - pi/2) -
        (1 - theta(10+obs+i)) * gamma * (x(i) - 5), 2) / (2*sigSqD);
      // vt
      v(i+1) = v(i) + theta(10+i);
      // x and y
      x(i+1) = x(i) + v(i+1)*cos(theta(11+2*obs+i));
      y(i+1) = y(i) + v(i+1)*sin(theta(11+2*obs+i));
      // data
      density += -0.5 * log(2*sigSqX) - pow(data(i, 0) - x(i+1), 2) / (2*sigSqX);
      density += -0.5 * log(2*sigSqY) - pow(data(i, 0) - y(i+1), 2) / (2*sigSqY);
    }
    return density + prior;
  }
};