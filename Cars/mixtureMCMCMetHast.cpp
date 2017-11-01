
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double logDensity (mat data, vec theta, vec mu, mat varInv){
  int T = data.n_rows;
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  for(int t = 2; t < T; ++t){
    dens +=  - 0.5 * (theta(0) + theta(1)) - 
      pow(data(t, 0) - theta(2) * data(t-1, 0) - theta(3) * data(t-2, 0), 2) / (2 * exp(theta(0))) -
      pow(data(t, 1) - theta(4) * data(t-1, 1) - theta(5) * data(t-2, 1), 2) / (2 * exp(theta(1)));
  }
  return dens;
}
  
  
  