// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
vec fullsimul(mat L, vec mu, int T){
  vec epsilon = randn<vec>(T);
  vec output(T);
  for(int i = 0; i < T; ++i){
    output[i] = mu[i];
    for(int j = 0; j <= i; ++j){
      output[i] += L(i, j) * epsilon[j];
    }
  }
  return output;
}

// [[Rcpp::export]]
vec partialsimul(mat L, vec mu, int T){
  double xtsd = sqrt(L.submat(T-1, 4, T-1, T-1) * L.submat(T-1, 4, T-1, T-1).t())[0];
  vec epsilon = randn<vec>(5);
  vec output(5);
  for(int i = 0; i < 4; ++i){
    output[i] = mu[i];
    for(int j = 0; j <= i; ++j){
      output[i] += L(i, j) * epsilon[j];
    }
  }
  output[4] = mu[T] + L(T-1, 0)*epsilon[0] + L(T-1, 1)*epsilon [1] + L(T-1, 2)*epsilon[2] + L(T-1, 3)*epsilon[3] + xtsd*epsilon[4];
  return output;
}