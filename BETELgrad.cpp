
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List matrixCalculations(vec y, vec z, double alpha, double beta, mat g, mat hessian, vec lambdaHat, vec exponent){
  int n = y.n_elem;
  int d = g.n_cols;
  
  cube dgdt(d, d, n, fill::zeros);
  for(int i = 0; i < n; ++i){
    dgdt.slice(i) = { 
      {-1, -z(i), 0},
      {-z(i), -pow(z(i), 2), 0},
      {-3 * pow(y(i) - alpha - beta * z(i), 2), -3 * z(i) * pow(y(i) - alpha - beta * z(i), 2), -1}
    };
  }
  
  mat dh2dtlam(d, d, fill::zeros);
  for(int i = 0; i < n; ++i){
    dh2dtlam += 1.0 / n * (exponent(i) * (eye(d, d) + g.row(i).t() * lambdaHat.t()) * dgdt.slice(i));
  }
    
  mat dlamdt = - hessian.i() * dh2dtlam;
  mat productRule(d, n, fill::zeros);
  for(int i = 0; i < n; ++i){
    productRule.col(i) = dlamdt.t() * g.row(i).t() + dgdt.slice(i).t() * lambdaHat;
  }
  
  vec dpdt(d, fill::zeros);
  for(int j = 0; j < d; ++j){
    dpdt(j) = sum(productRule.row(j)) - n / sum(exponent) * as_scalar(productRule.row(j) * exponent);
  }
  double logp = sum(log(exponent)) - n * log(sum(exponent));
  return Rcpp::List::create(Rcpp::Named("grad") = dpdt,
                            Rcpp::Named("val") = logp);
  
}
