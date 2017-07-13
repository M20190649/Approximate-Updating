// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include <math.h>
#include "sobol.h"

// [[Rcpp::export]]
arma::mat testSobol(int M, int N){
  arma::mat test = sobol_points(M, N);
  return test;
}

// [[Rcpp::export]]
arma::mat shuffle(arma::mat sobol){
  int N = sobol.n_rows;
  int D = sobol.n_cols;
  arma::mat output(N, D, arma::fill::zeros);
  // draw a random rule of: switch 1 and 0  /  do not switch for each binary digit.
  arma::vec rule = arma::randu<arma::vec>(16);
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < D; ++j){
      // grab element of the sobol sequence
      double x = sobol(i, j);
      // convert to a binary representation
      arma::uvec binary(16, arma::fill::zeros);
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