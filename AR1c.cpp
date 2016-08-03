// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include <vector>
#include <functional>   
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
vec genAR1(int n, rowvec theta) {
  double mu = theta[0];
  double rho = theta[1];
  double sigma = theta[2];
  vec eps(n);
  eps = randn<vec>(n);
  vec out(n);
  out[0] = mu + rho*(0-mu) + sigma*eps[0];
  for(int i = 1; i < n; ++i){
    out[i] = mu + rho*(out[i-1]-mu) + sigma*eps[i];
  }
  return out;
}

// [[Rcpp::export]]
vec gen_gamma(int n, int a, double scale) {
  //NumericVector unif;
  vec out(n);
  for(int i = 0; i < n; ++i){
    vec unif = randu<vec>(a);
    out[i] = -scale * sum(log(unif));
  }
  return out;
}

// [[Rcpp::export]]
arma::mat gen_mvrnorm(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
double autocor(vec x, int lag) {
  int T = x.size();  
  double out;
  for(int i = 0; i < T-lag; ++i) {
    out += x[i+lag]*x[i];
  }
  return out/T;
}

// [[Rcpp::export]]
vec sum_stats(vec x) { //Take three summary statistics of a vector
  vec out(3);
  out[0] = mean(x);
  out[1] = var(x);
  out[2] = autocor(x, 1);
  return out;
}

// [[Rcpp::export]]
cube transformsigma(vec x, mat Inv) { //Multiply each element of the Inverse matrix by each element of x
  int n = x.size();
  cube out(2,2,n);
  for(int i = 0; i < n; ++i){
    out(0,0,i) = Inv(0,0)*x[i];
    out(0,1,i) = Inv(0,1)*x[i];
    out(1,0,i) = Inv(1,0)*x[i];
    out(1,1,i) = Inv(1,1)*x[i];
  }
  return out;
}

// [[Rcpp::export]]
double eucdist(vec x, vec y) { //Euclidean distance
  double out = 0;
  for(int i = 0; i < 3; ++i){
    out += (x[i]-y[i])*(x[i]-y[i]);
  }
  return sqrt(out);
}

// [[Rcpp::export]]
vec ABC_Loop(int n, int T, mat theta, vec yss) { //Generate a model, take summary statistics, compare to yss.
  vec dist(n);
  for(int i = 0; i < n; ++i) {
    vec z(T);
    rowvec trow = theta.row(i);
    z = genAR1(T, trow); //Generate the model
    vec zss(3);
    zss = sum_stats(z); //Summary statistics of the new data
    dist[i] = eucdist(zss, yss); //comparision to y summary statistics
  }
  return dist;
}

// [[Rcpp::export]]
vec rao(mat chain, vec support, int M, double yt) {
  int N = support.size();
  vec out(N);
  double temp = 0;
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < M; ++j){
      temp += R::dnorm4(support(i), chain(j,0) + chain(j,1)*yt, chain(j,2), FALSE);
    }
    out(i) = temp/M;
    temp = 0;
  }
  return out;
}
