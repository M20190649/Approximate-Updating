// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include <vector>
#include <functional>   
#include <cmath>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double logpkernelfC(vec y_in, double mu_in, double sig_in){
  int n = y_in.size();
  double logkernel = (-51)*log(sig_in);
  vec terms(n);
  for(int i = 0; i < n; ++i){
    terms[i] = -0.5*(6)*log(1+(pow((y_in[i]-mu_in)/sig_in, 2.0)/3));
  }
  return(logkernel + sum(terms));
}

// [[Rcpp::export]]
double IGfC(double sigin, double vin, double sigsqhatin) {
  double out1 = 2/tgamma(vin/2);
  double out2 = pow(vin*sigsqhatin/2, vin/2);
  double out3 = exp(-vin*sigsqhatin/(2*pow(sigin, 2.0)));
  double out4 = pow(sigin, vin+1);
  return(out1*out2*out3/out4);
}


// [[Rcpp::export]]
double mu_MHfC (double mu_last, double sig_in, vec y_in) {
  double ybar = mean(y_in);
  double sy = sqrt(var(y_in));
  double mu_c = ybar + sy/sqrt(50)*randn<vec>(1)[1];
  double logalpha;
  logalpha = logpkernelfC(y_in, mu_c, sig_in);
  logalpha -= logpkernelfC(y_in, mu_last, sig_in);
  logalpha -= log(R::dnorm(mu_c, ybar, sy/sqrt(50), FALSE));
  logalpha += log(R::dnorm(mu_last, ybar, sy/sqrt(50), FALSE));
  double alpha = min(1.0, exp(logalpha));
  double u = randu<vec>(1)[0];
    if(u <= alpha)
      return mu_c;
    else
      return mu_last;
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
double sig_MHfC (double mu_in, double sig_last, vec y_in){
  double sighatsq = 0;
  int n = y_in.size();
  for(int i = 0; i < n; ++i) {
    sighatsq += pow(y_in[i]-mu_in, 2.0)/50;
  }
  double sig_c = 1/sqrt(gen_gamma(1, 25, 2/(50*sighatsq))[0]);
  double logalpha = logpkernelfC(y_in, mu_in, sig_c);
  logalpha -= logpkernelfC(y_in, mu_in, sig_last);
  logalpha -= log(IGfC(sig_c, 50, sighatsq));
  logalpha += log(IGfC(sig_last, 50, sighatsq));
  double alpha = min(1.0, exp(logalpha));
  double u = randu<vec>(1)[0];
  if(u <= alpha)
    return sig_c;
  else
    return sig_last;
}

// [[Rcpp::export]]
mat MHchain (int B, int NG, vec y){
  int iter = B+NG;
  mat MH(iter, 2);
  MH(0, 0) = mean(y);
  MH(0, 1) = sqrt(var(y));
  for(int i = 1; i < iter; ++i){
    MH(i, 0) = mu_MHfC(MH(i-1, 0), MH(i-1, 1), y);
    MH(i, 1) = sig_MHfC(MH(i-1, 0), MH(i-1, 1), y);
  }
  return MH;
}

