// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost::math;

// [[Rcpp::export]]
vec FFBS_C(vec y, vec theta){
  int T = y.size();
  //subtract the mean from y
  y -= theta[1];
  //forward filter steps
  vec att(T);
  vec ptt(T);
  vec draws(T+1);
  double at;
  double pt = pow(theta[0], 2)*theta[3] + theta[3];
  double vt = y[0];
  att[0] = pt * vt / (pt + theta[2]);
  ptt[0] = pt - pow(pt, 2) / (pt + theta[2]);
  for(int t = 1; t < T; ++t){
    at = theta[0]*att[t-1];
    pt = pow(theta[0], 2)*ptt[t-1] + theta[3];
    vt = y[t] - at;
    att[t] = pt * vt / (pt + theta[2]);
    ptt[t] = pt - pow(pt, 2) / (pt + theta[2]);
  }     

  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  double atT;
  double ptT;
  draws[T] = att[T] + sqrt(ptt[T])*randn<vec>(1)[0];
  for(int t = T - 1; t > 0; --t){
    vstar = draws[t+1] - theta[0]*att[t];
    fstar = pow(theta[0], 2)*ptt[t] + theta[3];
    mstar = ptt[t]*theta[0];
    atT = att[t] + mstar * vstar / fstar;
    ptT = ptt[t] - pow(mstar, 2) / fstar;
    draws[t] = atT + sqrt(ptT)*randn<vec>(1)[0];
  }
  double a0T = theta[0]*draws[1] / (pow(theta[0], 2) + 1);
  double p0T = theta[3] / (pow(theta[0], 2) + 1);
  draws[0] = a0T + sqrt(p0T) * randn<vec>(1)[0];
  return draws;
}