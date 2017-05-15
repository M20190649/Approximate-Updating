// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// Used to filter states forward in forecasting with extra data added
// [[Rcpp::export]]
vec FFUpdatercpp(vec y, double phi, double gamma, double sigmaSqY, double sigmaSqX, double xTmean, double xTvar){
  int J = y.n_elem;
  y -= gamma;
  vec at(J);
  vec pt(J);
  vec att(J);
  vec ptt(J);
  at[0] = phi*xTmean;
  pt[0] = pow(phi, 2)*xTvar + sigmaSqX;
  double vt = y[0] - at[0];
  att[0] = at[0] + pt[0] * vt / (pt[0] + sigmaSqY);
  ptt[0] = pt[0] - pow(pt[0], 2) / (pt[0] + sigmaSqY);
  for(int t = 1; t < J; ++t){
    at[t] = phi*att[t-1];
    pt[t] = pow(phi, 2)*ptt[t-1] + sigmaSqX;
    vt = y[t] - at[t];
    att[t] = at[t] + pt[t] * vt / (pt[t] + sigmaSqY);
    ptt[t] = pt[t] - pow(pt[t], 2) / (pt[t] + sigmaSqY);
  }
  vec output = {att[J-1], ptt[J-1]};
  
  return output;
}

// Kalman FFBS for MCMC
// [[Rcpp::export]]
rowvec FFBScpp(vec y, rowvec theta){
  int T = y.size();
  //subtract the mean from y
  y -= theta[3];
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T+3); //0 - T -> x0 - xT
  double at;
  double pt = pow(theta[2], 2)*theta[1] + theta[1];
  double vt = y[0];
  att[0] = pt * vt / (pt + theta[0]);
  ptt[0] = pt - pow(pt, 2) / (pt + theta[0]);
  for(int t = 1; t < T; ++t){
    at = theta[2]*att[t-1];
    pt = pow(theta[2], 2)*ptt[t-1] + theta[1];
    vt = y[t] - at;
    att[t] = at + pt * vt / (pt + theta[0]);
    ptt[t] = pt - pow(pt, 2) / (pt + theta[0]);
  }   
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T);
  vec ptT(T);
  draws[T] = att[T-1] + sqrt(ptt[T-1])*randn<vec>(1)[0];
  for(int t = T - 1; t > 0; --t){
    vstar = draws[t+1] - theta[2]*att[t-1];
    fstar = pow(theta[2], 2)*ptt[t-1] + theta[1];
    mstar = ptt[t-1]*theta[2];
    atT[t] = att[t-1] + mstar * vstar / fstar;
    ptT[t] = ptt[t-1] - pow(mstar, 2) / fstar;
    draws[t] = atT[t] + sqrt(ptT[t])*randn<vec>(1)[0];
  }
  atT[0] = theta[2]*draws[1] / (pow(theta[2], 2) + 1);
  ptT[0] = theta[1] / (pow(theta[2], 2) + 1);
  draws[0] = atT[0] + sqrt(ptT[0]) * randn<vec>(1)[0];
  draws[T+1] = att[T-1];
  draws[T+2] = ptt[T-1];
  
  return draws;
}

// Main MCMC algorithm
// [[Rcpp::export]]
Rcpp::List DLM_MCMC(vec y, int reps){
  int T = y.n_elem;
  mat theta(reps, 4);
  mat x(reps, T+3, fill::randn);
  rowvec initial = {1, 1, 0, 0};
  theta.row(0) = initial;
  
  double sigmaSqYShape = T/2 + 1;
  double sigmaSqXShape = T/2 + 1.5;
  double sigmaSqYScale;
  double sigmaSqXScale;
  double meanPhiNumer;
  double meanPhiDenom;
  double meanGammaNumer;
  double meanGammaDenom;

  for(int i = 1; i < reps; ++i){
    //sigmaSqY ~ IG(shape, scale)
    sigmaSqYScale = 1;
    for(int t = 0; t < T; ++t){
      sigmaSqYScale += pow(y[t] - x(i-1, t+1) - theta(i-1, 3), 2) / 2;
    }
    theta(i, 0) = (1 / randg<vec>(1, distr_param(sigmaSqYShape, 1/sigmaSqYScale)))[0];
    
    //sigmaSqX ~ IG(shape, scale)
    sigmaSqXScale = 1 + pow(x(i, 0), 2) / 2;
    for(int t = 0; t < T; ++t){
      sigmaSqXScale += pow(x(i-1, t+1) - theta(i-1, 2)*x(i-1, t), 2) / 2;
    }
    theta(i, 1) = (1 / randg<vec>(1, distr_param(sigmaSqXShape, 1/sigmaSqXScale)))[0];
    
    //phi ~ TruncNormal(mean, var, low = -1, high = 1) - Not truncated right now, working fine as is.
    // Should implement it properly
    meanPhiNumer = meanPhiDenom = 0;
    for(int t = 0; t < T; ++t){
      meanPhiNumer += x(i-1, t) * x(i-1, t+1);
      meanPhiDenom += pow(x(i-1, t), 2);
    }
    theta(i, 2) = (meanPhiNumer/meanPhiDenom + sqrt(theta(i, 1)/meanPhiDenom) * randn<vec>(1))[0];
    
    //gamma ~ Normal(mean, var)
    meanGammaNumer = 0;
    meanGammaDenom = 10*T + theta(i, 0);
    for(int t = 0; t < T; ++t){
      meanGammaNumer += 10 * y[t] - x(i-1, t+1);
    }
    theta(i, 3) = (meanGammaNumer/meanGammaDenom + sqrt(10*theta(i, 0)/meanGammaDenom) * randn<vec>(1))[0];
    
    x.row(i) = FFBScpp(y, theta.row(i));
  }
  
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("x") = x);
}