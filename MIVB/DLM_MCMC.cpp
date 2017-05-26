// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
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
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T+3); //0 - T -> x0 - xT
  double at = theta[3];  //E(x1 | E(x0) = gamma)
  double pt = pow(theta[2], 2)*theta[1] + theta[1];
  double vt = y[0];
  att[0] = pt * vt / (pt + theta[0]); //E(x1 | y1, x0)
  ptt[0] = pt - pow(pt, 2) / (pt + theta[0]);
  for(int t = 1; t < T; ++t){
    at = theta[3] + theta[2]*(att[t-1]-theta[3]); 
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

double logPhiDens(rowvec x, double sigmaSqX, double phi, double gamma){
  int T = x.n_elem-1;
  double term1 = 0.5 * log(1-pow(phi, 2));
  double term2 = - 1.0 / (2*sigmaSqX);
  double term3 = (1-pow(phi, 2)) * pow(x[0] - gamma, 2);
  double term4 = 0;
  for(int t = 1; t < T+1; ++t){
    term4 += pow(x[t] - gamma - phi * (x[t-1] - gamma), 2);
  }
  return term1 + term2 * (term3 + term4);
}

double logTruncDens (double x, double mu, double sd){
  boost::math::normal_distribution<> aux (mu, sd);
  double numer = log(pdf(aux, x));
  boost::math::normal_distribution<> min(mu, sd);
  boost::math::normal_distribution<> max(mu, sd);
  double denom = log(sd) + log(cdf(max, 1) - cdf(min, -1));
  return numer - denom;
}

// Main MCMC algorithm
// [[Rcpp::export]]
Rcpp::List DLM_MCMC(vec y, int reps){
  int T = y.n_elem;
  mat theta(reps, 4);
  mat x(reps, T+3, fill::randn);
  rowvec initial = {1, 1, 0, 0};
  theta.row(0) = initial;
  
  double sigmaSqYShape = 0.5*T + 1;
  double sigmaSqXShape = 0.5*T + 1.5;
  double accept = 0;

  for(int i = 1; i < reps; ++i){
    //sigmaSqY ~ IG(shape, scale)
    double sigmaSqYScale = 1;
    for(int t = 1; t < T+1; ++t){
      sigmaSqYScale += pow(y[t-1] - x(i-1, t), 2) / 2;
    }
    theta(i, 0) = (1.0 / randg<vec>(1, distr_param(sigmaSqYShape, 1.0/sigmaSqYScale)))[0];
    
    //sigmaSqX ~ IG(shape, scale)
    double sigmaSqXScale = 1 + (1-pow(theta(i-1, 2), 2)) * pow(x(i-1, 0) - theta(i-1, 3), 2) / 2;
    for(int t = 1; t < T+1; ++t){
      sigmaSqXScale += pow(x(i-1, t) - theta(i-1, 3) - theta(i-1, 2)*(x(i-1, t-1) - theta(i-1, 3)), 2) / 2;
    }
    theta(i, 1) = (1.0 / randg<vec>(1, distr_param(sigmaSqXShape, 1.0/sigmaSqXScale)))[0];
    
    //phi ~ MH with RW Truncnorm candidate
    bool flag = true;
    double candidate;
    while(flag){
      candidate = theta(i-1, 2) + 0.1 * randn<vec>(1)[0];
      if(candidate < 1 & candidate > -1){
        flag = false;
      }
    }
    double canQDens = logTruncDens (candidate, theta(i-1, 2), 0.1);
    double oldQDens = logTruncDens (theta(i-1, 2), candidate, 0.1);
    double canPDens = logPhiDens (x.row(i-1), theta(i, 1), candidate, theta(i-1, 3));
    double oldPDens = logPhiDens (x.row(i-1), theta(i, 1), theta(i-1, 2), theta(i-1, 3));
    double ratio = exp(canPDens - oldPDens - canQDens + oldQDens);
    if (randu<vec>(1)[0] < ratio){
      theta(i, 2) = candidate;
      accept += 1;
    } else {
      theta(i, 2) = theta(i-1, 2);
    }
    
    //gamma ~ Normal(mean, var)
    double denom = theta(i, 1) + 100 * (T + T*pow(theta(i, 2), 2) - 2*T*theta(i, 2) + 1 - pow(theta(i, 2), 2));
    double meanNumer = (1 - pow(theta(i, 2), 2)) * x(i-1, 0);
    double varNumer = 100 * theta(i, 1);
    for(int t = 1; t < T+1; ++t){
      meanNumer += x(i-1, t) + pow(theta(i, 2), 2)*x(i-1, t-1) - theta(i, 2) * (x(i-1, t-1) + x(i-1, t)); 
    }
    meanNumer *= 100;
    theta(i, 3) = meanNumer / denom + sqrt(varNumer / denom) * randn<vec>(1)[0];
    x.row(i) = FFBScpp(y, theta.row(i));
    if(i % 2000 == 0){
      Rcpp::Rcout << meanNumer / denom<< std::endl;
    }
  }
  Rcpp::Rcout << accept / reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("x") = x);
}
