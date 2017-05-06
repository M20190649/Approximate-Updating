// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
rowvec FFBSsvm(vec y, rowvec theta, rowvec s){
  int T = y.size();
  vec mixmeans = {-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -0.108819};
  vec mixvar = {5.79586, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261};
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T); //0 - T -> x0 - xT
  double at = theta[1] / (1 - theta[2]);
  double pt = theta[0] / (1 - pow(theta[2], 2));
  double vt = y[0] - at - mixmeans[s[0]] + 1.2704;
  att[0] = at + pt * vt / (pt + mixvar[s[0]]);
  ptt[0] = pt - pow(pt, 2) / (pt + mixvar[s[0]]);
  for(int t = 1; t < T; ++t){
    at = theta[1] + theta[2]*att[t-1];
    pt = pow(theta[2], 2)*ptt[t-1] + theta[0];
    vt = y[t] - at - at - mixmeans[s[t]] + 1.2704;
    att[t] = at + pt * vt / (pt + mixvar[s[t]]);
    ptt[t] = pt - pow(pt, 2) / (pt + mixvar[s[t]]);
  }   
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T);
  vec ptT(T);
  draws[T-1] = att[T-1] + sqrt(ptt[T-1])*randn<vec>(1)[0];
  for(int t = T - 1; t > 0; --t){
    vstar = draws[t] - theta[2]*att[t-1];
    fstar = pow(theta[2], 2)*ptt[t-1] + theta[0];
    mstar = ptt[t-1]*theta[2];
    atT[t] = att[t-1] + mstar * vstar / fstar;
    ptT[t] = ptt[t-1] - pow(mstar, 2) / fstar;
    draws[t-1] = atT[t] + sqrt(ptT[t])*randn<vec>(1)[0];
  }
  return draws;
}

double zdens(double x){
  return sqrt(2 * 3.14159) * exp(-pow(x, 2) / 2);
}

// [[Rcpp::export]]
Rcpp::List SVM_MCMC(vec y, int reps){
  int T = y.n_elem;
  mat theta(reps, 3);
  mat a(reps, T, fill::randn);
  rowvec initial = {1, 1, 0.7};
  theta.row(0) = initial;
  mat s(reps, T);
  s.fill(4);
  vec ps(7);
  vec cumsum(7);
  
  vec mixmeans = {-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -0.108819};
  vec mixvar = {5.79586, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261};
  vec mixweights = {0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24556, 0.2575};

  double sigmaSqShape = 0.5*T + 1;
  double sigmaSqScale;
  double meanGammaNumer;
  double meanGammaDenom;
  double meanPhiNumer;
  double meanPhiDenom;
  
  for(int i = 1; i < reps; ++i){
    //alpha ~ FFBS
    a.row(i) = FFBSsvm(y, theta.row(i-1), s.row(i-1));
    
    //S ~ Multinomial
    for(int t = 0; t < T; ++t){
      for(int j = 0; j < 7; ++j){
        ps[j] = mixweights[j] / mixvar[j] * zdens((y[t] - a(i, t) + 1.2704 - mixmeans[j]) / mixvar[j]);
        if(j == 0){
          cumsum[j] = ps[j];
        } else {
          cumsum[j] = ps[j] + cumsum[j-1]; 
        }
      }
      double u = randu<vec>(1)[0] * cumsum[6];
      if(u < cumsum[0]){
        s(i, t) = 0;
      } else if(u < cumsum[1]){
        s(i, t) = 1;
      } else if(u < cumsum[2]){
        s(i, t) = 2;
      } else if(u < cumsum[3]){
        s(i, t) = 3;
      } else if(u < cumsum[4]){
        s(i, t) = 4;
      } else if(u < cumsum[5]){
        s(i, t) = 5;
      } else {
        s(i, t) = 6;
      }
    }
    
    //sigmaSq ~ IG(shape, scale)
    sigmaSqScale = 1;
    for(int t = 1; t < T; ++t){
      sigmaSqScale += pow(a(i,t) - theta(i-1, 2)*a(i, t-1) - theta(i-1, 1), 2) / 2;
    }
    theta(i, 0) = (1 / randg<vec>(1, distr_param(sigmaSqShape, 1/sigmaSqScale)))[0];
    
    //gamma ~ Normal(mean, var)
    meanGammaNumer = 0;
    meanGammaDenom = 10*(T-1) + theta(i, 0);
    for(int t = 1; t < T; ++t){
      meanGammaNumer += 10 * (a(i, t) - theta(i-1, 2) * a(i, t-1));
    }
    theta(i, 1) = (meanGammaNumer/meanGammaDenom + sqrt(10*theta(i, 0)/meanGammaDenom) * randn<vec>(1))[0];
    
    //phi ~ TruncNormal(mean, var, low = -1, high = 1)
    meanPhiNumer = 0;
    meanPhiDenom = 0;
    for(int t = 1; t < T; ++t){
      meanPhiNumer += a(i, t)*a(i, t-1) - theta(i, 1)*a(i, t-1);
      meanPhiDenom += pow(a(i, t-1), 2);
    }
    theta(i, 2) = (meanPhiNumer/meanPhiDenom + sqrt(theta(i, 0)/meanPhiDenom) * randn<vec>(1))[0];
  }
  
  return Rcpp::List::create(Rcpp::Named("theta") = theta, 
                            Rcpp::Named("alpha") = a,
                            Rcpp::Named("mixture") = s);
}