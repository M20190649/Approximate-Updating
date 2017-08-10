// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
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
  rowvec draws(T+1); //0 - T -> x1 - xT+1
  double at = theta[1];
  double pt = theta[0] / (1 - pow(theta[2], 2));
  double vt = y[0]  -  at  -  mixmeans[s[0]]  +  1.2704;
  att[0] = at  +  pt * vt / (pt + mixvar[s[0]]);
  ptt[0] = pt  -  pow(pt, 2) / (pt + mixvar[s[0]]);
  for(int t = 1; t < T; ++t){
    at = theta[1]  +  theta[2] * (att[t-1] - theta[1]);
    pt = pow(theta[2], 2) * ptt[t-1]  +  theta[0];
    vt = y[t]  -  at  -  mixmeans[s[t]]  +  1.2704;
    att[t] = at  +  pt * vt / (pt + mixvar[s[t]]);
    ptt[t] = pt  -  pow(pt, 2) / (pt + mixvar[s[t]]);
  }   
  // x_T+1
  at = theta[1]  +  theta[2] * (att[T-1] - theta[1]);
  pt = pow(theta[2], 2) * ptt[T-1]  +  theta[0];
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T+1);
  vec ptT(T+1);
  draws[T] = at  +  sqrt(pt) * randn<vec>(1)[0];
  for(int t = T; t > 0; --t){
    vstar = draws[t]  -  theta[1]  -  theta[2] * (att[t-1] - theta[1]);
    fstar = pow(theta[2], 2) * ptt[t-1]  +  theta[0];
    mstar = ptt[t-1] * theta[2];
    atT[t] = att[t-1]  +  mstar * vstar / fstar;
    ptT[t] = ptt[t-1]  -  pow(mstar, 2) / fstar;
    draws[t-1] = atT[t]  +  sqrt(ptT[t]) * randn<vec>(1)[0];
  }
  return draws;
}

double zdens(double x){
  return sqrt(2 * 3.14159) * exp(-pow(x, 2) / 2);
}

// [[Rcpp::export]]
double phiDens(double phi, double sigSq, double mu, rowvec a){
  double density = 19.0 * log(1 + phi)  +  0.5 * log(1 - phi)  +  0.5 * log(1-pow(phi, 2))  -  (1 - pow(phi, 2)) * pow(a(0) - mu, 2) / (2 * sigSq);
  for(int t = 1; t < a.n_elem; ++t){
    density -= pow(a(t) - mu - phi * (a(t-1) - mu), 2) / (2 * sigSq);
  }
  return density;
}

double logTruncDens (double x, double mu, double sd){
  boost::math::normal_distribution<> aux (mu, sd);
  double numer = log(pdf(aux, x));
  double denom = log(sd) + log(cdf(aux, 1) - cdf(aux, -1));
  return numer - denom;
}

// [[Rcpp::export]]
Rcpp::List SVM_MCMC(vec y, int reps, double stepSize){
  int T = y.n_elem;
  mat theta(reps, 3), a(reps, T+1), s(reps, T);
  rowvec initial = {0.02, 0, 0.9};
  theta.row(0) = initial;
  s.fill(4);
  
  a(0 , 0) = initial(1)  +   sqrt(initial(0) / (1 - pow(initial(2), 2))) * randn<vec>(1)[0];
  for(int t = 1; t < T; ++t){
    a(0, t) = initial(1)  +  initial(2) * (a(0, t-1) - initial(1))  +  sqrt(initial(0)) * randn<vec>(1)[0];
  }
  
  vec ps(7), cumsum(7);
  
  vec mixmeans = {-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -0.108819};
  vec mixvar = {5.79586, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261};
  vec mixweights = {0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24556, 0.2575};
  
  double sigmaSqShape, sigmaSqScale, meanMuNumer, meanMuDenom, candidate, pOld, pCan, qOld, qCan, ratio, u;
  sigmaSqShape = 0.5 * T + 3.0;
  double accept = 0;
  
  for(int i = 1; i < reps; ++i){
     //sigmaSq ~ IG(shape, scale)
    sigmaSqScale = 0.025  +  (1 - pow(theta(i-1, 2), 2)) * pow(a(i-1, 0) - theta(i-1, 1), 2) / 2;
    for(int t = 1; t < T+1; ++t){
      sigmaSqScale += pow(a(i-1,t) - theta(i-1, 1) - theta(i-1, 2)*(a(i-1, t-1) - theta(i-1, 1)), 2) / 2;
    }
    if(sigmaSqScale < 0.001){
      sigmaSqScale = 0.001;
    }
    theta(i, 0) = 1.0 / randg<vec>(1, distr_param(sigmaSqShape, 1.0/sigmaSqScale))[0];
    
    //Mu ~ Normal(mean, var)
    meanMuNumer = a(i-1, 0) * (1 - pow(theta(i-1, 2), 2));
    meanMuDenom = theta(i, 0)  +  10 * (T + T*pow(theta(i-1, 2), 2)  - 2*T*theta(i-1, 2) + 1 - pow(theta(i-1, 2), 2));
    for(int t = 1; t < T+1; ++t){
      meanMuNumer += a(i-1, t)  +  pow(theta(i-1, 2), 2) * a(i-1, t-1)  -  theta(i-1, 2) * (a(i-1, t) + a(i-1, t-1)); 
    }
    meanMuNumer *= 10;
    theta(i, 1) = meanMuNumer/meanMuDenom + sqrt(10 * theta(i, 0) / meanMuDenom) * randn<vec>(1)[0];
    
    //phi ~ Truncated RWMH
    bool stationary = false;
    while(!stationary){
      candidate = theta(i-1, 2)  +  stepSize * randn<vec>(1)[0];
      if(candidate < 1 & candidate > -1){
        stationary = true;
      }
    }
    
    pCan = phiDens(candidate, theta(i, 0), theta(i, 1), a.row(i-1));
    pOld = phiDens(theta(i-1, 2), theta(i, 0), theta(i, 1), a.row(i-1));
    qCan = logTruncDens(candidate, theta(i-1, 2), stepSize);
    qOld = logTruncDens(theta(i-1, 2), candidate, stepSize);
    
    ratio = min(1.0, exp(pCan - pOld + qOld - qCan));
    u = randu<vec>(1)[0];
    if(u < ratio){
      accept += 1;
      theta(i, 2) = candidate;
    } else {
      theta(i, 2) = theta(i-1, 2);
    }
    
    //S ~ Multinomial
    for(int t = 0; t < T; ++t){
      for(int j = 0; j < 7; ++j){
        ps[j] = mixweights[j] / mixvar[j] * zdens((y[t] - a(i-1, t+1) + 1.2704 - mixmeans[j]) / mixvar[j]);
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
    
    //alpha ~ FFBS
    a.row(i) = FFBSsvm(y, theta.row(i), s.row(i));

  }
  Rcpp::Rcout << accept/reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta, 
                            Rcpp::Named("x") = a,
                            Rcpp::Named("mixture") = s);
}