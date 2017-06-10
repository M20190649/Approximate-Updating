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
vec FF(vec y, double phi, double gamma, double sigmaSqY, double sigmaSqX, double xTmean, double xTvar){
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
rowvec FFBS(vec y, rowvec theta){
  int T = y.size();
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T+3); //0 - T -> x0 - xT
  double at = theta[3];  //E(x1 | E(x0) = gamma)
  double pt = pow(theta[2], 2) * theta[1]  +  theta[1];
  double vt = y[0]  -  at;
  att[0] = at  +  pt * vt / (pt + theta[0]); //E(x1 | y1, x0)
  ptt[0] = pt  -  pow(pt, 2) / (pt + theta[0]);
  for(int t = 1; t < T; ++t){
    at = theta[3]  +  theta[2] * (att[t-1] - theta[3]); 
    pt = pow(theta[2], 2) * ptt[t-1]  +  theta[1];
    vt = y[t]  -  at;
    att[t] = at  +  pt * vt / (pt + theta[0]);
    ptt[t] = pt  -  pow(pt, 2) / (pt + theta[0]);
  }   
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T);
  vec ptT(T);
  draws[T] = att[T-1]  +  sqrt(ptt[T-1]) * randn<vec>(1)[0];
  for(int t = T - 1; t > 0; --t){
    vstar = draws[t+1]  -  theta[3]  -  theta[2] * (att[t-1] - theta[3]);
    fstar = pow(theta[2], 2) * ptt[t-1]  +  theta[1];
    mstar = ptt[t-1] * theta[2];
    atT[t] = att[t-1]  +  mstar * vstar / fstar;
    ptT[t] = ptt[t-1]  -  pow(mstar, 2) / fstar;
    draws[t] = atT[t]  +  sqrt(ptT[t]) * randn<vec>(1)[0];
  }
  atT[0] = theta[2] * draws[1]  -  theta[2] * theta[3] + pow(theta[2], 2) * theta[3]  +  1  -  pow(theta[2], 2);
  ptT[0] = theta[1];
  draws[0] = atT[0]  +  sqrt(ptT[0]) * randn<vec>(1)[0];
  draws[T+1] = att[T-1];
  draws[T+2] = ptt[T-1];
  return draws;
}

double logPhiDens(rowvec x, double sigmaSq, double phi, double gamma){
  int T = x.n_elem-1;
  
  double prior = 19.0 * log(1 + phi)  +  0.5 * log(1 - phi);
  double term1 = 0.5 * log(1-pow(phi, 2));
  double term2 = -1.0 / (2*sigmaSq);
  double term3 = (1 - pow(phi, 2)) * pow(x[0] - gamma, 2);
  double term4 = 0;
  for(int t = 1; t < T+1; ++t){
    term4 += pow(x[t]-gamma-phi*(x[t-1]-gamma), 2);
  }
  return prior + term1 + term2 * (term3 + term4);
}

double logTruncDens (double x, double mu, double sd){
  boost::math::normal_distribution<> aux (mu, sd);
  double numer = log(pdf(aux, x));
  double denom = log(sd) + log(cdf(aux, 1) - cdf(aux, -1));
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
    double canPDens = logPhiDens (x(i-1, span(0, T)), theta(i, 1), candidate, theta(i-1, 3));
    double oldPDens = logPhiDens (x(i-1, span(0, T)), theta(i, 1), theta(i-1, 2), theta(i-1, 3));
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
    x.row(i) = FFBS(y, theta.row(i));
  }
  Rcpp::Rcout << accept / reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("x") = x);
}

double logPhiDens_AR1 (vec x, double phi, double sigmaSq){
  double prior =  19.0 * log(1 + phi)  +  0.5 * log(1 - phi);
  double term1 = 0.5 * log(1-pow(phi, 2));
  double term2 = -1.0 / (2*sigmaSq);
  double term3 = - pow(phi*x[0], 2);
  double term4 = 0;
  for(int t = 1; t < x.n_elem; ++t){
    term3 += pow(phi * x[t-1], 2);
    term4 -= 2.0* phi  *  x[t]*x[t-1];
  }
  return prior + term1 + term2 * (term3 + term4);
}

// [[Rcpp::export]]
mat MCMC_AR1 (vec x, int reps = 10000){
  int T = x.n_elem - 1;
  mat draws(reps, 2, fill::zeros);
  int accept = 0;
  double sigShape = 0.5*(T+3);
  
  for(int i = 1; i < reps; ++i){
      
    double sigScale = 1 + 0.5 * (1- pow(draws(i-1, 1)*x[0], 2));
    for(int t = 1; t <= T; ++t){
      sigScale += 0.5 * pow(x[t] - draws(i-1, 1)*x[t-1], 2);
    }
    draws(i, 0) = 1.0 / randg<vec>(1, distr_param(sigShape, 1.0/sigScale))[0];
   
    double candidate;
    bool flag = true;
    while(flag){
      candidate = draws(i-1, 1) + 0.1*randn<vec>(1)[0];
      if(candidate < 1 & candidate > -1){
        flag = false;
      }
    }
    double canQDens = logTruncDens(candidate, draws(i-1, 1), 0.1);
    double oldQDens = logTruncDens(draws(i-1, 1), candidate, 0.1);
    double canPDens = logPhiDens_AR1(x, candidate, draws(i, 0));
    double oldPDens = logPhiDens_AR1(x, draws(i-1, 1), draws(i, 0));
    double ratio = exp(canPDens - oldPDens - canQDens + oldQDens);
    if(randu<vec>(1)[0] < ratio){
      draws(i, 1) = candidate;
      accept += 1;
    } else {
      draws(i, 1) = draws(i-1, 1);
    }
  }
  Rcpp::Rcout << accept * 1.0 / (reps * 1.0) << std::endl;
  return draws;
}

// [[Rcpp::export]]
mat MCMC_AR1M (vec x, int reps=10000){
  int T = x.n_elem -1;
  mat theta(reps, 3);
  rowvec initial = {1, 0, 0};
  theta.row(0) = initial;
  
  double sigmaSqShape = 0.5*T + 1.5;
  double accept = 0;
  
  for(int i = 1; i < reps; ++i){
    //sigmaSqX ~ IG(shape, scale)
    double sigmaSqScale = 1 + (1-pow(theta(i-1, 1), 2)) * pow(x(0) - theta(i-1, 2), 2) / 2;
    for(int t = 1; t < T+1; ++t){
      sigmaSqScale += pow(x(t) - theta(i-1, 2) - theta(i-1, 1)*(x(t-1) - theta(i-1, 2)), 2) / 2;
    }
    theta(i, 0) = (1.0 / randg<vec>(1, distr_param(sigmaSqShape, 1.0/sigmaSqScale)))[0];
    
    //phi ~ MH with RW Truncnorm candidate
    
    bool flag = true;
    double candidate;
    while(flag){
      candidate = theta(i-1, 1) + 0.1 * randn<vec>(1)[0];
      if(candidate < 1 & candidate > -1){
        flag = false;
      }
    }
    double canQDens = logTruncDens (candidate, theta(i-1, 1), 0.1);
    double oldQDens = logTruncDens (theta(i-1, 1), candidate, 0.1);
    double canPDens = logPhiDens (x.t(), theta(i, 0), candidate, theta(i-1, 2));
    double oldPDens = logPhiDens (x.t(), theta(i, 0), theta(i-1, 1), theta(i-1, 2));
    double ratio = exp(canPDens - oldPDens - canQDens + oldQDens);
    if (randu<vec>(1)[0] < ratio){
      theta(i, 1) = candidate;
      accept += 1;
    } else {
      theta(i, 1) = theta(i-1, 1);
    }
    //gamma ~ Normal(mean, var)
    double denom = theta(i, 0) + 100 * (T + T*pow(theta(i, 1), 2) - 2*T*theta(i, 1) + 1 - pow(theta(i, 1), 2));
    double meanNumer = (1 - pow(theta(i, 1), 2)) * x(0);
    double varNumer = 100 * theta(i, 0);
    for(int t = 1; t < T+1; ++t){
      meanNumer += x(t) + pow(theta(i, 1), 2)*x(t-1) - theta(i, 1) * (x(t-1) + x(t)); 
    }
    meanNumer *= 100;
    theta(i, 2) = meanNumer / denom + sqrt(varNumer / denom) * randn<vec>(1)[0];
  }
  Rcpp::Rcout << accept / reps << std::endl;
  return theta;
}

double zdens(double x){
  return sqrt(2.0 * 3.14159) * exp(-pow(x, 2) / 2);
}

rowvec svmFFBS(vec y, double sigmaSq, double phi, double gamma, double x0, rowvec s, vec mixmeans, vec mixvar){
  int T = y.size();
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T+2); //0 - (T-1) -> x1- xT, add xT mean and var to the end.
  
  double at = gamma  +  phi * (x0 - gamma);  //E(x1 | x0) = gamma + phi(x0 - gamma)
  double pt = sigmaSq; //Var(x1 | x0) = sigsq
  double vt = y[0]  +  1.2704  -  at  -  mixmeans(s(0));
  att[0] = at  +  pt * vt / (pt + mixvar(s(0))); //E(x1 | y1, x0)
  ptt[0] = pt  -  pow(pt, 2) / (pt + mixvar(s(0)));
  for(int t = 1; t < T; ++t){
    at = gamma  +  phi * (att[t-1] - gamma);    //E(x_t+1 | y_t, xt)
    pt = pow(phi, 2)*ptt[t-1] + sigmaSq;
    vt = y[t]  +  1.2704  -  at  -  mixmeans(s(t));
    att[t] = at  +  pt * vt / (pt + mixvar(s(t)));     //E(x_t+1 | y_t+1, x_t)
    ptt[t] = pt  -  pow(pt, 2) / (pt + mixvar(s(t)));
  }   
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T-1); // 0 - (T-2) -> x1 - x(t-1)
  vec ptT(T-1);
  draws[T-1] = att[T-1]  +  sqrt(ptt[T-1]) * randn<vec>(1)[0];
  for(int t = T - 2; t >= 0; --t){
    vstar = draws[t+1]  -  gamma  -  phi * (att[t] - gamma);
    fstar = pow(phi, 2) * ptt[t] + sigmaSq;
    mstar = ptt[t] * phi;
    atT[t] = att[t]  +  mstar * vstar / fstar;
    ptT[t] = ptt[t]  -  pow(mstar, 2) / fstar;
    draws[t] = atT[t]  +  sqrt(ptT[t]) * randn<vec>(1)[0];
  }
  draws[T] = att[T-1];
  draws[T+1] = ptt[T-1];

  return draws;
  
}

// [[Rcpp::export]]
Rcpp::List SVM_MCMC(vec y, int reps){
  // initialising things {
  int T = y.n_elem;
  mat theta(reps, 3);
  mat x = randn<mat>(reps, T+3);
  mat s (reps, T);
  vec ps(7);
  
  vec mixmeans = {-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -0.108819};
  vec mixvar = {5.79587, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261};
  vec mixweights = {0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.2575};
  
  vec sumweights = mixweights;
  for(int i = 1; i < 7; ++i){
    sumweights(i) += sumweights(i-1);
  }
  for(int t = 0; t < T; ++t){
    double u = randu<vec>(1)[0];
    if(u < sumweights(0)) s(0, t) = 0;
    else if(u < sumweights(1)) s(0, t) = 1;
    else if(u < sumweights(2)) s(0, t) = 2;
    else if(u < sumweights(3)) s(0, t) = 3;
    else if(u < sumweights(4)) s(0, t) = 4;
    else if(u < sumweights(5)) s(0, t) = 5;
    else s(0, t) = 6;
  }
  
  rowvec initial = {0.2, 0.9, 0};
  theta.row(0) = initial;
  
  double acceptPhi = 0;
  double sigmaSqShape = 0.5 * T  +  3.0;
  //}
  for(int i = 1; i < reps; ++i){
    //x0 - xT {
    //x0 ~ Normal
    double x0mean = theta(i-1, 1) * x(i-1, 1)  -  theta(i-1, 1) * theta(i-1, 2)  +  pow(theta(i-1, 1), 2) * theta(i-1, 2)  +  1.0  -  pow(theta(i-1, 1), 2);
    x(i, 0) = x0mean  +  sqrt(theta(i-1, 0)) * randn<vec>(1)[0];
    
    //x1 - xT ~ FFBS MH 
    x(i, span(1, T+2)) = svmFFBS(y, theta(i-1, 0), theta(i-1, 1), theta(i-1, 2), x(i, 0), s.row(i-1), mixmeans, mixvar);
    //}
    
    //s ~ multinomial {
    for(int t = 0; t < T; ++t){
      for(int j = 0; j < 7; ++j){
        ps(j) = mixweights(j) / sqrt(mixvar(j)) * zdens((y[t] - x(i, t+1) + 1.2704 - mixmeans(j)) / sqrt(mixvar(j)));
        if(j == 0){
          sumweights(j) = ps(j);
        } else {
          sumweights(j) = sumweights(j-1)  +  ps(j);
        }
      }
      double u = sumweights(6) * randu<vec>(1)[0];
      if(u < sumweights(0)) s(i, t) = 0;
      else if(u < sumweights(1)) s(i, t) = 1;
      else if(u < sumweights(2)) s(i, t) = 2;
      else if(u < sumweights(3)) s(i, t) = 3;
      else if(u < sumweights(4)) s(i, t) = 4;
      else if(u < sumweights(5)) s(i, t) = 5;
      else s(i, t) = 6;
    }
    // }
    
    // sigmaSq ~ Inverse Gamma {
    double sigmaSqScale = 0.025  +  (1 - pow(theta(i-1, 1), 2)) * pow(x(i, 0) - theta(i-1, 2), 2) / 2;
    for(int t = 1; t < T+1; ++t){
      sigmaSqScale += pow(x(i, t) - theta(i-1, 2) - theta(i-1, 1)*(x(i, t-1)-theta(i-1, 2)), 2) / 2;
    }
    theta(i, 0) = (1.0 / randg<vec>(1, distr_param(sigmaSqShape, 1.0/sigmaSqScale)))[0];
    //}

    //phi ~ MH with RW Truncnorm candidate {
    bool flag = true;
    double candidate;
    while(flag){
      candidate = theta(i-1, 1)  +  0.1 * randn<vec>(1)[0];
      if(candidate < 1 & candidate > -1){
        flag = false;
      }
    }
    double canQDens = logTruncDens (candidate, theta(i-1, 1), 0.1);
    double oldQDens = logTruncDens (theta(i-1, 1), candidate, 0.1);
    double canPDens = logPhiDens (x(i, span(0, T)), theta(i, 0), candidate, theta(i-1, 2));
    double oldPDens = logPhiDens (x(i, span(0, T)), theta(i, 0), theta(i-1, 1), theta(i-1, 2));
    double ratio = exp(canPDens - oldPDens - canQDens + oldQDens);
    if (randu<vec>(1)[0] < ratio){
      theta(i, 1) = candidate;
      acceptPhi += 1;
    } else {
      theta(i, 1) = theta(i-1, 1);
    } 
    // }
    
    //gamma ~ Normal(mean, var) {
    double denom = theta(i, 0)  +  100 * (T + T*pow(theta(i, 1), 2) - 2*T*theta(i, 1) + 1 - pow(theta(i, 1), 2));
    double meanNumer = (1 - pow(theta(i, 1), 2)) * x(i, 0);
    double varNumer = 100 * theta(i, 0);
    for(int t = 1; t < T+1; ++t){
      meanNumer += x(i, t)  +  pow(theta(i, 1), 2) * x(i, t-1)  -  theta(i, 1) * (x(i, t-1)  +  x(i, t)); 
    }
    meanNumer *= 100;
    theta(i, 2) = meanNumer / denom  +  sqrt(varNumer / denom) * randn<vec>(1)[0];
    // }
  }
  
  Rcpp::Rcout << acceptPhi / reps << std::endl;
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("x") = x, Rcpp::Named("s") = s);
}

// [[Rcpp::export]]
vec thetaWeights (vec y, mat theta, mat x, mat s){
  int T = y.n_elem;
  int reps = theta.n_rows;
  vec trueY (T);
  vec w(reps);
  double wSum = 0;
  vec mixmeans = {-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -0.108819};
  vec mixvar = {5.79587, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261};
  vec mixweights = {0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.2575};
  
  for(int t = 0; t < T; ++t){
    trueY(t) = sqrt(exp(y(t)));
  }
  
  for(int i = 0; i < reps; ++i){
    double logF = 0;
    double logK = 0;
    for(int t = 0; t < T; ++t){
      double sd = sqrt(exp(x(i, t+1)));
      logF += log( zdens(trueY(t)/sd) / sd );
      double kComp = 0;
      for(int j = 0; j < 7; ++j){
        kComp += mixweights(j) / sqrt(mixvar(j)) * zdens((y(t) - x(i, t+1) + 1.2704 - mixmeans(j)) / sqrt(mixvar(j)));
      }
      logK += log(kComp);
    }
    w(i) = logF - logK;
    wSum += exp(w(i));
  }
  
  vec weights(reps);
  for(int i = 0; i < reps; ++i){
    weights(i) = exp(w(i)) / wSum;
  }
  return weights;
}