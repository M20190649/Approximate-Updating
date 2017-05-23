// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/gamma.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost::math;

double PLogDens (vec y, double sigmaSq) {
  int n = y.n_elem;
  double density = -2 * log(sigmaSq) - 1 / sigmaSq;
  for(int i = 0; i < n; ++i){
    density += 1.0 / sqrt(2*3.14159*sigmaSq) - pow(y[i], 2) / (2*sigmaSq);
  }
  return density;
}

double qLogDens (double epsilon, double mu, double delta) {
  double density = -log(2*3.14159)*0.5 - pow(epsilon, 2)*0.5 - delta - mu - exp(delta)*epsilon;
  return density;
}

double ELBO (vec y, double mu, double delta, int n = 100){
  double elbo = 0;
  vec epsilon = randn<vec>(100);
  vec sigmaSq = exp(mu + delta*epsilon);
  for(int i = 0; i < n; ++i){
    elbo += PLogDens(y, sigmaSq[i]) - qLogDens(epsilon[i], mu, delta);
  }
  return elbo/n;
}

vec deriv (vec y, double mu, double delta, double epsilon) {
  int N = y.n_elem;
  double sigmaSq = exp(mu + exp(delta)*epsilon);
  double dpdsigmaSq = -(N+2.0) / sigmaSq + 1/pow(sigmaSq, 2);
  for(int i = 0; i < N; ++i){
    dpdsigmaSq += pow(y[i], 2) / (2 * pow(sigmaSq, 2));
  }
  
  double dqdmu = -1;
  double dqddelta = - epsilon*exp(delta) - 1;
  double dsdmu = sigmaSq;
  double dsddelta = epsilon*sigmaSq*exp(delta);
  
  vec derivs = {dpdsigmaSq * dsdmu - dqdmu, dpdsigmaSq * dsddelta - dqddelta};
  return derivs;
}

//[[Rcpp::export]]
Rcpp::List SGA_Var_LN (vec y, int M, int maxIter, double alpha, double initMu, double initDelta, bool adam = true) {
  int T = y.n_elem;
  double beta1 = 0.9;
  double beta2 = 0.999;
  double e = 0.0000001;
  
  vec params = {initMu, initDelta};
  vec Mt(2, fill::zeros);
  vec Vt(2, fill::zeros);
  int iter = 0;
  vec LB(maxIter+1);
  LB[0] = ELBO(y, params[0], params[1]);
  double diff = 1;
  while(diff > 0.01 | iter < 10){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    vec partials(2, fill::zeros);
    vec eps = randn<vec>(M);
    for(int i = 0; i < M; ++i){
      partials += deriv(y, params[0], params[1], eps[i])/M;
    }
    if(adam){
      Mt = beta1*Mt + (1-beta1)*partials;
      Vt = beta2*Vt + (1-beta2)*pow(partials, 2);
      vec MtHat = Mt / (1 - pow(beta1, iter));
      vec VtHat = Vt / (1 - pow(beta2, iter));
      params += alpha * MtHat / (sqrt(VtHat) + e);
    } else {
      for(int j = 0; j < 2; ++j){
        Mt[j] += pow(partials[j],2); //Each GT is the sum of the squared partial diffs. 
        Vt[j] = alpha * pow(Mt[j], -0.5); //transform to get pt value
        if(iter > 0){ //first step is always of size eta so we skip it
          params[j] += Vt[j] * partials[j];
        }
      }
    }
    LB[iter] = ELBO(y, params[0], params[1]);
    diff = abs(LB[iter] - LB[iter-1]);
  }
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  return Rcpp::List::create(Rcpp::Named("Params") = params,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}

double qLogDensIG (double sigmaSq, double alpha, double beta) {
  double density = alpha * log(beta) - lgamma(alpha) - (alpha+1)*log(sigmaSq) - beta / sigmaSq;
  return density;
}

double qigammac (double x, double alpha, double beta){
  boost::math::gamma_distribution<> G(alpha, 1.0/beta);
  return 1.0 / quantile(G, 1-x);
}

double ELBOIG (vec y, double alpha, double beta, int n = 100){
  double elbo = 0;
  vec sigmaSq = 1/randg<vec>(100, distr_param(alpha, 1.0/beta));
  for(int i = 0; i < n; ++i){
    elbo += PLogDens(y, sigmaSq[i]) - qLogDensIG(sigmaSq[i], alpha, beta);
  }
  return elbo/n;
}

vec derivIG (vec y, double alpha, double beta, double eps) {
  double sigmaSq = qigammac(eps, alpha, beta);
  int N = y.n_elem;
  double dpds = -(N+2.0) / sigmaSq + 1/pow(sigmaSq, 2);
  for(int i = 0; i < N; ++i){
    dpds += pow(y[i], 2) / (2 * pow(sigmaSq, 2));
  }
  double h = 0.000001;
  double alpha2 = alpha + h;
  double beta2 = beta + h;
  double sigmaSqA = qigammac(eps, alpha2, beta);
  double sigmaSqB = qigammac(eps, alpha, beta2);
  double dsda = (sigmaSqA - sigmaSq)/h;
  double dsdb = (sigmaSqB - sigmaSq)/h;
  double dqda = (qLogDensIG(sigmaSqA, alpha2, beta) - qLogDensIG(sigmaSq, alpha, beta))/h;
  double dqdb = (qLogDensIG(sigmaSqB, alpha, beta2) - qLogDensIG(sigmaSq, alpha, beta))/h;
  vec derivs = {dpds * dsda - dqda, dpds * dsdb - dqdb};
  return derivs;
}


//[[Rcpp::export]]
Rcpp::List SGA_Var_IG (vec y, int M, int maxIter, double alpha, double initAlpha, double initBeta, bool adam = true) {
  int T = y.n_elem;
  double beta1 = 0.9;
  double beta2 = 0.999;
  double e = 0.0000001;
  
  vec params = {initAlpha, initBeta};
  vec Mt(2, fill::zeros);
  vec Vt(2, fill::zeros);
  int iter = 0;
  vec LB(maxIter+1);
  LB[0] = ELBOIG(y, params[0], params[1]);
  double diff = 1;
  while(diff > 0.01 | iter < 10){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    vec partials(2, fill::zeros);
    vec eps = randu<vec>(100);
    for(int i = 0; i < M; ++i){
      partials += derivIG(y, params[0], params[1], eps[i])/M;
    }
    if(adam){
      Mt = beta1*Mt + (1-beta1)*partials;
      Vt = beta2*Vt + (1-beta2)*pow(partials, 2);
      vec MtHat = Mt / (1 - pow(beta1, iter));
      vec VtHat = Vt / (1 - pow(beta2, iter));
      params += alpha * MtHat / (sqrt(VtHat) + e);
    } else {
      for(int j = 0; j < 2; ++j){
        Mt[j] += pow(partials[j],2); //Each GT is the sum of the squared partial diffs. 
        Vt[j] = alpha * pow(Mt[j], -0.5); //transform to get pt value
        if(iter > 0){ //first step is always of size eta so we skip it
          params[j] += Vt[j] * partials[j];
        }
      }
    }
    LB[iter] = ELBOIG(y, params[0], params[1]);
    diff = abs(LB[iter] - LB[iter-1]);
  }
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  return Rcpp::List::create(Rcpp::Named("Params") = params,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}