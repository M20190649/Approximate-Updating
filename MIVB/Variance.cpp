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

//[[Rcpp::export]]
double pLogDens (vec y, double sigmaSq) {
  int n = y.n_elem;
  double density = -2 * log(sigmaSq) - 1.0 / sigmaSq;
  for(int i = 0; i < n; ++i){
    density +=  -log(sigmaSq)*0.5 - pow(y[i], 2) / (2*sigmaSq);
  }
  return density;
}

double qLogDens (double sigmaSq, double par1, double par2, bool lognormal) {
  double density;
  if(lognormal){
    density = - log(sigmaSq) - par2 - log(2*3.14159)*0.5 - pow(log(sigmaSq) - par1, 2) / (2*exp(2*par2));
  } else {
    par1 = exp(par1);
    par2 = exp(par2);
    density = par1 * log(par2) - lgamma(par1) - (par1+1)*log(sigmaSq) - par2 / sigmaSq;
  }
  return density;
}

double qigamma (double x, double alpha, double beta){
  boost::math::gamma_distribution<> G(alpha, 1.0/beta);
  return 1.0 / quantile(G, 1-x);
}

double ELBO (vec y, double par1, double par2, bool lognormal, int n = 100){
  double elbo = 0;
  vec epsilon;
  double sigmaSq;
  if(lognormal){
    epsilon = randn<vec>(n);
  } else {
    epsilon = randu<vec>(n);
  }
  for(int i = 0; i < n; ++i){
    if(lognormal){
      sigmaSq = exp(par1 + exp(par2)*epsilon[i]);
    } else {
      sigmaSq = qigamma(epsilon[i], par1, par2);
    }
    elbo += pLogDens(y, sigmaSq) - qLogDens(sigmaSq, par1, par2, lognormal);
  }
  return elbo/n;
}

vec deriv (vec y, double par1, double par2, double epsilon, bool lognormal) {
  int N = y.n_elem;
  if(lognormal){
    double sigmaSq = exp(par1 + exp(par2)*epsilon);
    double dldp1 = -N*0.5 - 1 + pow(sigmaSq, -1);
    double dldp2 = exp(par2)*epsilon * (-N*0.5 -1 + pow(sigmaSq, -1)) + 1;
    for(int i = 0; i < N; ++i){
      dldp1 += pow(y[i], 2) / (2 * sigmaSq);
      dldp2 += exp(par2) * epsilon * pow(y[i], 2) / (2 * sigmaSq);
    }
    vec derivs = {dldp1, dldp2};
    return derivs;
  } else {
    double sigmaSq = qigamma(epsilon, exp(par1), exp(par2));
    double h = 0.00000001;
    double par1h = par1 + h;
    double par2h = par2 + h;
    double sigmaSq1 = qigamma(epsilon, exp(par1h), exp(par2));
    double sigmaSq2 = qigamma(epsilon, exp(par1), exp(par2h));
    double dpds = -(0.5*N+2.0) * pow(sigmaSq, -1) + pow(sigmaSq, -2);
    for(int i = 0; i < N; ++i){
        dpds += pow(y[i], 2) * pow(sigmaSq, -2) * 0.5;
      }
    double dsdp1 = (sigmaSq1 - sigmaSq)/h;
    double dsdp2 = (sigmaSq2 - sigmaSq)/h;
    double dqdp1 = (qLogDens(sigmaSq1, par1h, par2, lognormal) - qLogDens(sigmaSq, par1, par2, lognormal))/h;
    double dqdp2 = (qLogDens(sigmaSq2, par1, par2h, lognormal) - qLogDens(sigmaSq, par1, par2, lognormal))/h;
    vec derivs = {dpds * dsdp1 - dqdp1, dpds * dsdp2 - dqdp2};
    return derivs;
  }
}

//[[Rcpp::export]]
Rcpp::List SGA_Var (vec y, bool lognormal, double initPar1, double initPar2, int M=50, int maxIter=5000, 
                    double threshold=0.01, double alpha=0.1, bool adam=true) {
  int T = y.n_elem;
  double beta1 = 0.9;
  double beta2 = 0.999;
  double e = 0.0000001;
  vec params = {initPar1, initPar2};
  vec Mt = zeros<vec>(2);
  vec Vt = zeros<vec>(2);
  vec partials = zeros<vec>(2);
  vec pSq = zeros<vec>(2);
  int iter = 0;
  vec LB(maxIter+1);
  LB[0] = ELBO(y, params[0], params[1], lognormal);
  double diff = 1;
  while(diff > threshold | iter < 10){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    partials.fill(0);
    pSq.fill(0);
    vec epsilon;
    if(lognormal){
      epsilon = randn<vec>(M);
    } else {
      epsilon = randu<vec>(M);
    }
    for(int i = 0; i < M; ++i){
      vec p = deriv(y, params[0], params[1], epsilon[i], lognormal);
      partials += p/M;
      pSq += pow(p, 2)/M;
    }
    if(adam){
      Mt = beta1*Mt + (1-beta1)*partials;
      Vt = beta2*Vt + (1-beta2)*pSq;
      vec MtHat = Mt / (1 - pow(beta1, iter));
      vec VtHat = Vt / (1 - pow(beta2, iter));
      if(iter > 1) {
        params += alpha * MtHat / (sqrt(VtHat) + e);
      }
    } else {
      Mt += pow(partials,2); 
      Vt = alpha * pow(Mt, -0.5);
      if(iter > 1){
        params += Vt % partials;
      }
    }
    LB[iter] = ELBO(y, params[0], params[1], lognormal);
    diff = abs(LB[iter] - LB[iter-1]);
  }
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  return Rcpp::List::create(Rcpp::Named("Params") = params,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}