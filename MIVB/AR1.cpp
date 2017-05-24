// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double pLogDens(vec x, double sigmaSq, double phi){
  int T = x.n_elem - 1;
  double density = -0.5*(T+5)*log(sigmaSq) + 0.5*log(1-pow(phi, 2)) - pow(sigmaSq, -1);
  density -= pow(x[0], 2) * (1 - pow(phi, 2)) / (2 * sigmaSq);
  for(int t = 1; t < T+1; ++t){
    density -= pow(x[t] - phi*x[t-1], 2) / (2*sigmaSq);
  }
  return density;
}

double qLogDens(vec epsilon, vec mu, mat L){
  double logdet = -(mu[0] + L(0, 0)*epsilon[0] + log(L(0, 0)) + log(L(1, 1)));
  return logdet;
}

double ELBO(vec x, vec mu, mat L, int n=100){
  double elbo = 0;
  for(int i = 0; i < n; ++i){
    vec epsilon = randn<vec>(2);
    double sigmaSq = exp(mu[0] + L(0, 0)*epsilon[0]);
    double phi = mu[1] + L(1, 0)*epsilon[0] + L(1, 1)*epsilon[1];
    elbo += pLogDens(x, sigmaSq, phi) - qLogDens(epsilon, mu, L);
  }
  return elbo/n;
}

vec pDeriv (vec x, double sigmaSq, double phi){
  double T = x.n_elem - 1;
  double dsigmaSq = -(T + 2.5) / sigmaSq + pow(sigmaSq, -1) + pow(x[0], 2) * (1 - pow(phi, 2)) / (2 * pow(sigmaSq, 2));
  double dphi = -phi / (1 - pow(phi, 2)) + phi * pow(x[0], 2) / sigmaSq;
  for(int t = 1; t < T+1; ++T){
    dsigmaSq -= pow(x[t] - phi*x[t-1], 2) / (2*pow(sigmaSq, 2));
    dphi += x[t-1] * (x[t] - phi*x[t-1]) / sigmaSq;
  }
  vec deriv = {dsigmaSq, dphi};
}

vec muDeriv (vec x, rowvec epsilon, vec mu, mat L){
  double sigmaSq = exp(mu[0] + L(0, 0)*epsilon[0]);
  double phi = mu[1] + L(1, 0)*epsilon[0] + L(1, 1)*epsilon[1];
  vec dpdt = pDeriv(x, sigmaSq, phi);
  vec dtdl = {sigmaSq, 1};
  vec djdl = {1, 0};
  return dpdt % dtdl + djdl;
}

mat LDeriv (vec x, rowvec epsilon, vec mu, mat L){
  double sigmaSq = exp(mu[0] + L(0, 0)*epsilon[0]);
  double phi = mu[1] + L(1, 0)*epsilon[0] + L(1, 1)*epsilon[1];
  vec dpdt = pDeriv(x, sigmaSq, phi);
  mat dtdL(2, 2);
  mat djdL(2, 2);
  
  dtdL(0, 0) = epsilon[0]*sigmaSq;
  dtdL.row(1) = epsilon;
  djdL(0, 0) = epsilon[0] + 1.0/L(0, 0);
  djdL(1, 1) = 1.0/L(1, 1);
  
  dtdL.col(0) = dtdL.col(0) % dpdt;
  dtdL.col(1) = dtdL.col(1) % dpdt;
  
  return dtdL + djdL;
}

// [[Rcpp::export]]
Rcpp::List SGA_AR1(vec x, int M, int maxIter, vec Mu, mat L, double threshold=0.01, 
                   double alpha=0.05, double beta1=0.9, double beta2=0.999){
  double e = pow(10, -8);
  
  vec MtMu(2);
  vec VtMu(2);
  vec MtMuHat(2);
  vec VtMuHat(2);
  vec pMu(2);
  vec pMuSq(2);
  
  mat MtL(2, 2);
  mat VtL(2, 2);
  mat MtLHat(2, 2);
  mat VtLHat(2, 2);
  mat pL(2, 2);
  mat pLSq(2, 2);
  
  int iter = 0;
  vec LB(maxIter+1, fill::zeros);
  LB[0] = ELBO(x, Mu, L);
  double diff = threshold + 1;
  
  while(diff > threshold | iter < 10){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    pMu.fill(0);
    pMuSq.fill(0);
    pL.fill(0);
    pLSq.fill(0);
    
    mat epsilon = randn<mat>(M, 2);
    for(int m = 0; m < M; ++m){
      vec dmu = muDeriv(x, epsilon.row(m), Mu, L);
      pMu += dmu/M;
      pMuSq += pow(dmu, 2)/M;
      mat dl = LDeriv(x, epsilon.row(m), Mu, L);
      pL += dl/M;
      pLSq += pow(dl, 2)/M;
    }
    MtMu = beta1*MtMu + (1-beta1)*pMu; // Creates biased estimates of first and second moment
    VtMu = beta2*VtMu + (1-beta2)*pMuSq;
    MtMuHat = MtMu / (1 - pow(beta1, iter)); // Corrects bias
    VtMuHat = VtMu / (1 - pow(beta2, iter));
    
    MtL = beta1*MtL + (1-beta1)*pL; // Creates biased estimates of first and second moment
    VtL = beta2*VtL + (1-beta2)*pLSq;
    MtLHat = MtL / (1 - pow(beta1, iter)); // Corrects bias
    VtLHat = VtL / (1 - pow(beta2, iter));
    
    Mu += alpha * MtMuHat / (sqrt(VtMuHat) + e);
    L += alpha * MtLHat / (sqrt(VtLHat) + e);
    LB[iter] = ELBO(x, Mu, L);
    diff = abs(LB[iter] - LB[iter-1]);
  }
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}
  
  