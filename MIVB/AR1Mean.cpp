// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double pLogDens(vec x, double sigmaSq, double phi, double gamma){
  int T = x.n_elem - 1;
  double density = -0.5*(T+5)*log(sigmaSq) + 0.5*log(1-pow(phi, 2)) - pow(sigmaSq, -1) - pow(gamma, 2)/200;
  density -= pow(x[0] - gamma/(1-phi), 2) * (1 - pow(phi, 2)) / (2 * sigmaSq);
  for(int t = 1; t < T+1; ++t){
    density -= pow(x[t] - gamma - phi*(x[t-1]-gamma), 2) / (2*sigmaSq);
  }
  return density;
}

double qLogDens(vec epsilon, vec mu, mat L){
  double logdet = -(mu[0] + L(0, 0)*epsilon[0] + log(L(0, 0)) + log(L(1, 1)) + log(L(2, 2)));
  double eps = 1.5*log(2.0*3.14159) - 0.5 * (pow(epsilon[0], 2) + pow(epsilon[1], 2) + pow(epsilon[2], 2));
  return logdet + eps;
}

double ELBO(vec x, vec mu, mat L, int n=100){
  double elbo = 0;
  for(int i = 0; i < n; ++i){
    vec epsilon = randn<vec>(3);
    vec trans = mu + L * epsilon;
    elbo += pLogDens(x, exp(trans[0]), trans[1], trans[2]) - qLogDens(epsilon, mu, L);
  }
  return elbo/n;
}

vec pDeriv (vec x, double sigmaSq, double phi, double gamma){
  double T = x.n_elem - 1;
  double dsigmaSq = -0.5*(T + 5) / sigmaSq + pow(sigmaSq, -1) + pow(x[0]-gamma/(1-phi), 2) * (1 - pow(phi, 2)) / (2 * pow(sigmaSq, 2));
  double dphi = -phi / (1 - pow(phi, 2)) - 1.0 / sigmaSq * ((gamma + x[0]*(phi-1))*(gamma-x[0]*(phi-1)*phi) / (pow(1-phi, 2)));
  double dgamma = -gamma/100 + 1.0 / sigmaSq * ((phi+1)*(gamma+x[0]*(phi-1))/(phi-1));
  for(int t = 1; t < T+1; ++T){
    dsigmaSq -= pow(x[t] - gamma - phi*(x[t-1]-gamma), 2) / (2*pow(sigmaSq, 2));
    dphi += (x[t-1] - gamma) * (x[t] - gamma - phi*(x[t-1]-gamma)) / sigmaSq;
    dgamma += (1 - phi) / sigmaSq * (x[t] - gamma - phi*(x[t-1]-gamma));
  }
  vec deriv = {dsigmaSq, dphi, dgamma};
  return deriv;
}

vec muDeriv (vec dpdt, double sigmaSq){
  vec dtdm = {sigmaSq, 1, 1};
  vec djdm = {1, 0, 0};
  return dpdt % dtdm + djdm;
}

mat LDeriv (vec x, vec epsilon, vec mu, mat L, vec dpdt){
  mat dtdL(3, 3);
  mat djdL(3, 3);
  
  dtdL(0, 0) = epsilon[0]*exp(mu[0]+L(0,0)*epsilon[0]);
  dtdL.row(1) = epsilon.head(2).t();
  dtdL.row(2) = epsilon.t();
  djdL(0, 0) = epsilon[0] + 1.0/L(0, 0);
  djdL(1, 1) = 1.0/L(1, 1);
  djdL(2, 2) = 1.0/L(1, 1);
  
  dtdL.col(0) = dtdL.col(0) % dpdt;
  dtdL.col(1) = dtdL.col(1) % dpdt;
  dtdL.col(2) = dtdl.col(2) % dpdt;
  
  return dtdL + djdL;
}

// [[Rcpp::export]]
Rcpp::List SGA_AR1(vec x, int M, int maxIter, vec Mu, mat L, double threshold=0.01, 
                   double alpha=0.05, double beta1=0.9, double beta2=0.999){
  double e = pow(10, -8);
  
  vec MtMu(3);
  vec VtMu(3);
  vec MtMuHat(3);
  vec VtMuHat(3);
  vec pMu(3);
  vec pMuSq(3);
  
  mat MtL(3, 3);
  mat VtL(3, 3);
  mat MtLHat(3, 3);
  mat VtLHat(3, 3);
  mat pL(3, 3);
  mat pLSq(3, 3);
  
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
    
    mat epsilon = randn<mat>(2, M);
    for(int m = 0; m < M; ++m){
      vec trans = Mu + L * epsilon.col(m);
      vec dpdt = pDeriv(x, exp(trans[0]), trans[1], trans[2]);
      vec dmu = muDeriv(dpdt, exp(trans[0]));
      pMu += dmu/M;
      pMuSq += pow(dmu, 2)/M;
      mat dl = LDeriv(x, epsilon.col(m), Mu, L, dpdt);
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

