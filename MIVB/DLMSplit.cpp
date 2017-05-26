// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double pLogDens(vec y, vec x, double sigmaSqY, double sigmaSqX, double phi, double gamma){
  int T = y.n_elem;
  double prior = - pow(gamma, 2) / 200 - 2 * log(sigmaSqY) - 2 * log(sigmaSqX) - 1.0/sigmaSqY - 1.0/sigmaSqX;
  double states = 0.5 * log(1-pow(phi, 2)) - (1 - pow(phi, 2)) * pow(x[0] - gamma, 2) / (2 * sigmaSqX);
  double data = 0;
  for(int t = 1; t < T+1; ++t){
    states -= pow(x[t] - gamma - phi * (x[t-1] - gamma), 2) / (2 * sigmaSqX);
    data -= pow(y[t-1] - x[t], 2) / (2 * sigmaSqY);
  }
  return data + states + prior;
}

double qLogDens(vec epsilon, vec mu, mat L){
  double logdet = mu[0] + L(0, 0)*epsilon[0] + mu[1] + L(1, 0)*epsilon[0] + L(1, 1)*epsilon[1];
  double eps = 0;
  for(int i = 0; i < epsilon.n_elem; ++i){
    logdet += log(abs(L(i, i)));
    eps -= 0.5 * pow(epsilon[i], 2);
  }
  return eps - logdet;
}

double ELBO(vec y,  vec mu, mat L, int n=25){
  double elbo = 0;
  int T = y.n_elem;
  for(int i = 0; i < n; ++i){
    vec epsilon = randn<vec>(T+5);
    vec transf = mu + L * epsilon;
    elbo += pLogDens(y, transf.tail(T+1), exp(transf[0]), exp(transf[1]), transf[2], transf[3]) - qLogDens(epsilon, mu, L);
  }
  return elbo/n;
}

vec pDeriv (vec y, vec x, double sigmaSqY, double sigmaSqX, double phi, double gamma, bool xderiv){
  double T = y.n_elem;
  vec derivs(T+5, fill::zeros);
  derivs[0] = -0.5*(T + 4) / sigmaSqY + pow(sigmaSqY, -2);
  derivs[1] = -0.5*(T + 5) / sigmaSqX + pow(sigmaSqX, -2) + 
    (1-pow(phi, 2)) * pow(x[0] - gamma, 2) / (2 * pow(sigmaSqX, 2));
  derivs[2] = -phi / (1 - pow(phi, 2)) +  phi / sigmaSqX * pow(x[0] - gamma, 2);
  derivs[3] = -gamma/100 + (1-pow(phi, 2)) * (x[0] - gamma) / sigmaSqX;
  for(int t = 1; t < T+1; ++t){
    derivs[0] += pow(y[t-1] - x[t], 2) / (2 * pow(sigmaSqY, 2));
    derivs[1] += pow(x[t] - gamma - phi*(x[t-1]-gamma), 2) / (2*pow(sigmaSqX, 2));
    derivs[2] += (x[t-1] - gamma) * (x[t] - gamma - phi*(x[t-1]-gamma)) / sigmaSqX;
    derivs[3] += (1 - phi) * (x[t] - gamma - phi*(x[t-1]-gamma)) / sigmaSqX;
  }
  if(!xderiv){
    return derivs;
  }
  derivs[4] = (phi*x[1] - x[0] + gamma*(1-phi)) / sigmaSqX;
  for(int i = 5; i < T+4; ++i){
    int t = i - 4;
    derivs[i] = (y[t-1] - x[t]) / sigmaSqY - (x[t] - gamma - phi*(x[t-1]-gamma)) / sigmaSqX +
      phi * (x[t+1] - gamma - phi*(x[t]-gamma)) / sigmaSqX;
  }
  derivs[T+4] = (y[T-1] - x[T]) / sigmaSqY - (x[T] - gamma - phi*(x[T-1]-gamma)) / sigmaSqX;
  return derivs;
}

vec muDeriv (vec dpdt, double sigmaSqY, double sigmaSqX, int T){
  vec dtdm(T+5, fill::ones);
  dtdm[0] = sigmaSqY;
  dtdm[1] = sigmaSqX;
  vec djdm(T+5, fill::zeros);
  djdm[0] = 1;
  djdm[1] = 1;
  return dpdt % dtdm + djdm;
}

mat LDeriv (vec dpdt, double sigmaSqY, double sigmaSqX, vec epsilon, mat L, int T, bool meanfield, bool xderiv){
  mat dtdL(T+5, T+5, fill::zeros);
  mat djdL(T+5, T+5, fill::zeros);
  mat dpdtM(T+5, T+5);
  dpdtM.each_col() = dpdt;
  
  for(int i = 0; i < T+5; ++i){
    if((i < 4) | xderiv){
      dtdL(i, i) = epsilon[i];
      djdL(i, i) = 1.0/L(i, i);
    }
  }
  dtdL(0, 0) *= sigmaSqY;
  dtdL(1, 1) *= sigmaSqX;
  djdL(0, 0) += epsilon[0];
  djdL(1, 1) += epsilon[1];
 
  if(!meanfield){
    for(int i = 1; i < T+5; ++i){
      for(int j = 0; j < i; ++j){
        dtdL(i, j) = epsilon[j];
      }
    }
    dtdL(1, 0) *= sigmaSqX;
    djdL(1, 0) += epsilon[0];
  }
  
  return dpdtM % dtdL + djdL;
}

// [[Rcpp::export]]
Rcpp::List SGA_DLM(vec y, int M, int maxIter, vec Mu, mat L, bool meanfield=false, bool xderiv=true,
	double threshold=0.01, double alpha=0.05, double beta1=0.9, double beta2=0.999){
  double e = pow(10, -8);
  int T = y.n_elem;
  
  vec MtMu(T+5, fill::zeros);
  vec VtMu(T+5, fill::zeros);
  vec MtMuHat(T+5);
  vec VtMuHat(T+5);
  vec pMu(T+5);
  vec pMuSq(T+5);
  
  mat MtL(T+5, T+5, fill::zeros);
  mat VtL(T+5, T+5, fill::zeros);
  mat MtLHat(T+5, T+5);
  mat VtLHat(T+5, T+5);
  mat pL(T+5, T+5);
  mat pLSq(T+5, T+5);
  
  int iter = 0;
  vec LB(maxIter+1, fill::zeros);
  LB[0] = ELBO(y, Mu, L);
  double diff = threshold + 1;
  
  while(diff > threshold | iter < 100){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    pMu.fill(0);
    pMuSq.fill(0);
    pL.fill(0);
    pLSq.fill(0);
    
    mat epsilon = randn<mat>(T+5, M);
    for(int m = 0; m < M; ++m){
      vec transf = Mu + L * epsilon.col(m);
      vec dpdt = pDeriv(y, transf.tail(T+1), exp(transf[0]), exp(transf[1]), transf[2], transf[3], xderiv);
      vec dmu = muDeriv(dpdt, exp(transf[0]), exp(transf[1]), T);
      pMu += dmu/M;
      pMuSq += pow(dmu, 2)/M;
      //mat dl = LDeriv(dpdt, exp(transf[0]), exp(transf[1]), epsilon.col(m), L, T, meanfield, xderiv);
      //pL += dl/M;
      //pLSq += pow(dl, 2)/M;
    }
    
    MtMu = beta1*MtMu + (1-beta1)*pMu; // Creates biased estimates of first and second moment
    VtMu = beta2*VtMu + (1-beta2)*pMuSq;
    MtMuHat = MtMu / (1 - pow(beta1, iter)); // Corrects bias
    VtMuHat = VtMu / (1 - pow(beta2, iter));
    
    MtL = beta1*MtL + (1-beta1)*pL; // Creates biased estimates of first and second moment
    VtL = beta2*VtL + (1-beta2)*pLSq;
    MtLHat = MtL / (1 - pow(beta1, iter)); // Corrects bias
    VtLHat = VtL / (1 - pow(beta2, iter));
    if(iter > 0){
      Mu += alpha * MtMuHat / (sqrt(VtMuHat) + e);
      L += alpha * MtLHat / (sqrt(VtLHat) + e);
    }
    LB[iter] = ELBO(y, Mu, L);
    diff = abs(LB[iter] - LB[iter-1]);
  }
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  Rcpp::Rcout << iter << std::endl;
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}

