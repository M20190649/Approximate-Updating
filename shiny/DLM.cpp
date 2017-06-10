// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double pLogDens(vec &y, vec &x, double &sigmaSqY, double &sigmaSqX, double &phi, double &gamma){
  int TS = y.n_elem;
  double prior = - pow(gamma, 2) / 200  -  2 * log(sigmaSqY)  -  2 * log(sigmaSqX)  -  1.0/sigmaSqY  -  1.0/sigmaSqX  +
    19.0 * log(1 + phi)  +  0.5 * log(1 - phi);
  double states = 0.5*log(1-pow(phi, 2))  -  0.5*log(sigmaSqX)  -  (1-pow(phi, 2)) * pow(x[0]-gamma, 2) / (2*sigmaSqX);
  double data = 0;
  for(int t = 1; t < TS+1; ++t){
    states -= pow(x[t]-gamma-phi*(x[t-1]-gamma), 2) / (2*sigmaSqX);
    data -= pow(y[t-1]-x[t], 2) / (2*sigmaSqY);
  }
  return data  +  states  +  prior;
}

double qLogDens(vec &epsilon, vec &mu, mat &L){
  double logdet = mu[0]  +  L(0, 0) * epsilon[0]  +  mu[1]  +  L(1, 0) * epsilon[0]  +  L(1, 1) * epsilon[1];
  double eps = 0;
  for(int i = 0; i < epsilon.n_elem; ++i){
    logdet += log(abs(L(i, i)));
    eps -= 0.5 * pow(epsilon[i], 2);
  }
  return eps  -  logdet;
}

double ELBO(vec &y,  vec &mu, mat &L, int n=50){
  double elbo = 0;
  int dim = y.n_elem + 5;
  for(int i = 0; i < n; ++i){
    vec epsilon = randn<vec>(dim);
    vec transf = mu + L * epsilon;
    vec x = transf.tail(dim-4);
    double sigmaSqY = exp(transf[0]);
    double sigmaSqX = exp(transf[1]);
    elbo += pLogDens(y, x, sigmaSqY, sigmaSqX, transf[2], transf[3]) - qLogDens(epsilon, mu, L);
  }
  return elbo / n;
}

vec pDeriv (vec &y, vec &x, double &sigmaSqY, double &sigmaSqX, double &phi, double &gamma, int &T, int &S, bool &xderiv){
  vec derivs(T+S+5, fill::zeros);
  derivs[0] = -0.5 * (T+S+4) / sigmaSqY  +  pow(sigmaSqY, -2);
  derivs[1] = -0.5 * (T+S+5) / sigmaSqX  +  pow(sigmaSqX, -2)  + 
    (1-pow(phi, 2)) * pow(x[0]-gamma, 2) / (2*pow(sigmaSqX, 2));
  derivs[2] = -phi / (1-pow(phi, 2))  +   phi / sigmaSqX * pow(x[0]-gamma, 2)  +
    19.0 * phi / (1 + phi)  -  0.5 * phi / (1 - phi);
  derivs[3] = -gamma/100  +  (1-pow(phi, 2)) * (x[0]-gamma) / sigmaSqX;
  for(int t = 1; t < T+S+1; ++t){
    derivs[0] += pow(y[t-1]-x[t], 2) / (2*pow(sigmaSqY, 2));
    derivs[1] += pow(x[t]-gamma-phi*(x[t-1]-gamma), 2) / (2*pow(sigmaSqX, 2));
    derivs[2] += (x[t-1]-gamma) * (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX;
    derivs[3] += (1-phi) * (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX;
  }
  if(!xderiv){
    return derivs;
  }
  if(S == 0){
    derivs[4] = (phi*x[1]-x[0]+gamma*(1-phi)) / sigmaSqX;
    for(int i = 5; i < T+4; ++i){
      int t = i - 4;
      derivs[i] = (y[t-1]-x[t]) / sigmaSqY  -  (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX  +
        phi * (x[t+1]-gamma-phi*(x[t]-gamma)) / sigmaSqX;
    }
  } else {
    for(int i = T+5; i < T+S+4; ++i){
      int t = i - 4;
      derivs[i] = (y[t-1]-x[t]) / sigmaSqY  -  (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX +
        phi * (x[t+1]-gamma-phi*(x[t]-gamma)) / sigmaSqX;
    }
  }
  derivs[T+S+4] = (y[T+S-1]-x[T+S]) / sigmaSqY  -  (x[T+S]-gamma-phi*(x[T+S-1]-gamma)) / sigmaSqX;
  return derivs;
}

vec muDeriv (vec dpdt, double &sigmaSqY, double &sigmaSqX, int dim){
  vec dtdm(dim, fill::ones);
  dtdm[0] = sigmaSqY;
  dtdm[1] = sigmaSqX;
  vec djdm(dim, fill::zeros);
  djdm[0] = 1;
  djdm[1] = 1;
  return dpdt % dtdm + djdm;
}

mat LDeriv (vec dpdt, double &sigmaSqY, double &sigmaSqX, vec &epsilon, mat &L, int T, int S,
            bool meanfield, bool xderiv){
  mat dtdL(T+S+5, T+S+5, fill::zeros);
  mat djdL(T+S+5, T+S+5, fill::zeros);

  for(int i = 0; i < 4; ++i){
    dtdL(i, i) = epsilon[i];
    djdL(i, i) = pow(L(i, i), -1);
  }
  if(xderiv){
    if(S == 0){
      for(int i = 4; i < T+5; ++i){
        dtdL(i, i) = epsilon[i];
        djdL(i, i) = pow(L(i, i), -1);
      }
    } else {
      for(int i = T+5; i < T+S+5; ++i){
        dtdL(i, i) = epsilon[i];
        djdL(i, i) = pow(L(i, i), -1);
      }
    }
  }
  
  dtdL(0, 0) *= sigmaSqY;
  dtdL(1, 0) *= sigmaSqX;
  dtdL(1, 1) *= sigmaSqX;
  djdL(0, 0) += epsilon[0];
  djdL(1, 0) += epsilon[0];
  djdL(1, 1) += epsilon[1];

  for(int i = 1; i < 4; ++i){
    for(int j = 0; j < i; ++j){
      dtdL(i, j) = epsilon[j];
    }
  }
  
  if(!meanfield & xderiv){
    if(S == 0){
      for(int i = 4; i < T+5; ++i){
        for(int j = 0; j < i; ++j){
          dtdL(i, j) = epsilon[j];
        }
      }
    } else {
      for(int i = T+5; i < T+S+5; ++i){
        for(int j = 0; j < i; ++j){
          dtdL(i, j) = epsilon[j];
        }
      }
    }
  }

  dtdL.each_col() %= dpdt;
  
  return dtdL + djdL;
}

// [[Rcpp::export]]
Rcpp::List SGA_DLM(vec y, int M, int maxIter, vec Mu, mat L, int S=0, bool meanfield=false, bool xderiv=true,
  bool variance=true, double threshold=0.01, double alpha=0.1, double beta1=0.9, double beta2=0.999){
  double e = pow(10, -8);
  int T = y.n_elem - S;
  
  vec MtMu(T+S+5, fill::zeros);
  vec VtMu(T+S+5, fill::zeros);
  vec MtMuHat(T+S+5);
  vec VtMuHat(T+S+5);
  vec pMu(T+S+5);
  vec pMuSq(T+S+5);
  
  mat MtL(T+S+5, T+S+5, fill::zeros);
  mat VtL(T+S+5, T+S+5, fill::zeros);
  mat MtLHat(T+S+5, T+S+5);
  mat VtLHat(T+S+5, T+S+5);
  mat pL(T+S+5, T+S+5);
  mat pLSq(T+S+5, T+S+5);
  
  int iter = 0;
  vec LB(maxIter+1, fill::zeros);
  LB[0] = ELBO(y, Mu, L);
  double diff = threshold + 1;
  double meanLB = 0;
  double omeanLB;
  
  while(diff > threshold | iter < 100){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    pMu.fill(0);
    pMuSq.fill(0);
    pL.fill(0);
    pLSq.fill(0);
    
    mat epsilon = randn<mat>(T+S+5, M);
    for(int m = 0; m < M; ++m){
      vec transf = Mu + L * epsilon.col(m);
      vec x = transf.tail(T+S+1);
      double sigmaSqY = exp(transf[0]);
      double sigmaSqX = exp(transf[1]);
      vec dpdt = pDeriv(y, x, sigmaSqY, sigmaSqX, transf[2], transf[3], T, S, xderiv);
      vec dmu = muDeriv(dpdt, sigmaSqY, sigmaSqX, T+S+5);
      pMu += dmu / M;
      pMuSq += pow(dmu, 2) / M;
      if(variance){
        vec eps = epsilon.col(m);
        mat dl = LDeriv(dpdt, sigmaSqY, sigmaSqX, eps, L, T, S, meanfield, xderiv);
        pL += dl / M;
        pLSq += pow(dl, 2) / M;
      }
    }
    MtMu = beta1*MtMu  +  (1-beta1)*pMu; // Creates biased estimates of first and second moment
    VtMu = beta2*VtMu  +  (1-beta2)*pMuSq;
    MtMuHat = MtMu / (1-pow(beta1, iter)); // Corrects bias
    VtMuHat = VtMu / (1-pow(beta2, iter));
  
    MtL = beta1*MtL  +  (1-beta1)*pL; // Creates biased estimates of first and second moment
    VtL = beta2*VtL  +  (1-beta2)*pLSq;
    MtLHat = MtL / (1-pow(beta1, iter)); // Corrects bias
    VtLHat = VtL / (1-pow(beta2, iter));
    if(iter > 0){
      Mu += alpha * MtMuHat / (sqrt(VtMuHat)+e);
      L += alpha * MtLHat / (sqrt(VtLHat)+e);
    }
    LB[iter] = ELBO(y, Mu, L);
    if(iter % 5 == 0){
      omeanLB = meanLB;
      meanLB = 0.2 * (LB[iter]+LB[iter-1]+LB[iter-2]+LB[iter-3]+LB[iter-4]);
      diff = abs(meanLB - omeanLB);
    }
    if(iter % 100 == 0 & diff > threshold){
      Rcpp::Rcout << iter << std::endl;
    }
  } // close while loop
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  Rcpp::Rcout << iter << std::endl;
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}

