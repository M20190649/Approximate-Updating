// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
//#include <distributions.h> //figure this out
using namespace Rcpp;
using namespace arma;
using namespace std;

#define Pi 3.14159265358979

class Distribution{
public:
  virtual void set_values(int, int) {};
  virtual double logdens(double x) {return -99.99;};
  virtual vec sample (int n_samples) {return zeros<vec>(n_samples);};
};

class Normal: public Distribution{
public:
  double mean, variance;
  Normal(double, double);
  void set_values(double m, double v){
    mean = m;
    variance = v;
  }
  void increase_values(double m, double st){
    mean += m;
    if(sqrt(variance) > -st){
      variance = pow(sqrt(variance) + st, 2);
    }
  }
  double logdens (double x) {
    return -0.5 * log(2*3.141593*variance) - pow(x - mean, 2) / (2 * variance);}
  vec sample (int n_samples){
    return  mean + sqrt(variance) * randn<vec>(n_samples);
  }
  double MeanDeriv(double x){
    return (x - mean) / variance;
  }
  double VarDeriv(double x){
    return (pow(x - mean, 2)/(2*pow(variance, 2)) - 1 / (2*variance));
  }
  double StdDeriv(double x){
    return pow(x - mean, 2)/(pow(variance, 1.5)) - 1 / sqrt(variance); 
  }
  double XDeriv(double x){
    return (mean - x) / variance;
  }
};

Normal::Normal(double m, double v){
  mean = m;
  variance = v;
}

class InverseGamma: public Distribution{
public:
  double shape, scale;
  InverseGamma(double, double);
  void set_values(double sh, double sc){
    shape = sh;
    scale = sc;
  }
  double logdens (double x) {
    double cons = shape * log(scale) - log(tgamma(shape));
    double kernel = -(shape + 1) * log(x) - scale / x;
    return cons + kernel;
  }
  vec sample (int n_samples){
    return 1.0 / randg<vec>(n_samples, distr_param(shape, 1/scale));
  }
};

InverseGamma::InverseGamma(double sh, double sc){
  shape = sh;
  scale = sc;
}

mat qSim (Distribution* qLatent[], int n_samples, int p){
  mat output(n_samples, p);
  for(int i = 0; i < p; ++i){
    output.col(i) = qLatent[i]->sample(n_samples);
    if(i < 2){
      output.col(i) = exp(output.col(i));
    }
  }
  return output;
}

double pdens (Distribution* Y[], vec y, Distribution* pLatent[], rowvec latent, int T, int p){
  double priorDens = 0;
  for(int i = 0; i < p; ++i){
    priorDens += pLatent[i]->logdens(latent[i]);
  }
  double dataDens = 0;
  for(int t= 0; t < T; ++t){
    dataDens += Y[t]->logdens(y[t]);
  }
  return priorDens + dataDens;
}
  
double qdens (Distribution* qLatent[], rowvec latent, int p){
  double approxDens = 0;
  for(int i = 0; i < p; ++i){
    approxDens += qLatent[i]->logdens(latent[i]);
  }
  return approxDens;
}
  
void updateP (Distribution* Y[], Distribution* pLatent[], rowvec latent, int T){
  for(int t = 0; t < T; ++t){
    dynamic_cast<Normal*>(pLatent[t+5])->set_values(latent[2] * latent[t+4], latent[1]);
    dynamic_cast<Normal*>(Y[t])->set_values(latent[3] + latent[t+5], latent[0]);
  }
} 

void updateQ (Distribution* qLatent[], mat increase, int p){
  for(int i = 0; i < p; ++i){
    dynamic_cast<Normal*>(qLatent[i])->increase_values(increase(i, 0), increase(i, 1));
  }
}
  
double ELBO(Distribution* Y[], vec y, Distribution* pLatent[], Distribution* qLatent[], int T, int p, int n = 100){
  mat latentSims = qSim(qLatent, n, p);
  double value = 0;
  for(int i = 0; i < n; ++i){
    updateP (Y, pLatent, latentSims.row(i), T);
    value += pdens(Y, y, pLatent, latentSims.row(i), T, p) - qdens(qLatent, latentSims.row(i), p);
  }
  return value / n;
}

double logjointDeriv (vec y, rowvec latent, int T, int i){
  double dpdf = 0;
  if(i == 0){ //sigmaSqY
    dpdf = -(T/2 + 2)/latent[i] + 1/pow(latent[i], 2);
    for(int t = 0; t < T; ++t){
      dpdf += pow(y[t] - latent[3] - latent[t+5], 2) / (2 * pow(latent[i], 2));
    } 
  } else if(i == 1) { //corresponds to sigmaSqX
    dpdf = -(T/2 + 5/2)/latent[i] + 1/pow(latent[i], 2) + (1-pow(latent[2],2))*pow(latent[4], 2) / (2 * pow(latent[i], 2));
    for(int t = 0; t < T; ++t){
      dpdf += pow(latent[t+5] - latent[2]*latent[t+4], 2) / (2 * pow(latent[i], 2));
    }
  } else if(i == 2) { //phi
    for(int t = 0; t < T; ++t){
      dpdf -= latent[i] * pow(latent[t+4], 2) + latent[t+4] * latent[t+5] / latent[1];
    }
  } else if(i == 3) { //mu
    dpdf = - latent[i] / 10;
    for(int t = 0; t < T; ++t){
      dpdf += - (latent[i] - y[t] + latent[t+5]) / latent[0];
    }
  } else if(i == 4) { //x0
    dpdf = - (latent[i]*(1 + pow(latent[2], 2)) - latent[2]*latent[5]) / latent[1];
  } else if(i == T + 4) { //xT
    dpdf = - (latent[T+4] - latent[2]*latent[T+3]) / latent[1] -
      (latent[T+4] + latent[3] - y[T]) / latent[0];
  } else { //xt, t = 1, 2, ... ,T-1
    dpdf = - ((1 + pow(latent[1],2))*latent[i] - latent[1]*(latent[i-1] + latent[i+1])) / latent[1] -
      (latent[i] + latent[3] - y[i-5]) / latent[0];
  }
  return dpdf;
}

rowvec elboDeriv (vec y, Distribution* qLatent[], rowvec latent, int T, int i){
  double dpdf = logjointDeriv(y, latent, T, i);
  double dfdlMu = 1;
  if(i < 2){
    dfdlMu = latent[i];
  }
  double dfdlSigma = sqrt(dynamic_cast<Normal*>(qLatent[i])->variance);
  double dqdSigma = 1/sqrt(dynamic_cast<Normal*>(qLatent[i])->variance);
  rowvec Derivs  = {dfdlMu * dpdf, dfdlSigma * dpdf - dqdSigma};
  return Derivs;
}

// [[Rcpp::export]]
mat DLM_SGA (vec y, int S, int maxIter){
  int T = y.n_elem;
  mat partials(T+5, 2);
  mat Gt(T+5, 2, fill::zeros);
  mat Pt(T+5, 2);
  double eta = 0.1;
  double iter = 0;
  double threshold = 0.01;
  double diff = threshold + 1;
  Distribution *Y[T];
  Distribution *pLatent[T+5];
  Distribution *qLatent[T+5];
  pLatent[0] = new InverseGamma(1, 1); qLatent[0] = new Normal(0, 1);
  pLatent[1] = new InverseGamma(1, 1); qLatent[1] = new Normal(0, 1);
  pLatent[2] = new Normal(0, 0.25); qLatent[2] = new Normal(0, 1);
  pLatent[3] = new Normal(0, 5); qLatent[3] = new Normal(0, 1);
  pLatent[4] = new Normal(0, 5); qLatent[4] = new Normal(0, 1);
  for(int t = 0; t < T; ++t){
    pLatent[t+5] = new Normal(0, 1);
    Y[t] = new Normal(0, 1);
    qLatent[t+5] = new Normal(0, 1);
  }
  double LBold;
  double LBnew = ELBO(Y, y, pLatent, qLatent, T, T+5);
  
  while(diff > threshold | iter < 10){
    iter = iter + 1;
    if(iter > maxIter){
      break;
    }
    partials.fill(0);
    mat latentSim = qSim(qLatent, S, T+5);
    for(int s = 0; s < S; ++s){
      updateP(Y, pLatent, latentSim.row(s), T);
      for(int i = 0; i < T+5; ++i){
        partials.row(i) += elboDeriv(y, qLatent, latentSim.row(s), T, i)/S;
      }
    }
    Gt += pow(partials, 2);
    Pt = eta * pow(Gt, -0.5);
    updateQ(qLatent, Pt % partials, T+5);
    LBold = LBnew;
    LBnew = ELBO(Y, y, pLatent, qLatent, T, T+5);
    diff = abs(LBnew - LBold);
  }
  Rcpp::Rcout << "Number of Iterations: " << iter << std::endl;
  Rcpp::Rcout << "Final ELBO: " << LBnew << std::endl;
  Rcpp::Rcout << "Final Difference in ELBO: " << diff << std::endl;
  mat lambda(T+5, 2);
  for(int i = 0; i < T+5; ++i){
    lambda(i, 0) = dynamic_cast<Normal*>(qLatent[i])->mean;
    lambda(i, 1) = dynamic_cast<Normal*>(qLatent[i])->variance;
  }
  return lambda;
}
  

