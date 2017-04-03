// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
//#include <distributions.h> //figure this out
using namespace Rcpp;
using namespace arma;
using namespace std;

class Distribution{
public:
  virtual void set_values(int, int) {};
  virtual double logdens(double x) {return -99.99;};
  virtual vec sample (int n_samples) {return zeros<vec>(n_samples);};
};

class Normal: public Distribution{
  // This just provides methods for the single variate normal distribution
public:
  double mean, variance;
  Normal(double, double);
  void set_values(double m, double v){
    mean = m;
    variance = v;
  }
  void increase_values(double m, double sd){
    mean += m;
    if(sqrt(variance) > -sd){
      variance = pow(sqrt(variance) + sd, 2);
    }
  }
  double logdens (double x) {
    return -0.5 * log(2*3.141593*variance) - pow(x - mean, 2) / (2 * variance);}
  double transform_epsilon(double x){
    return mean + sqrt(variance) * x;
  }
};

Normal::Normal (double m, double v){
  mean = m;
  variance = v;
}

class InverseGamma: public Distribution{
  // The distribution for 1/x  where x ~ gamma(shape, rate)
  double shape, scale;
public:
  InverseGamma(double, double);
  void set_values(double sh, double sc){
    shape = sh;
    scale = sc;
  }
  double logdens (double x) {
    // the input x should be log(sigmaSq)
    double constant = shape * log(scale) - log(tgamma(shape));
    double kernel = -(shape + 1) * x - scale / exp(x);
    return constant + kernel;
  }
  vec sample (int n_samples){
    // The sampler returns sigmaSq instead of log(sigmaSq)
    // We only need the sampler if q ~ IG, and in this case we don't need the log transform
    // I don't need this but it's here
    return 1.0 / randg<vec>(n_samples, distr_param(shape, 1/scale));
  }
};

InverseGamma::InverseGamma (double sh, double sc){
  shape = sh;
  scale = sc;
}

cube qSim (Normal* qLatent[], int n_samples, int p){
  cube output(n_samples, p, 2);
  mat epsilon = randn<mat>(n_samples, p);   
  mat theta (n_samples, p);
  for(int i = 0; i < n_samples; ++i){
    for(int j = 0; j < p; ++j){
      theta(i, j) = qLatent[j]->transform_epsilon(epsilon(i, j));
    }
  }
  output.slice(0) = theta;
  output.slice(1) = epsilon;
  return output;
}

double pdens (Normal* Y[], vec y, Distribution* pLatent[], rowvec latent, int T, int p){
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
  
double qdens (Normal* qLatent[], rowvec latent, int p){
  double approxDens = 0;
  for(int i = 0; i < p; ++i){
    approxDens += qLatent[i]->logdens(latent[i]);
  }
  return approxDens;
}
  
void updateP (Normal* Y[], Distribution* pLatent[], rowvec latent, int T){
  dynamic_cast<Normal*>(pLatent[4])->set_values(0.0, 4*exp(latent[1]));
  for(int t = 0; t < T; ++t){
    dynamic_cast<Normal*>(pLatent[t+5])->set_values(latent[2] * latent[t+4], exp(latent[1]));
    Y[t]->set_values(latent[3] + latent[t+5], exp(latent[0]));
  }
} 

void updateQ (Normal* qLatent[], mat increase, int p){
  for(int i = 0; i < p; ++i){
    qLatent[i]->increase_values(increase(i, 0), increase(i, 1));
  }
}

double ELBO(Normal* Y[], vec y, Distribution* pLatent[], Normal* qLatent[], int T, int p, int n = 100){
  mat latentSims = qSim(qLatent, n, p).slice(0);
  double value = 0;
  for(int i = 0; i < n; ++i){
    updateP (Y, pLatent, latentSims.row(i), T);
    value += pdens(Y, y, pLatent, latentSims.row(i), T, p) - qdens(qLatent, latentSims.row(i), p);
  }
  return value / n;
}

double logjointDeriv (vec y, rowvec latent, int T, int i){
  double dpdf = 0;
  latent[0] = exp(latent[0]); // transform to sigmaSq to make p(y, theta) a little bit easier to deal with
  latent[1] = exp(latent[1]);
  // Next is the derivatives of p(y, theta) wrt theta.
  if(i == 0){ //sigmaSqY
    dpdf = -(T/2 + 2) / latent[i] + 1 / pow(latent[i], 2);
    for(int t = 0; t < T; ++t){
      dpdf += pow(y[t] - latent[3] - latent[t+5], 2) / (2 * pow(latent[i], 2));
    }
  } else if(i == 1){  //sigmaSqX
    dpdf = -(T/2 + 5/2) / latent[i] + 1 / pow(latent[i], 2) + pow(latent[4], 2) / (2 * pow(latent[i], 2));
    for(int t = 5; t < T+5; ++t){
      dpdf += pow(latent[t] - latent[2]*latent[t-1], 2) / (2 * pow(latent[i], 2));
    }
  } else if(i == 2){ //Phi
    for(int t = 5; t < T+5; ++t){
      dpdf += (latent[t]*latent[t-1] - latent[i]*pow(latent[t], 2)) / latent[1];
    }
  } else if(i == 3){ //Mu
    dpdf = - latent[i] / 100;
    for(int t = 0; t < T; ++t){
      dpdf -= (latent[i] + latent[t+5] - y[t]) / latent[0];
    }
  } else if(i == 4){ //X0
    dpdf = (latent[2]*(latent[5]) - (latent[2] + 1)*latent[i]) / latent[1];
  } else if(i == T+5){ //XT
    dpdf = (y[T] - latent[3] - latent[i]) / latent[0] - (latent[i] - latent[2]*latent[i-1]) / latent[1];
  } else { //X1 ... XT-1
    dpdf = (y[i-5] - latent[3] - latent[i]) / latent[0] - (latent[i] - latent[2]*latent[i-1]) / latent[1] +
      latent[2] * (latent[i+1] - latent[2] * latent[i]) / latent[1];
  }
  return dpdf;
}

rowvec elboDeriv (vec y, Normal* qLatent[], rowvec latent, rowvec epsilon, int T, int i){
  
  double dpdf = logjointDeriv(y, latent, T, i); //Initialise to 0 
  double dfdm = 1; //theta = mu + sum(L*eps) for i > 2. Derivative of 1.
  double dfds = epsilon[i];
  // lnq = -( mu1 + mu2 + sum(log(Lii) for i = 1, ..., p) + L11eps1 + L21eps1 + L22eps2) + ln(p(eps))
  double dqdm = 0; // dq/dmu is zero except for i = 1, 2
  double dqds = -1 / sqrt(qLatent[i]->variance);
 
  if(i < 2){
    dqdm = -1; // dq/dmu is 1 for i = 1, 2
    dqds -= epsilon[i]; // dq/dL11 = 1/L11 + eps1, dq/dL21 = eps1, so adding eps1 to above
    // theta1 = exp(mu1 + L11 eps1), theta2 = exp(mu2 + L21 eps1 + L22 eps2)
    dfdm = latent[i]; // derivative wrt mu is theta
    dfds *= latent[i]; // derivative wrt Lij is theta_i * epsilon_j for j <= i, dfdL already contains relevant epsilon terms
  }
 
  rowvec Derivs  = {dfdm * dpdf - dqdm, dfds * dpdf - dqds};
  return Derivs;
}

// [[Rcpp::export]]
mat DLM_SGA (vec y, int S, int maxIter, double threshold, mat values, bool initial){
  int T = y.n_elem;
  mat partials(T+5, 2);
  mat Gt(T+5, 2, fill::zeros);
  mat Pt(T+5, 2);
  double eta = 0.1;
  double iter = 0;
  double diff = threshold + 1;
  Normal *Y[T];
  Distribution *pLatent[T+5];
  Normal *qLatent[T+5];
  if(initial){
    pLatent[0] = new InverseGamma(1, 1); qLatent[0] = new Normal(values(0, 0), values(0, 1));
    pLatent[1] = new InverseGamma(1, 1); qLatent[1] = new Normal(values(1, 0), values(1, 1));
    pLatent[2] = new Normal(0, 0.25); qLatent[2] = new Normal(values(2, 0), values(2, 1));
    pLatent[3] = new Normal(0, 5); qLatent[3] = new Normal(values(3, 0), values(3, 1));
    pLatent[4] = new Normal(0, 5); qLatent[4] = new Normal(values(4, 0), values(4, 1));
    for(int t = 0; t < T; ++t){
      pLatent[t+5] = new Normal(0, 1);
      Y[t] = new Normal(0, 1);
      qLatent[t+5] = new Normal(values(t+5, 0), values(t+5, 1));
    }
  } else {
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
  }
  double LBold;
  double LBnew = ELBO(Y, y, pLatent, qLatent, T, T+5);
  
  while(diff > threshold | iter < 10){
    iter = iter + 1;
    if(iter > maxIter){
      break;
    }
    partials.fill(0);
    cube Sims = qSim(qLatent, S, T+5);
    mat latentSims = Sims.slice(0); 
    mat epsilon = Sims.slice(1);
    for(int s = 0; s < S; ++s){
      updateP(Y, pLatent, latentSims.row(s), T);
      for(int i = 0; i < T+5; ++i){
        partials.row(i) += elboDeriv(y, qLatent, latentSims.row(s), epsilon.row(s), T, i)/S;
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
  

