// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double alphaTrueDens(double a, int t){
  
}

double normDens(double x, double mu, double sigmaSq){
  return 1 / sqrt(2*3.14159*sigmaSq) * exp(-(x - mu)*(x - mu) / (2*sigmaSq));
}

double alphaCanSample(){
  double u = randu<vec>(1)[0];
  double draw = randn<vec>(1)[0]
  if(u < 1){
    draw = m + sd*draw;
  } else if(u < 2){
    
  } else if(u < 3){
    
  } else if(u < 4){
    
  } else if(u < 5){
    
  } else if(u < 6){
    
  } else {
    
  }
  return draw;
}

double alphaCanDens(double x){
  double density = w1 * normDens(x, m1, var1) +
    w2 * normDens(x, m2, var2) + 
    w3 * normDens(x, m3, var3);
  return density;
}

double phiTrueDens(double phi, double sigmaSq, double gamma, rowvec alpha){
  density = 0;
  if(phi > -1 & phi < 1){
    
  }
  return density;
}


// [[Rcpp::export]]
Rcpp::List DLM_MCMC(vec y, int reps){
  int T = y.n_elem;
  mat theta(reps, 3);
  mat a(reps, T+1), fill::randn);
  theta.row(0) = {1, 0, 0};
  vec accept(T+2, fill::zeros);

  double sigmaSqShape = 0.5*T + 1.5;
  double sigmaSqScale;
  double meanGammaNumer;
  double meanGammaDenom;
  
  for(int i = 1; i < reps; ++i){
    //sigmaSq ~ IG(shape, scale)
    sigmaSqScale = 1;
    for(int t = 0; t < T; ++t){
      sigmaSqScale += pow(y[t] - x(i-1, t+1) - theta(i-1, 3), 2) / 2;
    }
    theta(i, 0) = (1 / randg<vec>(1, distr_param(sigmaSqShape, 1/sigmaSqScale)))[0];
    
    //gamma ~ Normal(mean, var)
    meanGammaNumer = 0;
    meanGammaDenom = 10*T + theta(i, 0);
    for(int t = 0; t < T; ++t){
      meanGammaNumer += 10 * y[t] - x(i-1, t+1);
    }
    theta(i, 1) = (meanGammaNumer/meanGammaDenom + sqrt(10*theta(i, 0)/meanGammaDenom) * randn<vec>(1))[0];
    
    //phi ~ TruncNormal(mean, var, low = -1, high = 1) Metropolis Hastings
    double candidate = theta(i-1, 2) + 0.1*randn<vec>(1);
    double oldp = phiTrueDens(theta(i-1, 2), theta(i, 0), theta(i, 1), a.row(i-1));
    double newp = phiTrueDens(candidate, theta(i, 0), theta(i, 1), a.row(i-1));
    double u = randu<vec>(1)[0];
    if(u < ratio){
      theta(i, 2) = candidate;
      accept(0) += 1;
    } else {
      theta(i, 2) = theta(i-1, 2);
    }
    
    for(int t = 0; t < T+1; ++t){
      double candidate = alphaCanSample();
      double oldp = alphaTrueDens(a(i-1, t));
      double oldq = alphaCanDens(a(i-1, t));
      double newp = alphaTrueDens(candidate);
      double newq = alphaCanDens(candidate);
      double ratio = newp + oldq - newq - oldp;
      double u = randu<vec>(1)[0];
      if(u < ratio){
        a(i, t) = candidate;
        accept(i+1) += 1;
      } else {
        a(i, t) = a(i-1, t);
      }
    }
  
    
  }
  
  return Rcpp::List::create(Rcpp::Named("theta") = theta, 
                            Rcpp::Named("alpha") = a,
                            Rcpp::Named('Acceptance') = accept/reps);
}