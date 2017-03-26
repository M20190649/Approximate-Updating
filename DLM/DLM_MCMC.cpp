// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
rowvec FFBScpp(vec y, rowvec theta){
  int T = y.size();
  //subtract the mean from y
  y -= theta[3];
  //forward filter steps
  vec att(T, fill::ones); //0 - (T-1) -> x1 - xT
  vec ptt(T, fill::ones);
  rowvec draws(T+1); //0 - T -> x0 - xT
  double at;
  double pt = pow(theta[2], 2)*theta[1] + theta[1];
  double vt = y[0];
  att[0] = pt * vt / (pt + theta[0]);
  ptt[0] = pt - pow(pt, 2) / (pt + theta[0]);
  for(int t = 1; t < T; ++t){
    at = theta[2]*att[t-1];
    pt = pow(theta[2], 2)*ptt[t-1] + theta[1];
    vt = y[t] - at;
    att[t] = at + pt * vt / (pt + theta[0]);
    ptt[t] = pt - pow(pt, 2) / (pt + theta[0]);
  }   
  //backwards sampler steps
  double vstar;
  double fstar;
  double mstar;
  vec atT(T);
  vec ptT(T);
  draws[T] = att[T-1] + sqrt(ptt[T-1])*randn<vec>(1)[0];
  for(int t = T - 1; t > 0; --t){
    vstar = draws[t+1] - theta[2]*att[t-1];
    fstar = pow(theta[2], 2)*ptt[t-1] + theta[1];
    mstar = ptt[t-1]*theta[2];
    atT[t] = att[t-1] + mstar * vstar / fstar;
    ptT[t] = ptt[t-1] - pow(mstar, 2) / fstar;
    draws[t] = atT[t] + sqrt(ptT[t])*randn<vec>(1)[0];
  }
  atT[0] = theta[2]*draws[1] / (pow(theta[2], 2) + 1);
  ptT[0] = theta[1] / (pow(theta[2], 2) + 1);
  draws[0] = atT[0] + sqrt(ptT[0]) * randn<vec>(1)[0];
  
  return draws;
}

// [[Rcpp::export]]
Rcpp::List DLM_MCMC(vec y, int reps){
  int T = y.n_elem;
  mat theta(reps, 4);
  mat x(reps, T+1);
  rowvec initial = {0, 0, 1, 1};
  theta.row(0) = initial;
  
  double sigmaSqYShape = T/2 + 1;
  double sigmaSqXShape = T/2 + 1.5;
  double sigmaSqYScale;
  double sigmaSqXScale;
  double meanPhiNumer;
  double meanPhiDenom;
  double meanMuNumer;
  double meanMuDenom;
  double draw;
  bool stationary;
  
  for(int i = 1; i < reps; ++i){
    x.row(0) = FFBScpp(y, theta.row(i-1));
    
    //sigmaSqY ~ IG(shape, scale)
    sigmaSqYScale = 1;
    for(int t = 0; t < T; ++t){
      sigmaSqYScale += pow(y[t] - x(i, t+1) - theta(i-1, 3), 2) / 2;
    }
    theta(i, 0) = (1 / randg<vec>(1, distr_param(sigmaSqYShape, 1/sigmaSqYScale)))[0];
    
    //sigmaSqX ~ IG(shape, scale)
    sigmaSqXScale = 1 + pow(x(i, 0), 2) / 2;
    for(int t = 0; t < T; ++t){
      sigmaSqXScale += pow(x(i, t+1) - theta(i-1, 2)*x(i, t), 2) / 2;
    }
    theta(i, 1) = (1 / randg<vec>(1, distr_param(sigmaSqXShape, 1/sigmaSqXScale)))[0];
    
    //phi ~ TruncNormal(mean, var, low = -1, high = 1)
    meanPhiNumer = meanPhiDenom = 0;
    for(int t = 0; t < T; ++t){
      meanPhiNumer += x(i, t) * x(i, t+1);
      meanPhiDenom += pow(x(i, t), 2);
    }
    stationary = FALSE;
    while(!stationary){
      draw = (meanPhiNumer/meanPhiDenom + sqrt(theta(i, 1)/meanPhiDenom) * randn<vec>(1))[0];
      if(draw > -1 & draw < 1){
        stationary = TRUE;
        theta(i, 2) = draw;
      }
    }
    
    //mu ~ Normal(mean, var)
    meanMuNumer = 0;
    meanMuDenom = 10*T + theta(i, 0);
    for(int t = 0; t < T; ++t){
      meanMuNumer += 10 * y[t] - x(i, t+1);
    }
    theta(i, 3) = (meanMuNumer/meanMuDenom + sqrt(10*theta(i, 0)/meanMuDenom) * randn<vec>(1))[0];
  }
  
  return Rcpp::List::create(theta = theta, x = x);
}

double normDens(double x, double mu, double sigmaSq){
  return 1 / sqrt(2*3.14159*sigmaSq) * exp(-(x - mu)*(x - mu) / (2*sigmaSq));
}

// This function seems pretty intensive.
// The aim is to calculate the log score for J many future values of y
// This is repeated at a range of different MCMC chain lengths supplied by the chainLength argument
// The vector y should be only the future values of y
// matrices x and theta should have at least max(chainlength) rows and contain the MCMC samples
// ySupport is a sequence of values to evaluate the density at
// burnRatio supplies the proportion of draws to be discarded
// [[Rcpp::export]]
mat ForecastEvalMCMC(vec y, vec chainLength, mat x, mat theta, vec ySupport, double burnRatio){
  int J = y.n_elem;
  int N = chainLength.n_elem;
  int M = 100;
  int P = ySupport.n_elem;
  mat output(J, N);
  
  // For different MCMC chain lengths
  for(int i = 0; i < N; ++i){ 
    int reps = chainLength[i];
    mat xfuture(M, J);
    // and J many forecast samples
    for(int t = 0; t < J; ++t){ 
      vec forecastDensity(P, fill::zeros);
      int yDensityIndex;
      // draw M many posterior forecast evaluations
      for(int m = 0; m < M; ++m){
        double xT;
        
        // Draw random samples from posterior, ignoring first 20% of draws
        int u = rand() % reps*(1-burnRatio) + burnRatio*reps;
        double sigmaSqY = theta(u, 0);
        u = rand() % reps*(1-burnRatio)+ burnRatio*reps;
          double sigmaSqX = theta(u, 1);
        u = rand() % reps*(1-burnRatio) + burnRatio*reps;
          double phi = theta(u, 2);
        u = rand() % reps*(1-burnRatio) + burnRatio*reps;
          double mu = theta(u, 3);
        // X_T+1 based off XT
        if(t == 0){ 
          u = rand() % reps*(1-burnRatio) + burnRatio*reps;
          xT = x(u, x.n_cols - 1);
        // X_T+j from X_T+j-1
        } else {
          u = rand() % M;
          xT = xfuture(u, t-1);
        }
        
        double x_Tj = phi*xT + sqrt(sigmaSqX) * randn<vec>(1)[0];
        
        // Now evaluate the forecast density over a range of values in the support of y
        for(int p = 0; p < P; ++p){
          forecastDensity[p] += normDens(ySupport[p], x_Tj+mu, sigmaSqY)/P;
        }
      }
      
      // Find the first index of the forecast density that is lower than observed y
      for(int p = 0; p < P; ++p){
        if(forecastDensity[p] < y[t]){
          yDensityIndex = p;
          break;
        }
      }
      // Linear interpolation between the two forecast density points 
      double linearInterpolate = forecastDensity[yDensityIndex] + (y[t] - ySupport[yDensityIndex]) *
          (forecastDensity[yDensityIndex+1] - forecastDensity[yDensityIndex]) / (ySupport[yDensityIndex+1] - ySupport[yDensityIndex]);
      
      // Finally save the log score for the t'th forecast and the i'th chain length
      output(t, i) = log(linearInterpolate);
    }
  }
  return output;
  
}
