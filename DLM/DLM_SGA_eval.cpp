// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double normDens(double x, double mu, double sigmaSq){
  return 1 / sqrt(2*3.14159*sigmaSq) * exp(-(x - mu)*(x - mu) / (2*sigmaSq));
}

// Update VB and take logscore of next observed y J many times
// Repeat for a range of max update iterations - see how cheap the update needs to be
// To add: code for restricted updates, will require the SGA algorithm itself to update
// Should add (or a new function): code for full rank approximation
// Major change is simulation with dependencies - could simulate whole vector, 
// or extract relevant parts from covariance matrix - requires an extra TXT matrix multiplication and 5x5 chol decomp
// [[Rcpp::export]]
mat ForecastEvalMF(vec y, vec ySupport, vec maxIter, int T, int J, string updates){
  List initialFit = DLM_SGA(y[T]); // fix this
  int M = 100;
  int P = ySupport.n_elem;
  int N = maxIter.n_elem;
  vec output(J, N);
  
  if(updates == "full"){
    for(int i = 0; i < N; ++i){
      for(int t = T; t < T+J; ++t){
        // Update the VB fit
        if(t > T){
          List updatedFit = DLM_SGA(y[T+t], maxIter[i]);
        }
        // Draw posterior samples
        vec phi = phimean + phisd * randn<vec>(M);
        vec mu = mumean + musd * randn<vec>(M);
        vec sigmaSqX = exp(sigXmean + sigXsd * randn<vec>(M));
        vec sigmaSqY =  exp(sigYmean + sigYsd * randn<vec>(M));
        vec xT = xTmean + xTsd * randn<vec>(M);
        vec xFuture(M);
        vec forecastDensity(P, fill::zeros);
      
        // Forecast future X value, use to form forecast density for future Y
        for(int m = 0; m < M; ++m){
        xFuture[m] = phi[m]*xT[m] + sqrt(sigmqSqX[m])*randn<vec>(1)[0]; 
        for(int p = 0; p < P; ++p){
          forecastDensity[p] += normDens(ySupport[p], xFuture[m] + mu[m], sqrt(sigmaSqY[m]));
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
      output(t-T, i) = log(linearInterpolate);
      }
    }
  }
  // Put in other updating methods
  // Will need to adjust VB algorithms
  return output;
}


