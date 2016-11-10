// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include <vector>
#include <functional>   
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat MCMC.load(vec y, mat x, int burn, int rep, int thin, int T) {
  
  // Prior parameters
  
  //Beta - 0 intercept, 1 temperature, 2-7 weekdays, 8-54 half hour effects, base - monday 12:00
  vec mu(55); mu.fill(0);
  vec tau(55); tau.fill(10);
  //Phi 
  vec gamma(3); gamma.fill(0);
  vec lambda(3); lambda.fill(10);
  //Sigma Squared
  double alphasig = 1;
  double betasig = 1;
  //Doesnt depend on parameters
  double alphabar = alphasig + T - 336;
  
  // Preliminary Data Calculations
  
  //Yt lags * Yt lags
  vec lags = {0, 1, 48, 336}; //lags used in calculations
  vec sumyy(9, fill::zeroes); //rows: yt yt -1 - yt yt-48 - yt yt-336 - yt-1 yt-1 - ... - yt-336 yt-336
  for(int i =  336; i < T; i++){ //iterate over time
    for(int k = 1; k < 4; k++){ //y t - all lags
      sumyy[k-1] += y[i]*y[i-lags[k]];
      sumyy[2+k] += y[i-1]*y[i-lags[k]];
    }
    sumyy[6] += y[i-48]*y[i-48];
    sumyy[7] += y[i-48]*y[i-336];
    sumyy[8] += y[i-336]*y[i-336];
  }

  //Xt lags * Yt lags
  mat sumyx(55, 16, fill::zeroes); //cols: yt xt - yt xt-1 - yt xt-2 - ... - yt-336 xt-336
  for(int i =  336; i < T; i++){ //iterate over time
    for(int k = 0; k < 4; k++){ //y t - all lags
      for(int l = 0; l < 4; l++){ //x t - all lags
        for(int j = 0; j < 55; j++){ 
          sumyx(j, k*4 + l) += y[i - lags[k]]*x(i - lags[l], j);
        }
      }
    }
  }
  
  //Xt lags * Xt lags
  mat sumxx(55, 10, fill::zeroes); //cols xt xt - xt xt-1 - ..  xt-336 xt-336
  for(int i = 336; i < T; i++){
    for(int j = 0; j < 55; j++){
      for(int k = 0; k < 4; k++){
        sumxx(j, k) += x(i, j) * x(i - lags[k], j);
        if(k > 0){
          sumxx(j, k + 3) += x(i - 1, j) * x(i - lags[k], j);
          if(k > 1){
            sumyx(j, k + 5) += x(i - 48, j) * x(i - lags[k], j);
          }
        }
      }
      sumxx(j, 9) += x(i - 336, j) * x(i - 336, j);
    }
  }
  
  double keep = (rep - burn)/thin; //number of kept draws
  mat store(floor(keep), 59); //one row per kept draw, eventual output matrix
  int keepiter = 0; //track kept row
  
  //initial conditions
  vec beta(55, fill::zeroes); 
  vec phi(3, fill::zeroes);
  double sigmasq = 1;
  
  //Main Loop
  
  for(int r = 1; r < rep; r++){
    
    //require xt b sums for a given beta
    vec sumbx(4, fill::zeroes)
    for(int i = 336; i < T; i++){
      for(int j = 0; j < 55; j++){
        for(int k = 0; k < 4; k++){
          sumbx[k] += beta[j]*x(i - lags[k], j);
        }
      }  
    }
    //Phi Conditionals
    double denom = 0;
    double numer = 0;
    double var = 0;
    //Phi 1
    numer = lambda[0]*(sumyy[0] + sumbx[0]*sumbx[1] - sum(beta * (sumyx.col(4) + sumyx.col(1))) + 
        phi[1]*(sumyy[4] + sumbx[1]*sumbx[2] - sum(beta * (sumyx.col(6) + sumyx.col(9)))) +
        phi[2]*(sumyy[5] + sumbx[1]*sumbx[3] - sum(beta * (sumyx.col(7) + sumyx.col(13))))) + sigmasq*gamma[0];
    denom = lambda[0]*(sumyy[3] + sumbx[1]*sumbx[1] - 2 * sum(beta * sumyx.col(5))) + sigmasq;
    var = (sigmasq*lambda[0])/denom;
    phi[0] = numer/denom + sqrt(var)*randn<vec>(1)[0];
      
    //Phi 2
    numer = lambda[1]*(sumyy[1] + sumbx[0]*sumbx[2] - sum(beta * (sumyx.col(8) + sumyx.col(2))) + 
      phi[0]*(sumyy[4] + sumbx[1]*sumbx[2] - sum(beta * (sumyx.col(6) + sumyx.col(9)))) +
      phi[2]*(sumyy[7] + sumbx[2]*sumbx[3] - sum(beta * (sumyx.col(11) + sumyx.col(15))))) + sigmasq*gamma[1];
    denom = lambda[1]*(sumyy[6] + sumbx[2]*sumbx[2] - 2 * sum(beta * sumyx.col(10))) + sigmasq;
    var = (sigmasq*lambda[1])/denom;
    phi[1] = numer/denom + sqrt(var)*randn<vec>(1)[0];
    
    //Phi 3
    numer = lambda[2]*(sumyy[2] + sumbx[0]*sumbx[3] - sum(beta * (sumyx.col(12) + sumyx.col(3))) + 
      phi[0]*(sumyy[5] + sumbx[1]*sumbx[3] - sum(beta * (sumyx.col(7) + sumyx.col(13)))) +
      phi[1]*(sumyy[7] + sumbx[2]*sumbx[3] - sum(beta * (sumyx.col(11) + sumyx.col(15))))) + sigmasq*gamma[2];
    denom = lambda[2]*(sumyy[8] + sumbx[3]*sumbx[3] - 2 * sum(beta * sumyx.col(15))) + sigmasq;
    var = (sigmasq*lambda[2])/denom;
    phi[2] = numer/denom + sqrt(var)*randn<vec>(1)[0];
    
    //Beta Conditionals
    double numera = 0;
    double numerb = 0;
    
    //Each j has a common component to the mean denominator
    double commondenom = sumxx[0] + phi[0]*phi[0]*sumxx[4] + phi[1]*phi[1]*sumxxx[7] - 
      2 * (phi[0] * sumxx[1] + phi[1] * sumxx[2] + phi[2] * sumxx[3] - phi[0]*phi[1]*sumxx[5] -
      phi[0]*phi[2]*sumxx[6] - phi[1]*phi[2]*sumxx[8]);
    
    //Draw for each beta j
    for(int j = 0; j < 55; j++){ 
      numera = sumyx(j, 0) - phi[0]*sumyx(j, 1) - phi[1]*sumyx(j, 2) - phi[2]*sumyx(j, 3) -
        phi[0]*sumyx(j, 4) + phi[0]*phi[0]*sumyx(j, 5) + phi[0]*phi[1]*sumyx(j, 6) + phi[0]*phi[2]*sumyx(j, 7) -
        phi[1]*sumyx(j, 8) + phi[1]*phi[0]*sumyx(j, 9) + phi[1]*phi[1]*sumyx(j, 10) + phi[1]*phi[2]*sumyx(j, 11) -
        phi[2]*sumyx(j, 12) + phi[2]*phi[0]*sumyx(j, 13) + phi[2]*phi[1]*sumyx(j, 14) + phi[2]*phi[2]*sumyx(j, 15);
      
      //Cant precalculate sum over t for the xkt beta k, k neq j sums
      for(int i = 336; i < T; i++){
        double sum1 = 0;
        double sum2 = 0;
        double sum3 = 0;
        double sum4 = 0;
        for(int k = 0; k < 55; k++){
          //sum is over all parts of beta except beta j
          if(k != j){
            sum1 += x(i, k) * beta(k);
            sum2 += phi[0] * x(i - 1, k) * beta(k);
            sum3 += phi[1] * x(i - 48, k) * beta(k);
            sum4 += phi[2] * x(i - 336, k) * beta(k);
          }
        }
        numerb += (x(i, j) - phi[0] * x(i - 1, j) - phi[1] * x(i - 48, j) - phi[2] * x(i - 336, j)) * (sum1 - sum2 - sum3 - sum4);
      }
        
      numer = tau[j] * (numera - numerb) + sigmasq*mu[j];
      denom = tau[j] * commondenom + sigmasq;
      var = tau[j] * sigmasq / denom;
      
      beta[j] = numer/denom + sqrt(var) * randn<vec>(1)[0];
    }
    
    //Sigma squared conditional
    betabar = betasig;
    double error;
    //Need to add half the sum of squared errors
    for(int i = 336; i < T; i++){
      error = y[i] - sum(beta * x.row(i)) + phi[0]* (y[i - 1] - sum(beta * x.row(i - 1))) + 
        phi[1]* (y[i - 48] - sum(beta * x.row(i - 48))) +  phi[2]* (y[i - 336] - sum(beta * x.row(i - 336)));
      betabar += 1/2 * error * error;
    }
    
    //Draw sigma squared from a gamma distribution here
    
    //store draws
    if(r > burn & r - burn % thin ==0){
      for(int j = 0; j < 55; j++){
        store(keepiter, j) = beta[j];
      }
      for(int j = 0; j < 3; j++){
        store(keepiter, j + 55) = phi[j];
      }
      store(keepiter, 58) = sigmasq[j];
      keepiter += 1
    }
  }
  return store;
}