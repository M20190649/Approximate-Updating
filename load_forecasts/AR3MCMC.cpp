// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

//Returns sum of element by element products of two vectors, R style
double vecprodsum(vec a, vec b){
  int N = a.n_elem;
  double out = 0;
  for(int i = 0; i < N; ++i){
    out += a[i]*b[i];
  }
  return out;
}

double rowprodsum(vec a, rowvec b){
  int N = a.n_elem;
  double out = 0;
  for(int i = 0; i < N; ++i){
    out += a[i]*b[i];
  }
  return out;
}

// [[Rcpp::export]]
vec yy(vec y, int T){
  vec lags = {0, 1, 48, 336}; //lags used in calculations
  vec sumyy(9, fill::zeros); //rows: yt yt -1 - yt yt-48 - yt yt-336 - yt-1 yt-1 - ... - yt-336 yt-336
  for(int i =  336; i < T; ++i){ //iterate over time
    for(int k = 1; k < 4; ++k){ //y t - all lags
      sumyy[k-1] += y[i]*y[i-lags[k]];
      sumyy[2+k] += y[i-1]*y[i-lags[k]];
    }
    sumyy[6] += y[i-48]*y[i-48];
    sumyy[7] += y[i-48]*y[i-336];
    sumyy[8] += y[i-336]*y[i-336];
  }
  return sumyy;
}

// [[Rcpp::export]]
mat yx(vec y, mat x, int T, int dimx){
  vec lags = {0, 1, 48, 336}; 
  mat sumyx(dimx, 16, fill::zeros); //cols: yt xt - yt xt-1 - yt xt-2 - ... - yt-336 xt-336
  for(int i =  336; i < T; ++i){ //iterate over time
    for(int k = 0; k < 4; ++k){ //y t - all lags
      for(int l = 0; l < 4; ++l){ //x t - all lags
        for(int j = 0; j < dimx; ++j){ 
          sumyx(j, k*4 + l) += y[i - lags[k]]*x(i - lags[l], j);
        }
      }
    }
  }
  return sumyx;
}

// [[Rcpp::export]]
cube xx(mat x, int T, int dimx){
  vec lags = {0, 1, 48, 336}; 
  cube sumxx(dimx, dimx, 10, fill::zeros); //slices xt xt - xt xt-1 - ..  xt-336 xt-336
  for(int i = 336; i < T; ++i){ //time
    for(int j = 0; j < dimx; ++j){ //x1
      for(int k = 0; k < dimx; ++k){ //x2
        for(int l = 0; l < 4; ++l){ //lags
          sumxx(j, k, l) += x(i, j) * x(i - lags[l], k);
          if(l > 0){
            sumxx(j, k, l + 3) += x(i - 1, j) * x(i - lags[l], k);
            if(l > 1){
              sumxx(j, k, l + 5) += x(i - 48, j) * x(i - lags[l], k);
            }
          }
        }
        sumxx(j, k, 9) += x(i - 336, j) * x(i - 336, k);
      }
    }
  }
  return sumxx;
}

// [[Rcpp::export]]
mat MCMC_looponly(int rep, vec y, mat x, int T, int dimx, vec sumyy, mat sumyx, cube sumxx,  vec beta, vec phi, double sigmasq,
                 vec mu, vec tau, vec lambda, vec gamma, double alphabar, double betasig){
  
  mat store(rep,  dimx + 4);
  
  for(int r = 0; r < rep; ++r){
    
    //require xt b sums for a given beta
    vec bxsq(9, fill::zeros); //xt xt-1 ---- xt - 336 xt -336 (omit xt xt term)
    for(int k = 0; k < dimx; ++k){
      for(int j = 0; j < dimx; ++j){
        for(int i = 0; i < 9; ++i){
          bxsq[i] += beta[k] * beta[j] * sumxx(k, j, i+1);
        }
      }  
    }
    
    //Phi Conditionals
    
    double numer = 0;
    double denom = 0;
    double var = 0;
    
    //Phi 1
    numer = lambda[0] * (sumyy[0] + bxsq[0] - vecprodsum(beta,(sumyx.col(4) + sumyx.col(1))) - 
      phi[1] * (sumyy[4] + bxsq[4] - vecprodsum(beta, (sumyx.col(6) + sumyx.col(9)))) -
      phi[2] * (sumyy[5] + bxsq[5] - vecprodsum(beta, (sumyx.col(7) + sumyx.col(13))))) + sigmasq * gamma[0];
    denom = lambda[0] * (sumyy[3] + bxsq[3] - 2 * vecprodsum(beta, sumyx.col(5))) + sigmasq;
    var = sigmasq * lambda[0] / denom;
    phi[0] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Phi 2
    numer = lambda[1] * (sumyy[1] + bxsq[1] - vecprodsum(beta, (sumyx.col(8) + sumyx.col(2))) - 
      phi[0] * (sumyy[4] + bxsq[4] - vecprodsum(beta, (sumyx.col(6) + sumyx.col(9)))) -
      phi[2] * (sumyy[7] + bxsq[7] - vecprodsum(beta, (sumyx.col(11) + sumyx.col(14))))) + sigmasq * gamma[1];
    denom = lambda[1] * (sumyy[6] + bxsq[6] - 2 * vecprodsum(beta, sumyx.col(10))) + sigmasq;
    var = sigmasq * lambda[1] / denom;
    phi[1] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Phi 3
    numer = lambda[2] * (sumyy[2] + bxsq[2] - vecprodsum(beta, (sumyx.col(12) + sumyx.col(3))) - 
      phi[0] * (sumyy[5] + bxsq[5] - vecprodsum(beta, (sumyx.col(7) + sumyx.col(13)))) -
      phi[1] * (sumyy[7] + bxsq[7] - vecprodsum(beta, (sumyx.col(11) + sumyx.col(14))))) + sigmasq * gamma[2];
    denom = lambda[2] * (sumyy[8] + bxsq[8] - 2 * vecprodsum(beta, sumyx.col(15))) + sigmasq;
    var = sigmasq * lambda[2] / denom;
    phi[2] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Beta Conditionals
    
    for(int j = 0; j < dimx; ++j){ 
      
      double numera = 0;
      double numerb = 0;
      
      numera = sumyx(j, 0) - phi[0] * sumyx(j, 1) - phi[1] * sumyx(j, 2) - phi[2] * sumyx(j, 3) -
        phi[0] * sumyx(j, 4) + phi[0] * phi[0] * sumyx(j, 5) + phi[0] * phi[1] * sumyx(j, 6) + phi[0] * phi[2] * sumyx(j, 7) -
        phi[1] * sumyx(j, 8) + phi[1] * phi[0] * sumyx(j, 9) + phi[1] * phi[1] * sumyx(j, 10) + phi[1] * phi[2] * sumyx(j, 11) -
        phi[2] * sumyx(j, 12) + phi[2] * phi[0] * sumyx(j, 13) + phi[2] * phi[1] * sumyx(j, 14) + phi[2] * phi[2] * sumyx(j, 15);
      
      for(int k = 0; k < dimx; ++k){
        //sum is over all parts of beta except beta j
        if(k != j){
          numerb += beta[k] * (sumxx(j, k, 0) - phi[0] * (sumxx(j, k, 1) + sumxx(k, j, 1)) - phi[1] * (sumxx(j, k, 2) + sumxx(k, j, 2)) -
            phi[2] * (sumxx(j, k, 3) + sumxx(k, j, 3)) + phi[0] * phi[0] * sumxx(j, k, 4) + phi[0] * phi[1] * (sumxx(j, k, 5) + sumxx(k, j, 5)) +
            phi[0] * phi[2] * (sumxx(j, k, 6) + sumxx(k, j, 6)) + phi[1] * phi[1] * sumxx(j, k, 7) + 
            phi[1] * phi[2] * (sumxx(j, k, 8) + sumxx(k, j, 8)) + phi[2] * phi[2] * sumxx(j, k, 9));
        }
      }
      
      numer = tau[j] * (numera - numerb) + sigmasq*mu[j];
      denom = tau[j] * (sumxx(j, j, 0) + phi[0] * phi[0]*sumxx(j, j, 4) + phi[1] * phi[1] * sumxx(j, j, 7) + phi[2] * phi[2] * sumxx(j, j, 9) - 
        2 * (phi[0] * sumxx(j, j, 1) + phi[1] * sumxx(j, j, 2) + phi[2] * sumxx(j, j, 3) - phi[0] * phi[1] * sumxx(j, j, 5) -
        phi[0] * phi[2] * sumxx(j, j, 6) - phi[1] * phi[2] * sumxx(j, j, 8))) + sigmasq;
      var = tau[j] * sigmasq / denom;
      
      beta[j] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    }
    
    //Sigma squared conditional
    double betabar;
    double error = 0;
    double sse = 0;
    //Need to add half the sum of squared errors to betabar
    for(int i = 336; i < T; ++i){
      error = y[i] - rowprodsum(beta, x.row(i)) - phi[0]* (y[i - 1] - rowprodsum(beta, x.row(i - 1))) - 
        phi[1] * (y[i - 48] - rowprodsum(beta, x.row(i - 48))) -  phi[2] * (y[i - 336] - rowprodsum(beta, x.row(i - 336)));
      sse += error * error;
    }
    betabar = betasig + 1/2 * sse;
    //betabar is gamma rate, randg takes scale 
    sigmasq = randg<vec>(1, distr_param(alphabar, 1 / betabar))[0];
    
    for(int j = 0; j < dimx; ++j){
      store(r, j) = beta[j];
    }
    for(int j = 0; j < 3; ++j){
      store(r, j + dimx) = phi[j];
    }
    store(r, dimx + 3) = sigmasq;
  }
  return store;
}

// [[Rcpp::export]]
mat MCMC_load(vec y, mat x, int rep, int burn = 0, int thin = 1) {
  int dimx = x.n_cols;
  int T = y.n_elem;
  
  // Prior parameters
  
  //Beta - 0 intercept,  1-47 half hour effects, 48-53 day of week, 54 temp, dimx linear time trend, base - sunday 12:00
  vec mu(dimx); mu.fill(0);
  vec tau(dimx); tau.fill(100000000);
  //Phi 
  vec gamma(3); gamma.fill(0);
  vec lambda(3); lambda.fill(10);
  //Sigma Squared
  double alphasig = 1;
  double betasig = 1;
  //Doesnt depend on parameters
  double alphabar = alphasig + T/2 - 168;
  double betabar;

  // Preliminary Data Calculations
  
  //Yt lags * Yt lags
  vec lags = {0, 1, 48, 336}; //lags used in calculations
  vec sumyy(9, fill::zeros); //rows: yt yt -1 - yt yt-48 - yt yt-336 - yt-1 yt-1 - ... - yt-336 yt-336
  for(int i =  336; i < T; ++i){ //iterate over time
    for(int k = 1; k < 4; ++k){ //y t - all lags
      sumyy[k-1] += y[i]*y[i-lags[k]];
      sumyy[2+k] += y[i-1]*y[i-lags[k]];
    }
    sumyy[6] += y[i-48]*y[i-48];
    sumyy[7] += y[i-48]*y[i-336];
    sumyy[8] += y[i-336]*y[i-336];
  }

  //Xt lags * Yt lags
  mat sumyx(dimx, 16, fill::zeros); //cols: yt xt - yt xt-1 - yt xt-2 - ... - yt-336 xt-336
  for(int i =  336; i < T; ++i){ //iterate over time
    for(int k = 0; k < 4; ++k){ //y t - all lags
      for(int l = 0; l < 4; ++l){ //x t - all lags
        for(int j = 0; j < dimx; ++j){ 
          sumyx(j, k*4 + l) += y[i - lags[k]]*x(i - lags[l], j);
        }
      }
    }
  }
  
  //Xt lags * Xt lags
  //Beta X square terms require all x cross products. This is a better approach than a somewhat simpler calculation inside the main loop
  cube sumxx(dimx, dimx, 10, fill::zeros); //slices xt xt - xt xt-1 - ..  xt-336 xt-336
  for(int i = 336; i < T; ++i){ //time
    for(int j = 0; j < dimx; ++j){ //x1
      for(int k = 0; k < dimx; ++k){ //x2
        for(int l = 0; l < 4; ++l){ //lags
          sumxx(j, k, l) += x(i, j) * x(i - lags[l], k);
          if(l > 0){
            sumxx(j, k, l + 3) += x(i - 1, j) * x(i - lags[l], k);
            if(l > 1){
              sumxx(j, k, l + 5) += x(i - 48, j) * x(i - lags[l], k);
            }
          }
        }
      sumxx(j, k, 9) += x(i - 336, j) * x(i - 336, k);
      }
    }
  }
  
  double keep = (rep - burn)/thin; //number of kept draws
  mat store(floor(keep), dimx + 4); //one row per kept draw, eventual output matrix
  int keepiter = 0; //track kept row
  
  //initial conditions
  vec beta(dimx, fill::zeros); 
  vec phi(3, fill::zeros);
  double sigmasq = 1;
  
  //Main Loop
  
  for(int r = 0; r < rep; ++r){
    
    //require xt b sums for a given beta
    vec bxsq(9, fill::zeros); //xt xt-1 ---- xt - 336 xt -336 (omit xt xt term)
    for(int k = 0; k < dimx; ++k){
      for(int j = 0; j < dimx; ++j){
        for(int i = 0; i < 9; ++i){
          bxsq[i] += beta[k] * beta[j] * sumxx(k, j, i+1);
        }
      }  
    }
    
    //Phi Conditionals
    
    double numer = 0;
    double denom = 0;
    double var = 0;
    
    //Phi 1
    numer = lambda[0] * (sumyy[0] + bxsq[0] - vecprodsum(beta,(sumyx.col(4) + sumyx.col(1))) - 
      phi[1] * (sumyy[4] + bxsq[4] - vecprodsum(beta, (sumyx.col(6) + sumyx.col(9)))) -
      phi[2] * (sumyy[5] + bxsq[5] - vecprodsum(beta, (sumyx.col(7) + sumyx.col(13))))) + sigmasq * gamma[0];
    denom = lambda[0] * (sumyy[3] + bxsq[3] - 2 * vecprodsum(beta, sumyx.col(5))) + sigmasq;
    var = sigmasq * lambda[0] / denom;
    phi[0] = numer / denom + sqrt(var) * randn<vec>(1)[0];
      
    //Phi 2
    numer = lambda[1] * (sumyy[1] + bxsq[1] - vecprodsum(beta, (sumyx.col(8) + sumyx.col(2))) - 
      phi[0] * (sumyy[4] + bxsq[4] - vecprodsum(beta, (sumyx.col(6) + sumyx.col(9)))) -
      phi[2] * (sumyy[7] + bxsq[7] - vecprodsum(beta, (sumyx.col(11) + sumyx.col(14))))) + sigmasq * gamma[1];
    denom = lambda[1] * (sumyy[6] + bxsq[6] - 2 * vecprodsum(beta, sumyx.col(10))) + sigmasq;
    var = sigmasq * lambda[1] / denom;
    phi[1] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Phi 3
    numer = lambda[2] * (sumyy[2] + bxsq[2] - vecprodsum(beta, (sumyx.col(12) + sumyx.col(3))) - 
      phi[0] * (sumyy[5] + bxsq[5] - vecprodsum(beta, (sumyx.col(7) + sumyx.col(13)))) -
      phi[1] * (sumyy[7] + bxsq[7] - vecprodsum(beta, (sumyx.col(11) + sumyx.col(14))))) + sigmasq * gamma[2];
    denom = lambda[2] * (sumyy[8] + bxsq[8] - 2 * vecprodsum(beta, sumyx.col(15))) + sigmasq;
    var = sigmasq * lambda[2] / denom;
    phi[2] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Beta Conditionals

    for(int j = 0; j < dimx; ++j){ 
      
      double numera = 0;
      double numerb = 0;
      
      numera = sumyx(j, 0) - phi[0] * sumyx(j, 1) - phi[1] * sumyx(j, 2) - phi[2] * sumyx(j, 3) -
        phi[0] * sumyx(j, 4) + phi[0] * phi[0] * sumyx(j, 5) + phi[0] * phi[1] * sumyx(j, 6) + phi[0] * phi[2] * sumyx(j, 7) -
        phi[1] * sumyx(j, 8) + phi[1] * phi[0] * sumyx(j, 9) + phi[1] * phi[1] * sumyx(j, 10) + phi[1] * phi[2] * sumyx(j, 11) -
        phi[2] * sumyx(j, 12) + phi[2] * phi[0] * sumyx(j, 13) + phi[2] * phi[1] * sumyx(j, 14) + phi[2] * phi[2] * sumyx(j, 15);

      for(int k = 0; k < dimx; ++k){
        //sum is over all parts of beta except beta j
        if(k != j){
          numerb += beta[k] * (sumxx(j, k, 0) - phi[0] * (sumxx(j, k, 1) + sumxx(k, j, 1)) - phi[1] * (sumxx(j, k, 2) + sumxx(k, j, 2)) -
            phi[2] * (sumxx(j, k, 3) + sumxx(k, j, 3)) + phi[0] * phi[0] * sumxx(j, k, 4) + phi[0] * phi[1] * (sumxx(j, k, 5) + sumxx(k, j, 5)) +
            phi[0] * phi[2] * (sumxx(j, k, 6) + sumxx(k, j, 6)) + phi[1] * phi[1] * sumxx(j, k, 7) + 
            phi[1] * phi[2] * (sumxx(j, k, 8) + sumxx(k, j, 8)) + phi[2] * phi[2] * sumxx(j, k, 9));
        }
      }
      
      numer = tau[j] * (numera - numerb) + sigmasq*mu[j];
      denom = tau[j] * (sumxx(j, j, 0) + phi[0] * phi[0]*sumxx(j, j, 4) + phi[1] * phi[1] * sumxx(j, j, 7) + phi[2] * phi[2] * sumxx(j, j, 9) - 
        2 * (phi[0] * sumxx(j, j, 1) + phi[1] * sumxx(j, j, 2) + phi[2] * sumxx(j, j, 3) - phi[0] * phi[1] * sumxx(j, j, 5) -
        phi[0] * phi[2] * sumxx(j, j, 6) - phi[1] * phi[2] * sumxx(j, j, 8))) + sigmasq;
      var = tau[j] * sigmasq / denom;
      
      beta[j] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    }
    
    //Sigma squared conditional
    
    double error = 0;
    double sse = 0;
    //Need to add half the sum of squared errors to betabar
    for(int i = 336; i < T; ++i){
      error = y[i] - rowprodsum(beta, x.row(i)) - phi[0]* (y[i - 1] - rowprodsum(beta, x.row(i - 1))) - 
        phi[1] * (y[i - 48] - rowprodsum(beta, x.row(i - 48))) -  phi[2] * (y[i - 336] - rowprodsum(beta, x.row(i - 336)));
      sse += error * error;
    }
    betabar = betasig + 1/2 * sse;
    //betabar is gamma rate, randg takes scale 
    sigmasq = randg<vec>(1, distr_param(alphabar, 1 / betabar))[0];
    
    //Store draws
    if(r >= burn & (r - burn) % thin == 0){
      for(int j = 0; j < dimx; ++j){
        store(keepiter, j) = beta[j];
      }
      for(int j = 0; j < 3; ++j){
        store(keepiter, j + dimx) = phi[j];
      }
      store(keepiter, dimx + 3) = sigmasq;
      keepiter += 1;
    }
  }
  return store;
}