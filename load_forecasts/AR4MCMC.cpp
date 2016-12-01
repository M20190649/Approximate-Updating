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
mat MCMC_load4(vec y, mat x, int rep,  vec lags, int burn = 0, int thin = 1) {
  int dimx = x.n_cols;
  int T = y.n_elem;
  
  // Prior parameters
  
  //Beta - 0 intercept,  1-47 half hour effects, 48-53 day of week, 54 temp, 55 time trend, base - sunday 12:00
  vec mu(dimx); mu.fill(0);
  vec tau(dimx); tau.fill(1000000);
  //Phi 
  vec gamma(4); gamma.fill(0);
  vec lambda(4); lambda.fill(10);
  //Sigma Squared
  double alphasig = 1;
  double betasig = 1;
  //Doesnt depend on parameters
  double alphabar = alphasig + T/2 - 8928;
  double betabar;
  
  // Preliminary Data Calculations
  
  //Yt lags * Yt lags
  vec sumyy(14, fill::zeros); //rows: yt yt -1 - yt yt-2 - yt yt-48 - ... - yt-336 yt-336
  for(int i = 17856; i < T; ++i){ //iterate over time
    for(int k = 1; k < 5; ++k){ //y t - all lags
      sumyy[k-1] += y[i]*y[i-lags[k]];
      sumyy[3+k] += y[i-1]*y[i-lags[k]];
      if(k > 1){
        sumyy[6+k] +=y[i-2]*y[i-lags[k]];
      }
    }
    sumyy[11] += y[i-336]*y[i-336];
    sumyy[12] += y[i-336]*y[i-17520];
    sumyy[13] += y[i-17520]*y[i-17520];
  }
  
  //Xt lags * Yt lags
  mat sumyx(dimx, 25, fill::zeros); //cols: yt xt - yt xt-1 - yt xt-2 - ... - yt-336 xt-336
  for(int i = 17856; i < T; ++i){ //iterate over time
    for(int k = 0; k < 5; ++k){ //y t - all lags
      for(int l = 0; l < 5; ++l){ //x t - all lags
        for(int j = 0; j < dimx; ++j){ 
          sumyx(j, k*5 + l) += y[i - lags[k]]*x(i - lags[l], j);
        }
      }
    }
  }
  
  //Xt lags * Xt lags
  //Beta X square terms require all x cross products. This is a better approach than a somewhat simpler calculation inside the main loop
  cube sumxx(dimx, dimx, 15, fill::zeros); //slices xt xt - xt xt-1 - ..  xt-336 xt-336
  for(int i = 17856; i < T; ++i){ //time
    for(int j = 0; j < dimx; ++j){ //x1
      for(int k = 0; k < dimx; ++k){ //x2
        for(int l = 0; l < 5; ++l){ //lags
          sumxx(j, k, l) += x(i, j) * x(i - lags[l], k);
          if(l > 0){
            sumxx(j, k, l + 4) += x(i - 1, j) * x(i - lags[l], k);
            if(l > 1){
              sumxx(j, k, l + 7) += x(i - 2, j) * x(i - lags[l], k);
            }
          }
        }
        sumxx(j, k, 12) += x(i - 336, j) * x(i - 336, k);
        sumxx(j, k, 13) += x(i - 336, j) * x(i - 17520, k);
        sumxx(j, k, 14) += x(i - 17520, j) * x(i - 17520, k);
      }
    }
  }
  
  double keep = (rep - burn)/thin; //number of kept draws
  mat store(floor(keep), dimx + 5); //one row per kept draw, eventual output matrix
  int keepiter = 0; //track kept row
  
  //initial conditions
  vec beta(dimx, fill::zeros); 
  vec phi(4, fill::zeros);
  double sigmasq = 1;
  
  //Main Loop
  
  for(int r = 0; r < rep; ++r){
    
    //require xt b sums for a given beta
    vec bxsq(14, fill::zeros); //xt xt-1 ---- xt - 336 xt -336 (omit xt xt term)
    for(int k = 0; k < dimx; ++k){
      for(int j = 0; j < dimx; ++j){
        for(int i = 0; i < 14; ++i){
          bxsq[i] += beta[k] * beta[j] * sumxx(k, j, i+1);
        }
      }  
    }
    
    //Phi Conditionals
    
    double numer = 0;
    double denom = 0;
    double var = 0;
    
    //Phi 1 yt -1 
    numer = lambda[0] * (sumyy[0] + bxsq[0] - vecprodsum(beta,(sumyx.col(5) + sumyx.col(1))) - 
      phi[1] * (sumyy[5] + bxsq[5] - vecprodsum(beta, (sumyx.col(7) + sumyx.col(11)))) -
      phi[2] * (sumyy[6] + bxsq[6] - vecprodsum(beta, (sumyx.col(8) + sumyx.col(16)))) -
      phi[3] * (sumyy[7] + bxsq[7] - vecprodsum(beta, (sumyx.col(9) + sumyx.col(21))))) + sigmasq * gamma[0];
    denom = lambda[0] * (sumyy[4] + bxsq[4] - 2 * vecprodsum(beta, sumyx.col(6))) + sigmasq;
    var = sigmasq * lambda[0] / denom;
    phi[0] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Phi 2 yt - 2
    numer = lambda[1] * (sumyy[1] + bxsq[1] - vecprodsum(beta, (sumyx.col(10) + sumyx.col(2))) - 
      phi[0] * (sumyy[5] + bxsq[5] - vecprodsum(beta, (sumyx.col(7) + sumyx.col(11)))) -
      phi[2] * (sumyy[9] + bxsq[9] - vecprodsum(beta, (sumyx.col(17) + sumyx.col(13)))) -
      phi[3] * (sumyy[10] + bxsq[10] - vecprodsum(beta, (sumyx.col(22) + sumyx.col(14))))) + sigmasq * gamma[1];
    denom = lambda[1] * (sumyy[8] + bxsq[8] - 2 * vecprodsum(beta, sumyx.col(12))) + sigmasq;
    var = sigmasq * lambda[1] / denom;
    phi[1] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Phi 3 yt-48
    numer = lambda[2] * (sumyy[2] + bxsq[2] - vecprodsum(beta, (sumyx.col(15) + sumyx.col(3))) - 
      phi[0] *  (sumyy[6] + bxsq[6] - vecprodsum(beta, (sumyx.col(8) + sumyx.col(16)))) -
      phi[1] * (sumyy[9] + bxsq[9] - vecprodsum(beta, (sumyx.col(17) + sumyx.col(13)))) -
      phi[3] * (sumyy[12] + bxsq[12] - vecprodsum(beta, (sumyx.col(19) + sumyx.col(23))))) + sigmasq * gamma[2];
    denom = lambda[2] * (sumyy[11] + bxsq[11] - 2 * vecprodsum(beta, sumyx.col(18))) + sigmasq;
    var = sigmasq * lambda[2] / denom;
    phi[2] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Phi 4 yt-336
    numer = lambda[3] * (sumyy[3] + bxsq[3] - vecprodsum(beta, (sumyx.col(20) + sumyx.col(4))) - 
      phi[0] * (sumyy[7] + bxsq[7] - vecprodsum(beta, (sumyx.col(9) + sumyx.col(21)))) -
      phi[1] * (sumyy[10] + bxsq[10] - vecprodsum(beta, (sumyx.col(22) + sumyx.col(14)))) -
      phi[2] * (sumyy[12] + bxsq[12] - vecprodsum(beta, (sumyx.col(19) + sumyx.col(23)))))+ sigmasq * gamma[3];
    denom = lambda[3] * (sumyy[13] + bxsq[13] - 2 * vecprodsum(beta, sumyx.col(24))) + sigmasq;
    var = sigmasq * lambda[2] / denom;
    phi[3] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    
    //Beta Conditionals
    
    for(int j = 0; j < dimx; ++j){ 
      
      double numera = 0;
      double numerb = 0;
      vec phicoefs = {1, -phi[0], -phi[1], -phi[2], -phi[3]};

      for(int m = 0; m < 5; ++m){
        for(int l = 0; l < 5; ++l){
          numera += phicoefs[l] * phicoefs[m] * sumyx(j, 5*m + l);
        }
      }
      
      for(int k = 0; k < dimx; ++k){
        //sum is over all parts of beta except beta j
        if(k != j){
          int counter = 0;
          for(int m = 0; m < 5; ++m){
            for(int l = 0; l < 5; ++l){
              if(l > m){
                numerb += beta[k] * phicoefs[l] * phicoefs[m] * (sumxx(j, k, counter) + sumxx(k, j, counter));
                counter += 1;                       
              } else if (l == m){                                                               
                numerb += beta[k] * phicoefs[l] * phicoefs[m] * sumxx(k, j, counter);
                counter += 1;       
              }
            }
          }
        }
      }
      
      numer = tau[j] * (numera - numerb) + sigmasq*mu[j];
      
      denom = 0;
      int counter = 0;
      for(int m = 0; m < 5; ++m){
        for(int l = 0; l < 5; ++l){
          if(l > m){
            denom += 2 * tau[j] * phicoefs[l] * phicoefs[m] * (sumxx(j, j, counter));
            counter += 1;                       
          } else if (l == m){                                                               
            denom += tau[j] * phicoefs[l] * phicoefs[m] * sumxx(j, j, counter);
            counter += 1;       
          }
        }
      }
      
      denom += sigmasq;
      var = tau[j] * sigmasq / denom;
      
      beta[j] = numer / denom + sqrt(var) * randn<vec>(1)[0];
    }
    
    //Sigma squared conditional
    
    double error = 0;
    double sse = 0;
    //Need to add half the sum of squared errors to betabar
    for(int i = 17856; i < T; ++i){
      error = y[i] - rowprodsum(beta, x.row(i)) - phi[0] * (y[i - 1] - rowprodsum(beta, x.row(i - 1))) - phi[1]*(y[i-48] - rowprodsum(beta, x.row(i-48))) -
        phi[2] * (y[i - 336] - rowprodsum(beta, x.row(i - 336))) - phi[3] * (y[i - 17520] - rowprodsum(beta, x.row(i - 17520)));
      sse += error * error;
    }
    betabar = betasig + 1/2 * sse;
    //betabar is gamma rate, randg takes scale 
    sigmasq = 1 / randg<vec>(1, distr_param(alphabar, 1 / betabar))[0];
    
    //Store draws
    if(r >= burn & (r - burn) % thin == 0){
      for(int j = 0; j < dimx; ++j){
        store(keepiter, j) = beta[j];
      }
      for(int j = 0; j < 4; ++j){
        store(keepiter, j + dimx) = phi[j];
      }
      store(keepiter, dimx + 4) = sigmasq;
      keepiter += 1;
    }
  }
  return store;
}