// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
vec colmean(mat y){ //input matrix, return vector of column means
  int n = y.n_rows;
  int m = y.n_cols;
  vec out(m);
  for(int i = 0; i < m; ++i){  //for each column i calculate 1/m sum_j^m (x_ij) 
    out(i) = 0;
    for(int j = 0; j < n; ++j){
      out(i) += y(j,i)/n;
    }
  }
  return out;
}

// [[Rcpp::export]] //input matrix of data, retained and burn in draws, hyperparams
mat gibbschain(mat y, int M, int B, double mubarbar, double lambda, double a, double b){ 
  vec ybar = colmean(y);     
  int J = y.n_cols; //number of columns in y
  int N = y.n_rows; //number of rows in y
  
  double var; //various parameters
  vec mujbar(100);
  double muhat;
  double shape = J/2+a;
  double rate;
  double mubarmean;
  double mubarvar;

  mat draws(B+M,102); //create matrix of draws, 100 muj, mubar, tau2
  for(int i = 0; i < 102; ++i){
    draws(0,i) = 1; //start all variables at 1
  }
  
  for(int j = 1; j < B+M; ++j){ //main gibbs loop
    var = 1/(N+pow(draws(j-1, 101), -2)); //posterior variance parameter
    for(int i = 0; i < 100; ++i){
      mujbar(i) = (N*ybar(i)+draws(j-1, 100)*pow(draws(j-1, 101), -2))*var; //individual posterior means
      draws(j, i) = mujbar(i) + sqrt(var)*randn<vec>(1)[0];  //draw from muj conditional
    }
    
    muhat = 0;
    for(int i = 0; i < 100; ++i){
      muhat += draws(j,i)/J; //calculate mean of muj draws
    }
    
    mubarvar = 1/(J*pow(draws(j-1, 101), -2)+pow(lambda, -2));  //steps for mubar conditional
    mubarmean = (J*muhat*pow(draws(j-1, 101), -2)+mubarbar*pow(lambda, -2))*mubarvar;
    draws(j, 100) = mubarmean + sqrt(mubarvar)*randn<vec>(1)[0];
    
    rate = b;   //steps for tau conditional
    for(int i = 0; i < 100; ++i){
      rate += pow(draws(j,i)-draws(j,100), 2)/2;
    }
    draws(j, 101) = 1/(randg<vec>(1, distr_param(shape, 1/rate))[0]);
  }
  return draws;
}

// [[Rcpp::export]]
mat varfit(mat y, int n, double mubarbar, double lambda, double a, double b) {
  bool flag = FALSE;
  mat var(n, 105, fill::ones); //0-99 mustar, 100 lstar, 101 mubar, 102 lbar, 103 alpha, 104 beta
  int J = y.n_cols; //number of columns in y
  int N = y.n_rows; //number of rows in y
  vec ybar = colmean(y);
  double bona;
  double temp;
  
  for(int i = 1; i < n; ++i){
    var(i, 103) = J/2 + a; //alpha is constant
    bona = var(i-1, 104)/var(i, 103); //beta divided by alpha
    
    var(i, 100) = 1/(N+bona); //lambda star
    
    for(int j = 0; j < 100; ++j) { //mu star
      var(i, j) = (N*ybar[j]+var(i-1, 101)*bona)*var(i, 100);
    }
    
    var(i, 102) = 1/(J*bona+pow(lambda, -2.0)); //lambda bar
   
    temp = 0; //sum for mubar
    for(int j = 0; j < 100; ++j) {
      temp += var(i, j)*bona;
    }
    var(i, 101) = (temp + mubarbar*pow(lambda, -2.0))*var(i, 102); //mubar
    
    var(i, 104) = b;  //beta has a sum component
    for(int j = 0; j < 100; ++j) {
      var(i, 104) += (pow(var(i, j), 2) + var(i, 100) + pow(var(i, 101), 2) + var(i, 102) - 2*var(i, j)*var(i, 101))/2;
    }
  }
  return var;
}
