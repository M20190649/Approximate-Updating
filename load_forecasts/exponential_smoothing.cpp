// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

mat statevector(vec y, vec alpha, vec initial, int T){
  mat states(385, T+1);
  states.row(0) = initial;
  for(int t = 1; t < T+1; ++t){
    states(0, t) = (1 - alpha[0])*states(0, t-1) - alpha[0]*states(48,t-1) - alpha[0]*states(384, t-1) + alpha[0]*y[t];
    states(1, t) = (1 - alpha[1])*states(48, t-1) - alpha[1]*states(0,t-1) - alpha[1]*states(384, t-1) + alpha[1]*y[t];
    states(49, t) = (1 - alpha[2])*states(384, t-1) - alpha[2]*states(0,t-1) - alpha[2]*states(48, t-1) + alpha[2]*y[t];
    for(int j = 2; j < 385; ++j){
      if(j != 49) {
        states(j, t) = states(j-1, t-1);
      }
    }
  }
  return states;
}

// [[Rcpp::export]]
mat m_alpha(vec alpha, vec y){
  int T = y.n_elem;
  mat D(385, 385);
  sp_mat Trans(385, 385);
  Trans(0, 0) = Trans(1, 48) = Trans(49, 384) = 1;
  for(int i = 2; i < 385; ++i){
    if(i != 49) {
      Trans(i, i-1) = 1;
    }
  }
  mat x_tilde(T, 385);
  vec y_tilde(T);
  mat b_bar(385, T);
  mat xtxinv(385, 385);
  vec x(385, fill::zeros);
  x[0] = x[48] = x[384] = 1;
  vec alpha_full(385, fill::zeros);
  alpha_full[0] = alpha[0];
  alpha_full[1] = alpha[1];
  alpha_full[49] = alpha[2];
  D = Trans - x * alpha_full.t();
  return D;
  //vec s_squared(1);
  //vec s(T);
  //vec beta_0_hat(385);
  
  //Calculate tilde data
  //x_tilde.row(0) = x.t() * D;
  //b_bar.col(0) = alpha_full * y[0];
  //y_tilde[0] = y[0];
  //for(int t = 1; t < T; ++t){
  //  b_bar.col(t) = D * b_bar.col(t-1) + alpha_full * y[t];
  //  x_tilde.row(t) = x_tilde.row(t-1) * D; 
  //  y_tilde[t] = (y[t] - x.t() * b_bar.col(t-1))[0]; 
  //}

  //xtxinv = inv(x_tilde.t() * x_tilde);
  //beta_0_hat = xtxinv * x_tilde.t() * y_tilde;
  //s = y_tilde - x_tilde * beta_0_hat;
  //s_squared = s.t() * s;
  
  //double density = (pow(det(xtxinv), 0.5) * pow(s_squared/(T-3), (T-3)/2))[0];
  //return density;
}



// [[Rcpp::export]]
mat exponential_smoothing(vec y, int rep){
  //Setup
  int T = y.n_elem;
  vec x(385, fill::zeros);
  x[0] = x[48] = x[384] = 1;
  sp_mat Trans(385, 385);
  Trans(0, 0) = Trans(1, 48) = Trans(49, 384) = 1;
  for(int i = 2; i < 49; ++i){
    Trans(i, i-1) = 1;
  }
  for(int i = 50; i < 385; ++i){
    Trans(i, i-1) = 1;
  }

  sp_mat D(385, 385);
  vec alpha(3);
  mat draws(rep, 7);
  mat x_tilde(T, 385);
  vec y_tilde(T);
  mat b_bar(385, T);
  mat xtxinv(385, 385);
  double c = 0; //replace with max m(alpha)
  
  double IG_shape = (T - 385)/2;
  double IG_scale;
  vec s_squared(1);
  vec s(T);
  vec beta_0_hat(385);
  
  //main loop
  for(int r = 0; r < rep; ++r){
    
    //Draw alpha candidate from uniform, apply accept-reject sampling
    bool flag = FALSE;
    while(!flag){
      alpha = randu<vec>(3);
      vec alpha_full(385, fill::zeros);
      alpha_full[0] = alpha[0];
      alpha_full[1] = alpha[1];
      alpha_full[49] = alpha[2];
      D = Trans - x * alpha_full.t();
      
      //Calculate tilde data
      x_tilde.row(0) = x.t() * D;
      b_bar.col(0) = alpha_full * y[0];
      y_tilde[0] = y[0];
      for(int t = 1; t < T; ++t){
        b_bar.col(t) = D * b_bar.col(t-1) + alpha_full * y[t];
        x_tilde.row(t) = x_tilde.row(t-1) * D; 
        y_tilde[t] = (y[t] - x.t() * b_bar.col(t-1))[0]; 
      }
      
      xtxinv = inv(x_tilde.t() * x_tilde);
      beta_0_hat = xtxinv * x_tilde.t() * y_tilde;
      s = y_tilde - x_tilde * beta_0_hat;
      s_squared = s.t() * s;
      
      double acceptance_ratio = (pow(det(xtxinv), 0.5) * pow(s_squared/(T-385), (T-385)/2))[0]/c;
      double u = randu<vec>(1)[0];
      if(u < acceptance_ratio){
        flag = TRUE;
      }
    }
    draws(r, 0) = alpha[0];
    draws(r, 1) = alpha[1];
    draws(r, 2) = alpha[2];

    //Draw sigma squared from IG
    IG_scale = 2/(s_squared[0]);
    
    draws(r, 3) = randg<vec>(1, distr_param(IG_shape, 1/IG_scale))[0];
    
    //Draw beta_0 from MVN
    mat L = chol(draws(r, 3) * xtxinv, "lower");
    vec eps = randn<vec>(3);
    draws(r, 4) = beta_0_hat[0] + L(0, 0) * eps[0];
    draws(r, 5) = beta_0_hat[1] + L(1, 0) * eps[0] + L(1, 1) * eps[1];
    draws(r, 6) = beta_0_hat[49];
    for(int i = 0; i < 50; ++i){
      draws(r, 6) += L(49, i) * eps[i];
    }
  } 
  return draws;
}