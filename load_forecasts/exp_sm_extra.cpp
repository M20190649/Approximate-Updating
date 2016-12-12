// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

alpha_details log_alpha_conditional(vec y, vec alpha, mat Trans, vec x, vec lags, double sigmasq, vec initial){
  int size = lags.n_elem; //Number of parameters in alpha
  int T = y.n_elem; //Length of time series
  int dim = sum(lags) - size + 1; //Sum of lags is dimension of Trans/D/states
  mat x_tilde(T, dim, fill::zeros); //X_tilde_t vector for t = 1.. T
  vec y_tilde(T);  //Y_tilde_t scalar for t = 1...T
  mat b_bar(dim, T ); //b predictions
  vec alpha_full(dim, fill::zeros); //entire parameter vector including zeroes for b_t = T b_t-1 + alpha_full*et
  vec csum_lags(size, fill::zeros); //cumulative sum of lags, starting at 0, omitting last term
  //cumsum is to indicate which elements of alpha_full are nonzero, should have indices for lt and first value of each seasonal effect
  //In comparision x typically had lt and last element of each seasonal effect (in old parameterisation)
  
  //Create D matrix
  for(int i = 0; i < size; ++i){
    if(i > 0) { //first term is already set to 0
      csum_lags[i] = csum_lags[i-1] + lags[i-1]; 
    }
    alpha_full[csum_lags[i]] = alpha[i]; //nonzero element of alpha_full at cumsum values (-1 for count starting at 0)
  }
  mat D(dim, dim); //D matrix
  D = Trans - alpha_full * x.t();
  //Calculate tilde data
  x_tilde.row(0) = x.t(); //x_tilde_1
  b_bar.col(0) = alpha_full * y[0]; //b_bar_1
  y_tilde[0] = y[0]; //y_tilde_1
  for(int t = 1; t < T; ++t){ //Loop for other values, recursion in paper
    b_bar.col(t) = D * b_bar.col(t-1) + alpha_full * y[t];
    x_tilde.row(t) = x_tilde.row(t-1) * D; 
    y_tilde[t] = (y[t] - x.t() * b_bar.col(t-1))[0]; 
  }
  //Evaluate likelihood (proportional to conditional density)
  double log_density = 0;
  for(int t = 0; t < T; ++t){
    log_density += -1/(2*sigmasq) * pow((y_tilde[t] - x_tilde.row(t) * initial)[0], 2); 
  }
  //Final computations
  mat xtxinv(dim, dim);
  xtxinv = (x_tilde.t() * x_tilde).i();
  vec beta0hat(dim);
  beta0hat = xtxinv * x_tilde.t() * y_tilde;
  double ssquared;
  ssquared = (((y_tilde - x_tilde * beta0hat).t() * (y_tilde - x_tilde * beta0hat)) / (T - size))[0];
  
  //Returning results
  alpha_details results;
  results.density = log_density;
  results.xtxinv = xtxinv;
  results.ssquared = ssquared;
  results.beta0hat = beta0hat;
  return results;
}

// [[Rcpp::export]]
mat MCMC_exp_smoothing(vec y, int rep, vec x, mat Trans, vec lags, double MHsd){
  //Inputs are data, MCMC draws, x vector, T matrix, and lag structure, eg. lags = (1, 48, 336) for level, daily and weekly seasonals
  //MHsd is a tuning parameter for MHRW algorithm
  int T = y.n_elem; //Length of time series
  int size = lags.n_elem; //Number of seasonal effects plus one for level
  int dim = sum(lags) - size + 1; //Parameterisation includes 1 term for l, then max lag - 1 for each seasonal effect.
  mat draws(rep, dim + size + 1); //size many smoothing parameters, dim many initial states, plus sigma squared
  //First columns are for smoothing parameters, then initial states, then sigma squared
  
  //Initialising various parameters to be drawn
  vec initial = randn<vec>(dim); //initial states of latent variables
  double sigmasq = 1; //variance
  vec alpha(size); alpha.fill(0.5); //smoothing parameters
  vec old_alpha(size);
  //Used in density calculations
  alpha_details alpha_result; //Marginal density result plus all the data calculations required stored in a structure
  alpha_details old_alpha_result;
  alpha_result = log_alpha_conditional(y, alpha, Trans, x, lags, sigmasq, initial);
  double IG_shape = (T - dim)/2;
  double IG_scale;
  mat var(dim, dim); //Initial states have MVN conditional, store the variance here
  mat L(dim, dim); //For the lower triangle of var
  vec eps(dim); //For standard normal variables in the MVN
  
  int reject = 0; //Count rejection rate
  bool inside_range; //Check alpha parameters inside(0, 1)
  double MH_ratio;
  double u;
  
  //Main loop
  for(int r = 0; r < rep; ++r){
    
    //Alpha drawn from its conditional via Metropolis Hastings
    old_alpha = alpha;
    old_alpha_result = alpha_result;
    //RW constrained to (0,1) - Find how truncation effects accept/reject step, probably not symmetric density
    for(int i = 0; i < size; ++i){
      inside_range = FALSE;
      while(!inside_range){
        alpha = old_alpha + MHsd * randn<vec>(size); //Random Walk Metropolis Hastings, can tune sd parameter
        inside_range = alpha[i] > 0 & alpha[i] < 1;
      }
    }
    
    alpha_result = log_alpha_conditional(y, alpha, Trans, x, lags, sigmasq, initial);
    MH_ratio = exp(alpha_result.density - old_alpha_result.density);
    u = randu<vec>(1)[0];
    if(u > MH_ratio){
      reject += 1;
      alpha = old_alpha;
      alpha_result = old_alpha_result;
    }
    
    //Sigmasquared from Inverse Gamma Density
    IG_scale = (T - size) / 2 * alpha_result.ssquared;
    sigmasq = 1 / randg<vec>(1, distr_param(IG_shape, 1/IG_scale))[0]; //Arma uses gamma(shape scale), convert IG scale to rate
    
    //Initial states from Gaussian Density, MVN via cholesky decomposition
    var = sigmasq * alpha_result.xtxinv;
    L = chol(var, "lower");
    eps = randn<vec>(dim); //Standard normal variables
    for(int i = 0; i < dim; ++i){
      initial[i] = alpha_result.beta0hat[i]; //Add means to MVN
      for(int j = 0; j <= i; ++j){
        initial[i] += L(i, j) * eps[j]; //Add sum of lower triangle + epsilon
      }
    }
    
    //Storage
    for(int i = 0; i < size; ++i){ 
      draws(r, i) = alpha[i];
    }
    draws(r, size) = sigmasq;
    for(int i = 0; i < dim; ++i){
      draws(r, size + 1 + i) = initial[i];
    }
  }
  return draws;
}

// [[Rcpp::export]]
double cond_density_only(vec y, vec alpha, mat Trans, vec x, vec lags, double sigmasq, vec initial){
  int size = lags.n_elem; //Number of parameters in alpha
  int T = y.n_elem; //Length of time series
  int dim = sum(lags) - size + 1; //Sum of lags is dimension of Trans/D/states
  mat x_tilde(T, dim, fill::zeros); //X_tilde_t vector for t = 1.. T
  vec y_tilde(T);  //Y_tilde_t scalar for t = 1...T
  mat b_bar(dim, T ); //b predictions
  vec alpha_full(dim, fill::zeros); //entire parameter vector including zeroes for b_t = T b_t-1 + alpha_full*et
  vec csum_lags(size, fill::zeros); //cumulative sum of lags, starting at 0, omitting last term
  //cumsum is to indicate which elements of alpha_full are nonzero, should have indices for lt and first value of each seasonal effect
  //In comparision x typically had lt and last element of each seasonal effect (in old parameterisation)
  
  //Create D matrix
  for(int i = 0; i < size; ++i){
    if(i > 0) { //first term is already set to 0
      csum_lags[i] = csum_lags[i-1] + lags[i-1]; 
    }
    alpha_full[csum_lags[i]] = alpha[i]; //nonzero element of alpha_full at cumsum values (-1 for count starting at 0)
  }
  mat D(dim, dim); //D matrix
  D = Trans - alpha_full * x.t();
  //Calculate tilde data
  x_tilde.row(0) = x.t(); //x_tilde_1
  b_bar.col(0) = alpha_full * y[0]; //b_bar_1
  y_tilde[0] = y[0]; //y_tilde_1
  for(int t = 1; t < T; ++t){ //Loop for other values, recursion in paper
    b_bar.col(t) = D * b_bar.col(t-1) + alpha_full * y[t];
    x_tilde.row(t) = x_tilde.row(t-1) * D; 
    y_tilde[t] = (y[t] - x.t() * b_bar.col(t-1))[0]; 
  }
  //Evaluate likelihood (proportional to conditional density)
  double log_density = 0;
  for(int t = 0; t < T; ++t){
    log_density += -1/(2*sigmasq) * pow((y_tilde[t] - x_tilde.row(t) * initial)[0], 2); 
  }
  return log_density;
}

// [[Rcpp::export]]
double density_only(vec y, vec alpha, mat Trans, vec x, vec lags){
  int size = lags.n_elem; //Number of parameters in alpha
  int T = y.n_elem; //Length of time series
  int dim = sum(lags) - size + 1; //Sum of lags is dimension of Trans/D/states
  mat x_tilde(T, dim, fill::zeros); //X_tilde_t vector for t = 1.. T
  vec y_tilde(T);  //Y_tilde_t scalar for t = 1...T
  mat b_bar(dim, T ); //b predictions
  vec alpha_full(dim, fill::zeros); //entire parameter vector including zeroes for b_t = T b_t-1 + alpha_full*et
  vec csum_lags(size, fill::zeros); //cumulative sum of lags, starting at 0, omitting last term
  //cumsum is to indicate which elements of alpha_full are nonzero, should have indices for lt and first value of each seasonal effect
  //In comparision x typically had lt and last element of each seasonal effect (in old parameterisation)
  
  //Create D matrix
  for(int i = 0; i < size; ++i){
    if(i > 0) { //first term is already set to 0
      csum_lags[i] = csum_lags[i-1] + lags[i-1]; 
    }
    alpha_full[csum_lags[i]] = alpha[i]; //nonzero element of alpha_full at cumsum values (-1 for count starting at 0)
  }
  mat D(dim, dim); //D matrix
  D = Trans - alpha_full * x.t();
  
  //Calculate tilde data
  x_tilde.row(0) = x.t(); //x_tilde_1
  b_bar.col(0) = alpha_full * y[0]; //b_bar_1
  y_tilde[0] = y[0]; //y_tilde_1
  for(int t = 1; t < T; ++t){ //Loop for other values, recursion in paper
    b_bar.col(t) = D * b_bar.col(t-1) + alpha_full * y[t];
    x_tilde.row(t) = x_tilde.row(t-1) * D; 
    y_tilde[t] = (y[t] - x.t() * b_bar.col(t-1))[0]; 
  }
  
  //Final computations
  mat xtxinv(dim, dim);
  xtxinv = (x_tilde.t() * x_tilde).i();
  vec beta0hat(dim);
  beta0hat = xtxinv * x_tilde.t() * y_tilde;
  double ssquared;
  ssquared = ((y_tilde - x_tilde * beta0hat).t() * (y_tilde - x_tilde * beta0hat) / (T - size))[0];
  double density;
  density = pow(det(xtxinv), 0.5) * pow(ssquared, (-(T - size)/2));
  return density;
}