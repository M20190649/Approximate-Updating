// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

struct alpha_details {
  double density; //marginal density
  mat xtxinv; //x transpose x inverse
  double ssquared; //data variance
  vec beta0hat; //OLS estimator 
};

//Inverse CDF method to draw from alpha marginal distribution. Approach marginalises out one variable, draws via inverse CDF, then conditions on that draw to draw second via inverse CDF
//Inputs: Precalculated marginal density for first parameter, support vector density is calculated over, bivariate density matrix
vec draw_alpha(vec XCDF, vec support, cube densityarray){
  vec out(2);
  double u = randu<vec>(1)[0];
  //I want to find the index of the first component of XCDF that is greater than u
  bool flag = FALSE;
  int i = 0;
  while(!flag){      
    flag = XCDF[i] > u;
    if(!flag){
      i += 1;
    }
  }
  out[0] = support[i];
  //Choose the column of the bivariate density that matches the first draw as the conditional density
  mat densitymat = densityarray.slice(1);
  vec conditionaly = densitymat.col(i);
  double normalisation = sum(conditionaly);
  conditionaly = conditionaly / normalisation;
  //Calculate the conditional CDF as a cumulative sum
  int N = support.n_elem;
  vec YCDF(N);
  YCDF[0] = conditionaly[0];
  for(int j = 1; j < N; ++j){
    YCDF[j] = YCDF[j-1] + conditionaly[j];
  }
  //Draw from the CDF via the same method used for X
  double v = randu<vec>(1)[0];
  flag = FALSE;
  i = 0;
  while(!flag){
    flag = YCDF[i] > v;
    if(!flag){
      i += 1;
    }
  }
  out[1] = support[i];
  return out;
}

//Three dimensional version of above. Input XCDF, draw, condition to make p(yz|x), marginalise into p(y | x), draw, condition to make p(z|yx), draw.
vec draw_alpha_3D(vec XCDF, vec support, cube densityarray){
  vec out(3);
  int N = support.n_elem;
  double u = randu<vec>(1)[0];
  //I want to find the index of the first component of XCDF that is greater than u
  bool flag = FALSE;
  int i = 0;
  while(!flag){      
    flag = XCDF[i] > u;
    if(!flag){
      i += 1;
    }
  }
  out[0] = support[i];
  //Condition on this value to make p(yz|x)
  mat densitymat = densityarray.slice(i);
  double nsum = 0;
  for(int k = 0; k < N; ++k){
    for(int j = 0; j < N; ++j){
      nsum += densitymat(k, j);
    }
  }
  densitymat = densitymat / nsum;
  //Marginalise over z to make p(y | x)
  vec ycond(N, fill::zeros);
  for(int k = 0; k < N; ++k){
    for(int j = 0; j < N; ++j){
      ycond[k] += densitymat(k, j);
    }
  }
  //Create the CDF as a cumulative sum
  vec YCDF(N);
  YCDF[0] = ycond[0];
  for(int j = 1; j < N; ++j){
    YCDF[j] = YCDF[j-1] + ycond[j];
  }
  //Draw via inverse CDF
  double v = randu<vec>(1)[0];
  flag = FALSE;
  i = 0;
  while(!flag){
    flag = YCDF[i] > v;
    if(!flag){
      i += 1;
    }
  }
  out[1] = support[i];
  //Condition on this to create p(z | xy)
  vec zcond = densitymat.row(i);
  nsum = sum(zcond);
  zcond = zcond / nsum;
  //Calculate the conditional CDF 
  vec ZCDF(N);
  ZCDF[0] = zcond[0];
  for(int j = 1; j < N; ++j){
    ZCDF[j] = ZCDF[j-1] + zcond[j];
  }
  //Draw from the CDF
  double w = randu<vec>(1)[0];
  flag = FALSE;
  i = 0;
  while(!flag){
    flag = ZCDF[i] > w;
    if(!flag){
      i += 1;
    }
  }
  out[2] = support[i];
  return out;
}


alpha_details alpha_marginal(vec y, vec alpha, mat Trans, vec x, vec lags){
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
  beta0hat = xtxinv * x_tilde.t() * y;
  double ssquared;
  ssquared = ((y_tilde - x_tilde * beta0hat).t() * (y_tilde - x_tilde * beta0hat) / (T - size))[0];
  double density;
  density = pow(det(xtxinv), 0.5) * pow(ssquared, (-(T - size)/2));
  
  //Returning results
  alpha_details results;
  results.density = density;
  results.xtxinv = xtxinv;
  results.ssquared = ssquared;
  results.beta0hat = beta0hat;
  return results;
}


// [[Rcpp::export]]
mat general_exp_sm_MCMC(vec y, int rep, vec x, mat Trans, vec lags, vec XCDF, vec support, cube densityarray){
  //Inputs are data, MCMC draws, x vector, T matrix, and lag structure, eg. lags = (1, 48, 336) for level, daily and weekly seasonals
  int T = y.n_elem; //Length of time series
  int size = lags.n_elem; //Number of seasonal effects plus one for level
  int dim = sum(lags) - size + 1; //Parameterisation includes 1 term for l, then max lag - 1 for each seasonal effect.
  mat draws(rep, dim + size + 1); //size many smoothing parameters, dim many initial states, plus sigma squared
  //First columns are for smoothing parameters, then initial states, then sigma squared
  
  //Initialising various parameters to be drawn
  vec initial(dim); //initial states of latent variables
  double sigmasq; //variance
  vec alpha(size); //smoothing parameters
  
  //Used in density calculations
  alpha_details alpha_result; //Marginal density result plus all the data calculations required stored in a structure
  double IG_shape = (T - dim)/2;
  double IG_scale;
  mat var(dim, dim); //Initial states have MVN conditional, store the variance here
  mat L(dim, dim); //For the lower triangle of var
  vec eps(dim); //For standard normal variables in the MVN
  
  //Main loop
  for(int r = 0; r < rep; ++r){

    //Alpha drawn from its marginal provided in paper
    if(size == 2){
      alpha = draw_alpha(XCDF, support, densityarray);
    } else if(size == 3){
      alpha = draw_alpha_3D(XCDF, support, densityarray); 
    }
   
    alpha_result = alpha_marginal(y, alpha, Trans, x, lags);
    
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
  beta0hat = xtxinv * x_tilde.t() * y;
  double ssquared;
  ssquared = ((y_tilde - x_tilde * beta0hat).t() * (y_tilde - x_tilde * beta0hat) / (T - size))[0];
  double density;
  density = pow(det(xtxinv), 0.5) * pow(ssquared, (-(T - size)/2));
  return density;
}


