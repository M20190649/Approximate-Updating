// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

//Calculate ELBO to check convergence
// [[Rcpp::export]]
double ELBO(vec lambda, vec y, int reps = 10){}

//Calculate partial derivatives
double deriv(vec theta, vec lold, vec y, int j, mat Trans, int dim, int T, vec lags, int m){}

//Evaluate log joint density
double log_joint_dens(vec y, vec alpha, double sigmasq, vec initial, vec x, mat Trans, int dim, int T, vec lags){
  mat x_tilde(T, dim, fill::zeros); //X_tilde_t vector for t = 1.. T
  vec y_tilde(T);  //Y_tilde_t scalar for t = 1...T
  mat b_bar(dim, T ); //b predictions
  vec alpha_full(dim, fill::zeros); //entire parameter vector including zeroes for b_t = T b_t-1 + alpha_full*et
  vec csum_lags(size, fill::zeros); //cumulative sum of lags, starting at 0, omitting last term
  for(int i = 0; i < size; ++i){
    if(i > 0) { //first term is already set to 0
      csum_lags[i] = csum_lags[i-1] + lags[i-1]; 
    }
    alpha_full[csum_lags[i]] = alpha[i]; //nonzero element of alpha_full at cumsum values (-1 for count starting at 0)
  }
  mat D(dim, dim); //D matrix
  D = Trans - alpha_full * x.t();
  x_tilde.row(0) = x.t(); //x_tilde_1
  b_bar.col(0) = alpha_full * y[0]; //b_bar_1
  y_tilde[0] = y[0]; //y_tilde_1
  for(int t = 1; t < T; ++t){ //Loop for other values, recursion in paper
    b_bar.col(t) = D * b_bar.col(t-1) + alpha_full * y[t];
    x_tilde.row(t) = x_tilde.row(t-1) * D; 
    y_tilde[t] = (y[t] - x.t() * b_bar.col(t-1))[0]; 
  }
  double density = 0;
  for(int i = 0; i < T; ++i){
    density += (-1/(2*sigmasq) * (y_tilde[i] - x_tilde.row(i) * initial) * (y_tilde[i] - x_tilde.row(i) * initial))[0]
  }
  return density;
}

vec SGA_exp_smoothing(vec y, vec initial, int s, double eta, double threshold, int M, vec lags, mat Trans, vec x){
  int n = initial.n_elem; //Number of parameters to optimise over
  int m = lags.n_elem;
  int dim = sum(lags) - m + 1;
  vec lnew = initial; //Update paramaeters
  vec lold(n); //Store current parameters
  vec partials(n); //Store partial derivatives in calculation 
  vec Gt = zeros<vec>(n); //Step size intermediate calculation
  vec pt(n); //Step size
  double k = 0; //Number of iterations
  double LBold; //Current lower bound
  double LBnew = ELBO(lnew, 10, y); //Updated lower bound
  double diff; //Difference between the two
  mat theta(s, 10); //Simulated variables
  
  while(diff > threshold | k < 50){
    if(k > M){ //Ensure algorithm stops if it fails to converge
      // Print some failure message to console
      break; 
    }
    lold = lnew; //Push last iteration parameters to 'old'
    LBold = LBnew;
    partials.fill(0); //Sum for partial diffs starts at 0
    theta = simulate(s, lnew); //Write this simulation function
    for(int i = 0; i < s; i++){ //for each sample in s
      for(int j = 0; j < n; j++){  //and for each parameter in lambda
        partials[j] += deriv(theta.row(i), lold, y, j)/s; //partial deriv = mean of s monte carlo estimates
      }
    }
    for(int j = 0; j < n; j++){
      Gt[j] += pow(partials[j],2); //Each GT is the sum of the squared partial diffs. 
      pt[j] = eta * pow(Gt[j], -0.5); //transform to get pt value
      lnew[j] = lold[j] + pt[j]*partials[j];  //new lambda is old + stepsize * partial diff
    }
    LBnew = ELBOc(lnew, y); //calculate lowerbound with new variables
    diff = abs(LBnew - LBold); //absolute value of the difference
    k = k + 1; //count iterations so far
  } 
  //Print k and LBnew (or add to return vector)
  return lnew;
}
