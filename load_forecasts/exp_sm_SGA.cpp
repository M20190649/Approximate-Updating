// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
#define my_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062


//Evaluate log joint density
double log_joint_dens(vec y, vec theta){
  int T = y.n_elem;
  double sigmasq = theta[2];
  vec initial = {theta[3], theta[4], theta[5], theta[6], theta[7], theta[8], theta[9]};
  int dim = 7;
  vec lags = {1, 7};
  vec alpha_full = {theta[0], theta[1], 0, 0, 0, 0, 0};
  vec x = {1, -1, -1, -1, -1, -1, -1};
  mat Trans(7, 7, fill::zeros);
  Trans(0, 0) = 1;
  Trans.row(1) = {0, -1, -1, -1, -1, -1, -1};
  for(int i = 2; i < 7; ++i){
    Trans(i, i-1) = 1;
  }
  mat x_tilde(T, dim, fill::zeros); //X_tilde_t vector for t = 1.. T
  vec y_tilde(T);  //Y_tilde_t scalar for t = 1...T
  mat b_bar(dim, T ); //b predictions
  vec alpha_full(dim, fill::zeros); //entire parameter vector including zeroes for b_t = T b_t-1 + alpha_full*et
  vec csum_lags(size, fill::zeros); //cumulative sum of lags, starting at 0, omitting last term
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

double densigammac(double x, double alpha, double beta){
  double cons = pow(beta, alpha) / tgamma(alpha); //constant
  double kernel = pow(x, -(alpha + 1)) * exp(-beta/x); //kernel of distribution
  return cons * kernel;
}

double dnormc(double x, double mu, double sigmasq){
  double cons = 1 / sqrt(2*my_PI*sigmasq);
  double kernel = exp( - 1 / (2 * sigmasq) * pow(x - mu, 2));
  return cons * kernel;
}

double logq(vec theta, vec lambda, mat par1, mat par2, mat copula, mat type){
  vec marginals(10);
  vec copulas(45); //Remove elements for each independence copula
  int type;
  for(int i = 0; i < 10; ++i){
    if(i == 0 | i == 1){
      marginals[i] = log(lambda[5*i + 4] * dnormc(theta[i], lambda[5*i], lambda[5*i + 1]) + 
        (1-lambda[5*i] + 4) * dnormc(theta[i], lambda[5*i + 2], lambda[5*i + 3]));
    } else if (i == 2) {
      marginals[i] = log(densigammac(theta[i], lambda[10], lambda[11]));
    } else {
      marginals[i] = log(dnormc(theta[i], lambda[6 + 2*i], lambda[7 + 2*i]));
    }
  }
  for(int i = 0; i < 45; ++i){
    type = type(copula(i, 0), copula(i, 1));
    if(type == 1) {
      
    } else if(type == 2) {
      
    } //Continue for all types of copula
    //May need to code copula density functions themselves
  }
  return sum(marginals) + sum(copulas);
}

double deriv_lambda(vec theta, vec lambda, mat par1, mat par2, vec y, mat copula, mat type, int j){}

double deriv_eta1(vec theta, vec lambda, mat par1, mat par2, vec y, mat copula, mat type, int j){}

double deriv_eta2(vec theta, vec lambda, mat par1, mat par2, vec y, mat copula, mat type, int j){}

//Calculate ELBO to check convergence
// [[Rcpp::export]]
double ELBO(vec y, vec lambda, mat par1, mat par2, mat theta, int reps = 10){
  double output = 0;
  for(int i = 0; i < reps; ++i){
    output += (log_joint_dens(y, theta.row(i)) - logq(theta.row(i), lambda, par1, par2))/reps;
  }
  return output;
}

// [[Rcpp::export]]
mat SGA_lambda(vec y, vec theta, vec lambda, mat par1, mat par2, double tuning, vec Gt, mat copula, mat type, int s){
  int n = lambda.n_elem;
  vec partials(n); partials.fill(0); //Store partial derivatives in calculation 
  vec pt(n); //Step size
  for(int i = 0; i < s; i++){ //for each sample in s
    for(int j = 0; j < n; j++){  //and for each parameter in lambda
      partials[j] += deriv_lambda(theta.row(i), lambda, par1, par2, y, copula, type, j)/s; //partial deriv = mean of s monte carlo estimates
    }
  }
  for(int j = 0; j < n; j++){
    Gt[j] += pow(partials[j],2); //Each GT is the sum of the squared partial diffs. 
    pt[j] = tuning * pow(Gt[j], -0.5); //transform to get pt value
    lambda[j] = lambda[j] + pt[j]*partials[j];  //new lambda is old + stepsize * partial diff
  }
  mat output(2, n);
  output.row(0) = lambda;
  output.row(1) = Gt;
  return output;
}

// [[Rcpp::export]]
cube SGA_eta1(vec y, vec theta, vec lambda, mat par1, mat par2, double tuning, mat Gt, mat copula, mat type, int s){
  int n = copula.n_rows;
  mat partials(10, 10, fill::zeros); //Store partial derivatives in calculation 
  mat pt(10, 10, fill::zeros); //Step size
  for(int i = 0; i < s; i++){ //for each sample in s
    for(int j = 0; j < n; ++j){ //and for each parameter in eta
        partials(copula(j, 0), copula(j, 1)) += 
          deriv_eta1(theta.row(i), lambda, par1, par2, y, copula, type, j)/s;
          //partial deriv = mean of s monte carlo estimates
    }
  }
  for(int j = 0; j < n; ++j){
    Gt(copula(j, 0), copula(j, 1)) += pow(partials(copula(j, 0), copula(j, 1)),2); //Each GT is the sum of the squared partial diffs. 
    pt(copula(j, 0), copula(j, 1)) = tuning * pow(Gt(copula(j, 0), copula(j, 1)), -0.5); //transform to get pt value
    par1(copula(j, 0), copula(j, 1)) = par1(copula(j, 0), copula(j, 1)) + 
    pt(copula(j, 0), copula(j, 1)) * partials(copula(j, 0), copula(j, 1));  //new lambda is old + stepsize * partial diff
  }
  cube output(10, 10, 2);
  output.slice(0) = par1;
  output.slice(1) = Gt;
  return output;
}

// [[Rcpp::export]]
cube SGA_eta2(vec y, vec theta, vec lambda, mat par1, mat par2, double tuning, mat Gt, mat copula, mat type, int s){
  int n = copula.n_rows;
  mat partials(10, 10, fill::zeros); //Store partial derivatives in calculation 
  mat pt(10, 10, fill::zeros); //Step size
  for(int i = 0; i < s; i++){ //for each sample in s
    for(int j = 0; j < n; ++j){ //and for each parameter in eta
      partials(copula(j, 0), copula(j, 1)) += 
        deriv_eta2(theta.row(i), lambda, par1, par2, y, copula, type, j)/s;
      //partial deriv = mean of s monte carlo estimates
    }
  }
  for(int j = 0; j < n; ++j){
    Gt(copula(j, 0), copula(j, 1)) += pow(partials(copula(j, 0), copula(j, 1)),2); //Each GT is the sum of the squared partial diffs. 
    pt(copula(j, 0), copula(j, 1)) = tuning * pow(Gt(copula(j, 0), copula(j, 1)), -0.5); //transform to get pt value
    par2(copula(j, 0), copula(j, 1)) = par2(copula(j, 0), copula(j, 1)) + 
      pt(copula(j, 0), copula(j, 1)) * partials(copula(j, 0), copula(j, 1));  //new lambda is old + stepsize * partial diff
  }
  cube output(10, 10, 2);
  output.slice(0) = par1;
  output.slice(1) = Gt;
  return output;
}