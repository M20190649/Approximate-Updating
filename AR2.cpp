// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]


#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions/inverse_gamma.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;
using boost::math::inverse_gamma;
using namespace boost::math;

//arma = bulk of vectors/matrices used (R package = "RcppArmadillo")
//boost = has IG quantile function (R package = "BH")

//Pi needed for gaussian density calculation
#define my_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

//Only the SGA function needs to be called in R. Everything else is called from SGA.

//Takes the inverse cdf transform of a uniform(0,1) random variable
// [[Rcpp::export]]
double qigammac(double x, double alpha, double beta){
  inverse_gamma IG(alpha, beta); //creates an IG object
  double out = quantile(IG, x); //quantile function is the inv. CDF
  return out;
}

//Simulates n many trios of N(0,I) and U(0,1) variables
// [[Rcpp::export]]
mat simc(int n){
  mat out(n, 3);
  for(int i = 0; i < n; i++) {
    out(i, 0) = randn<vec>(1)[0];
    out(i, 1) = randn<vec>(1)[0];
    out(i, 2) = randu<vec>(1)[0];
  }
  return(out);
}

//Evaluates the log joint density. Hyperparameters are included with default values
// [[Rcpp::export]]
double logjointc(vec y, double phi1, double phi2, double sig2, double phi1b = 0, double phi2b = 0, double phi1l = 10, double phi2l = 10, double a = 1, double b = 1){
  int T = y.size();
  double loglike = -(T-2)*log(sig2)/2; //loglikelihood first term
  for(int i = 2; i < T; i++){
    loglike = loglike - pow(y[i] - phi1*y[i-1] - phi2*y[i-2], 2)/(2*sig2); //the sum of y - mean(y) from t = 3 to T
  }
  double p1 = pow(phi1 - phi1b, 2)/phi1l + pow(phi2 - phi2b, 2)/phi2l; //phi priors
  double p2 = -(a+1)*log(sig2) - b/sig2; //sigma squared prior
  return loglike - p1/2 + p2;
} 

//Returns the log density of a bivariate normal distribution
// [[Rcpp::export]]
double logbvnc(vec x, vec mu, mat Sigma){ 
  double rho = Sigma(0, 1) / pow(Sigma(0,0) * Sigma(1,1), 0.5); //calculates correlation from Sigma matrix
  double cons = - log(2*my_PI*pow(Sigma(0,0) * Sigma(1,1), 0.5)*sqrt(1-pow(rho, 2))); //bvn constant
  double z = pow(x[0] - mu[0], 2)/Sigma(0,0) + pow(x[1] - mu[1], 2)/Sigma(1,1) - 2*rho*(x[0]-mu[0])*(x[1] - mu[1])/ pow(Sigma(0,0) * Sigma(1,1), 0.5); //inside exp term
  double out = cons - z/(2*(1-rho*rho)); //full log distribution
  return out;
}

//Returns the density of an inverse gamma distribution
// [[Rcpp::export]]
double densigammac(double x, double alpha, double beta){
  double cons = pow(beta, alpha) / tgamma(alpha); //constant
  double kernel = pow(x, -(alpha + 1)) * exp(-beta/x); //kernel of distribution
  return cons * kernel;
}

//Evaluates the log of the entire q distribution
//No longer needed as the score has expectation zero - still part of ELBO
// [[Rcpp::export]]
double logqc(vec theta, vec lambda){
  mat L(2,2); //steps to create lower triangle from parameter vector 
  L(0, 0) = lambda[2];
  L(1, 0) = lambda[3];
  L(1, 1) = lambda[4];
  mat Lt = L.t(); //transpose of lower triangle
  mat Sigma = L * Lt; //Variance matrix
  vec x1 = {theta[0], theta[1]}; //Only pass the first two elements of theta into bvn
  vec mu = {lambda[0], lambda[1]}; //Extract mean from the parameter vector
  double bvn = logbvnc(x1, mu, Sigma); //calculate log bvn density
  double ig = log(densigammac(theta[2], lambda[5], lambda[6])); //calculates then takes log of IG density
  return bvn + ig;
}

//Calculates numerical partial derivative of ELBO with respect to parameter j
// [[Rcpp::export]]
double allderivc(rowvec theta, vec lambda, vec y, double param){
  double h = 0.000001; //set to what you want, lower is better
  vec lambda2 = lambda; //create x+h vector
  lambda2[param] = lambda2[param] + h;
  
  vec params1(3); //actual parameters are a transform of simulated
  params1[0] = lambda[2] * theta[0] + lambda[0]; //transform standard normal to correlated bivariate, param =  mu + L * theta 
  params1[1] = lambda[3] * theta[0] + lambda[4] * theta[1] + lambda[1];
  params1[2] = qigammac(theta[2], lambda[5], lambda[6]); //inverse transform of uniform to make an IG
  
  vec params2(3); //repeated for x+h
  params2[0] = lambda2[2] * theta[0] + lambda2[0];
  params2[1] = lambda2[3] * theta[0] + lambda2[4] * theta[1] + lambda2[1];
  params2[2] = qigammac(theta[2], lambda2[5], lambda2[6]);
  
  double p1 = logjointc(y, params1[0], params1[1], params1[2]); //f(x)
  double p2 = logjointc(y, params2[0], params2[1], params2[2]); //f(x+h)
  return (p2 - p1)/h; //deriv approximately equals f(x+h) - f(x) / h
}

//evaluate the ELBO for a given set of parameters as a monte carlo estimate
//found to have low variance so a small n is fine
// [[Rcpp::export]]
double ELBOc(vec lambda, int n, vec y) {
  vec theta(3); //transformed variables
  vec eps1(2); //bivariate normal
  double out = 0; //output
  double eps2; //uniform
  for(int i = 0; i < n ; i++){ //for each i, simulate a trio of variables, transform, calculate logjoint - logq, then average over n
    eps1 = randn<vec>(2); //bivariate standard normal
    eps2 = randu<vec>(1)[0]; //uniform
    theta[0] = lambda[2]*eps1[0] + lambda[0]; //this is mu + L * epsilon 
    theta[1] = lambda[3]*eps1[0] + lambda[4]*eps1[1] + lambda[1];
    theta[2] = qigammac(eps2, lambda[5], lambda[6]);
    out += (logjointc(y, theta[0], theta[1], theta[2]) - logqc(theta, lambda))/n; //calculate mean easily by summing terms/n 
  }
  return out;
}

//put everything together in the SGA algorithm
//arguments: y = data, initial = starting lambda, s = size of monte carlo simulations, seems to work fine with s = 1 
//threshold = diff in ELBO level, try <= 0.0001, M = maximum number of iterations, 10,000 seems to work, adjust with lower thresholds, 
//eta = stepsize tuning parameter, 0.1 works fine
// [[Rcpp::export]]
vec SGA(vec y, vec initial, int s, double threshold, int M, double eta){
  
  //initial setup
  vec lnew = initial; //we have the old lambda from last iteration and the new lambda created from this iteration
  vec lold(7);
  vec partials(7); //store partial diffs
  vec Gt = zeros<vec>(7); //for the step size calculation
  vec pt(7); //final step size
  double k = 0; //counts iterations
  double LBold; //old lambdas ELBO value
  double LBnew = ELBOc(lnew, 10, y); //new lambdas ELBO value
  double diff = 1; //difference between two ELBO's, just needs to start above threshold
  mat theta(s, 3); //store simulated values for monte carlo
  
  while(diff > threshold | k < 50){ //repeat until convergence, make sure you get a few iterations in to start off as it sometimes stopped after 1
    if(k > M){ //stops it from going on forever
      break;
    }
    lold = lnew; //last set of lambdas is now saved as old
    partials = zeros<vec>(7); //reset partial derivative sums to zero
    theta = simc(s);  //simulate s trios of variables
    for(int i = 0; i < s; i++){ //for each sample in s
      for(int j = 0; j < 7; j++){  //and for each parameter in lambda
        partials[j] += allderivc(theta.row(i), lold, y, j)/s; //partial deriv = mean of s monte carlo estimates
      }
    }
    for(int j = 0; j < 7; j++){
      Gt[j] += pow(partials[j],2); //Each GT is the sum of the squared partial diffs. 
      pt[j] = eta * pow(Gt[j], -0.5); //transform to get pt value
      lnew[j] = lold[j] + pt[j]*partials[j];  //new lambda is old + stepsize * partial diff
    }
    LBold = LBnew; //last lowerbound is saved as old to match change of lambda above
    LBnew = ELBOc(lnew, 10, y); //calculate lowerbound with new variables
    diff = abs(LBnew - LBold); //absolute value of the difference
    k = k + 1; //count iterations so far
  } //end of while loop
  
  vec out = {lnew[0], lnew[1], lnew[2], lnew[3], lnew[4], lnew[5], lnew[6], k, LBnew}; //output final set of parameters, number of iterations and max LB
  return out;
}



