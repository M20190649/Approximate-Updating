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
#define my_PI 3.14159265358979324

//Only the SGA function needs to be called in R. Everything else is called from SGA.

//Takes the Inv.Gamma inverse cdf transform of a uniform(0,1) random variable
// [[Rcpp::export]]
double qigammac(double x, double alpha, double beta){
  inverse_gamma IG(alpha, beta); //creates an IG object from BH
  double out = quantile(IG, x); //quantile function is the inv. CDF
  return out;
}

//Simulates n many trios of two N(0,1) and one U(0,1) variables
// [[Rcpp::export]]
mat simc(int n){
  mat out(n, 3);
  for(int i = 0; i < n; i++) {
    out(i, 0) = randn<vec>(1)[0];
    out(i, 1) = randn<vec>(1)[0];
    out(i, 2) = randu<vec>(1)[0]; //uniform
  }
  return(out);
}

//Evaluates the log joint density. Hyperparameters are included with default values
// [[Rcpp::export]]
double logjointc(vec y, double phi1, double phi2, double sig2, double phi1b = 0, double phi2b = 0, double phi1l = 3, double phi2l = 3, double a = 1, double b = 1){
  int T = y.size();
  double loglike = -(T-2)*log(sig2)/2; //loglikelihood first term
  for(int i = 2; i < T; i++){
    loglike -=  pow(y[i] - phi1*y[i-1] - phi2*y[i-2], 2)/(2*sig2); //the sum of y - mean(y) from t = 3 to T
  }
  double p1 = -pow(phi1 - phi1b, 2)/(2*phi1l) - pow(phi2 - phi2b, 2)/(2*phi2l); //phi priors
  double p2 = -(a+1)*log(sig2) - b/sig2; //sigma squared prior
  return loglike + p1 + p2;
} 

//Returns the log density of a bivariate normal distribution
// [[Rcpp::export]]
double logbvnc(vec x, vec mu, mat Sigma){ 
  double rho = Sigma(0, 1) / pow(Sigma(0,0) * Sigma(1,1), 0.5); //calculates correlation from Sigma matrix
  double cons = - log(2*my_PI*pow(Sigma(0,0) * Sigma(1,1), 0.5)*sqrt(1-pow(rho, 2))); //bvn constant
  double z = pow(x[0] - mu[0], 2)/Sigma(0,0) + pow(x[1] - mu[1], 2)/Sigma(1,1) - 2*rho*(x[0]-mu[0])*(x[1] - mu[1])/ pow(Sigma(0,0) * Sigma(1,1), 0.5); //inside exp term
  double out = cons - z/(2*(1-pow(rho, 2))); //full log distribution
  return out;
}

//Returns the density of an inverse gamma distribution (IG(shape, scale) = 1/gamma(shape, rate))
// [[Rcpp::export]]
double logdigammac(double x, double alpha, double beta){
  double cons = alpha*log(beta) - log(tgamma(alpha)); //constant, tgamma is the cpp gamma function
  double kernel = -(alpha + 1)*log(x) - beta/x; //kernel of distribution
  return cons + kernel;
}

//Evaluates the log of the entire q distribution
// [[Rcpp::export]]
double logqc(vec theta, vec lambda){
  mat L(2,2); //steps to create lower triangle from parameter vector 
  L(0, 0) = lambda[2];
  L(1, 0) = lambda[3];
  L(1, 1) = lambda[4];
  mat Sigma = L * L.t(); //Variance matrix
  vec x1 = {theta[0], theta[1]}; //The part of theta that comes from the BVN distribution
  vec mu = {lambda[0], lambda[1]}; //Mean vector of the BVN
  double bvn = logbvnc(x1, mu, Sigma); //calculate log bvn density
  double ig = logdigammac(theta[2], lambda[5], lambda[6]); //calculates then takes log of IG density
  return bvn + ig;
}

// Takes the derivative of the log joint density with respect to either phi1 or phi2 (denoted by whichphi)
// [[Rcpp::export]]
double phideriv(vec y, double phi1, double phi2, double sigmasq, int whichphi, double lambda = 3, double phibar = 0){
  int T = y.size();
  double prior;
  if(whichphi == 1) { //derivative of log prior distribution
    prior = - 1 / lambda * (phi1 - phibar);
  } else {
    prior = - 1 / lambda * (phi2 - phibar);
  }
  double likelihood = 0;
  for(int t = 2; t < T; ++t){ //derivative of log likelihood as a sum over t
    if(whichphi == 1) {
      likelihood -= phi1*pow(y[t-1], 2) - y[t]*y[t-1] + phi2*y[t-1]*y[t-2];
    } else {
      likelihood -= phi2*pow(y[t-2], 2) - y[t]*y[t-2] + phi1*y[t-1]*y[t-2];
    }
  }
  likelihood = likelihood / sigmasq;
  double out = prior + likelihood; 
  return out;
}

//Derivative of log joint with respect to sigma squared
// [[Rcpp::export]]
double sigderiv(vec y, double phi1, double phi2, double sigmasq, double a = 1, double b = 1){
  double T = y.size();
  double part1 = -(T + a - 1)/sigmasq; //components divided by sigma squared (originally out front of distribution)
  double part2 = b; //components divided by sigma^4 (originally in the exponents)
  for(int t = 2; t < T; ++t){
    part2 += 1/2 * pow(y[t] - phi1*y[t-1] - phi2*y[t-2], 2);
  }
  part2 = part2 / pow(sigmasq, 2);
  return part1 + part2;
}

//Calculates partial derivative of ELBO with respect to parameter j
// [[Rcpp::export]]
double allderivc(rowvec epsilon, vec lambda, vec y, double param){
  
  vec theta(3); //actual parameters are a transform of simulated epsilons
  theta[0] = lambda[2]*epsilon[0] + lambda[0]; //transform standard normal to correlated bivariate, theta =  mu + L * epsilon
  theta[1] = lambda[3]*epsilon[0] + lambda[4]*epsilon[1] + lambda[1];
  theta[2] = qigammac(epsilon[2], lambda[5], lambda[6]); //inverse transform of uniform to make an IG
  double deriv1;
  //Phi derivatives are easily obtained via the chain rule
  if(param == 0){
    deriv1 = phideriv(y, theta[0], theta[1], theta[2], 1);
  } else if(param == 1) {
    deriv1 = phideriv(y, theta[0], theta[1], theta[2], 2);
  } else if(param == 2) {
    deriv1 = epsilon[0] * phideriv(y, theta[0], theta[1], theta[2], 1);
  } else if(param == 3) {
    deriv1 = epsilon[0] * phideriv(y, theta[0], theta[1], theta[2], 2);
  } else if(param == 4) {
    deriv1 = epsilon[1] * phideriv(y, theta[0], theta[1], theta[2], 2);
  }
  //dsigmasq/dlambda is unknown as qigamma is complex. Do this part numerically but the d(log(p(theta, y)))/dsigmasq analytically.
  
  double h = 0.000001; //set to what you want, lower is better
  vec lambda2 = lambda; //x+h in f(x+h)
  lambda2[param] = lambda2[param] + h; //Works for either sigmasq parameter without if functions
  double sigmasq2;
  sigmasq2 = qigammac(epsilon[2], lambda2[5], lambda2[6]); //f in f(x+h)
  double dSigmasqDLambda = (sigmasq2 - theta[2])/h; //(f(x+h)-f(x))/h
  if(param == 5 | param == 6){
    deriv1 = dSigmasqDLambda *  sigderiv(y, theta[0], theta[1], theta[2]);
  }
  vec theta2(3);
  theta2[0] = lambda2[2]*epsilon[0] + lambda2[0]; 
  theta2[1] = lambda2[3]*epsilon[0] + lambda2[4]*epsilon[1] + lambda2[1];
  theta2[2] = qigammac(epsilon[2], lambda2[5], lambda2[6]);
  double deriv2 = (logqc(theta2, lambda2)-logqc(theta, lambda))/h;
  return deriv1 - deriv2;
}

//evaluate the ELBO for a given set of parameters as a monte carlo estimate
//need a large n to be a consistent enough estimator to work with our low threshold
// [[Rcpp::export]]
double ELBOc(vec lambda, vec y, int n = 1000) {
  mat epsilon = simc(n); //to transform to theta
  vec theta(3); //store transformed variables
  vec out(n); //output
  for(int i = 0; i < n ; ++i){ //for each i, transform a row of epsilon to theta, calculate logjoint - logq, then average over n
    theta[0] = lambda[2]*epsilon(i,0) + lambda[0]; //this is mu + L * epsilon 
    theta[1] = lambda[3]*epsilon(i,0) + lambda[4]*epsilon(i,1) + lambda[1];
    theta[2] = qigammac(epsilon(i,2), lambda[5], lambda[6]); //inverse transform
    out[i] = (logjointc(y, theta[0], theta[1], theta[2]) - logqc(theta, lambda)); //calculate by summing terms/n 
  }
  return mean(out);
}


//put everything together in the SGA algorithm
//arguments: y = data, lambda = starting parameters, s = size of monte carlo simulations, seems to work fine with s = 1 
//threshold = difference in ELBO level required for convergence, try <= 0.0001,
//M = maximum number of iterations, 10,000 seems to work, adjust with lower thresholds, 
//eta = stepsize tuning parameter, 0.1 works fine, don't want to restrict movement too much
// [[Rcpp::export]]
vec SGA(vec y, vec lambda, int s, double threshold, int M, double eta){
  //initial setup
  vec partials(7); //store partial derivatives
  vec Gt = zeros<vec>(7); //for the adagrad step size calculation 
  vec pt(7); //final step size
  double k = 0; //counts iterations
  double LBold; //old lambdas ELBO value
  double LBnew = ELBOc(lambda, y); //new lambdas ELBO value
  double diff = threshold + 1; //difference between two ELBO's, just needs to start above threshold
  mat epsilon(s, 3); //store simulated values for monte carlo
  
  while(diff > threshold | k < 50){ //repeat until convergence, make sure you get a few iterations in to start off as it sometimes got stuck after 1
    if(k >= M){ //stops it from going on forever
      break;
    }
    partials = zeros<vec>(7); //reset partial derivative sums to zero
    epsilon = simc(s);  //simulate s trios of variables
    for(int i = 0; i < s; ++i){ //for each sample in s
      for(int j = 0; j < 7; ++j){  //and for each parameter in lambda
        partials[j] += allderivc(epsilon.row(i), lambda, y, j)/s; //partial deriv = mean of s monte carlo estimates
      }
    }
    for(int j = 0; j < 7; ++j){
      Gt[j] += pow(partials[j],2); //Each GT is the sum of the squared partial diffs. 
      pt[j] = eta * pow(Gt[j], -0.5); //transform to get pt value
      lambda[j] += pt[j] * partials[j];  //new lambda is old + stepsize * partial diff
    }
    LBold = LBnew; //last lowerbound (ELBO) is saved as old
    LBnew = ELBOc(lambda, y); //calculate lowerbound with new variables
    diff = abs(LBnew - LBold); //absolute value of the difference
    k = k + 1; //count iterations so far
  } //end of while loop
  
  vec out = {lambda[0], lambda[1], lambda[2], lambda[3], lambda[4], lambda[5], lambda[6], k, LBnew}; //output final set of parameters, number of iterations and maximised LB
  return out;
}

// [[Rcpp::export]]
vec ELBOtest(vec lambda, vec y, int n = 100) {
  mat epsilon = simc(n);
  vec theta(3); //transformed variables
  vec out(n); //output
  for(int i = 0; i < n ; ++i){ //for each i, simulate a trio of variables, transform, calculate logjoint - logq, then average over n
    theta[0] = lambda[2]*epsilon(i,0) + lambda[0]; //this is mu + L * epsilon 
    theta[1] = lambda[3]*epsilon(i,0) + lambda[4]*epsilon(i,1) + lambda[1];
    theta[2] = qigammac(epsilon(i,2), lambda[5], lambda[6]); //inverse transform
    out[i] =  (logjointc(y, theta[0], theta[1], theta[2]) - logqc(theta, lambda)); //calculate by summing terms/n 
  }
  return out;
}

