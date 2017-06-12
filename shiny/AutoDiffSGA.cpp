// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <stan/math.hpp>
#include <vector>
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// This structure is initialised with the data y, the random noise epsilon and an indicator for a meanfield approximation to the latent states
// Once initialised, it acts as a functor (?), with the operator() allowing it to be evaluated as a function that takes an argument lambda
// It then evaluates the ELBO at that point lambda.
// The use of T types allows the variables to be of type stan::math::var which can be differentiated automatically
struct ELBO {
  const Matrix<double, Dynamic, 1> y;
  const Matrix<double, Dynamic, 1> epsilon;
  const bool meanfield;
  const int dim;
  // The constructor takes the arguments: data, epsilon, meanfield
  ELBO(const MatrixXd& yin, const MatrixXd& epsin, const bool& mf, const int& di) :
    y(yin), epsilon(epsin), meanfield(mf), dim(di) {}
  template <typename T>
  // The functor takes the argument lambda
  T operator()(const Matrix<T, Dynamic, 1>& lambda) 
    const {
    // Require these to be called explicitly to allow them to take stan::math::var as an argument instead of double
    // fabs is a version of abs that consistently returns a double (abs int = int) so is preferred apparently.
    using std::pow; using std::log; using std::exp; using std::sqrt; using std::fabs;
    // AR1 based models start with y0, so Obs = length(y) - 1. (T is used as the variable type so calling it Obs)
    int Obs = y.size() - 1;
    // Various parameters that will be used in more than one case
    T sigmaSq;
    T phi;
    T gamma;
    T plogDens;
    T qlogDens;
    // The value of dim (and number of random variables) defines the model to be evaluted. dim = 2 and 3 correspond to AR1 models,
    // dim = T+5 corresponds to the state space model
    // AR 1 model without a mean
    if(dim == 2){
      // Transform epsilon
      sigmaSq = exp(lambda(0) + lambda(2)*epsilon(0));
      phi = lambda(1) + lambda(4)*epsilon(0) + lambda(5)*epsilon(1);
      // p log density as a function of theta and y
      // If I enforce some kind of stationary restriction on phi it tends to sample from outside of the stationary region at some point
      // This will return an NA for log(1-phi^2), which then NA's everything.
      // Truncated normal sampling of phi is required to correct this
      // I could just keep re-sampling epsilon until -1 < phi < 1, but not implmented yet.
      plogDens = -0.5 * (Obs+5) * log(sigmaSq)  -  1.0 / sigmaSq  -  pow(y(0) - gamma, 2) / (2*sigmaSq);
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(y(t) - phi*y(t-1), 2) / (2*sigmaSq);
      }
      // q log density as a function of epsilon and lambda
      qlogDens = lambda(0)  +  lambda(2) * epsilon(0)  +  log(fabs(lambda(2)))  +  log(fabs(lambda(5)));
      qlogDens -= 0.5 * (pow(epsilon(0), 2) + pow(epsilon(1), 2));
    // AR 1 model with a mean 
    } else if(dim==3){
      // Transform
      sigmaSq = exp(lambda(0) + lambda(3)*epsilon(0));
      phi = lambda(1) + lambda(6)*epsilon(0) + lambda(7)*epsilon(1);
      gamma = lambda(2) + lambda(9)*epsilon(0) + lambda(10)*epsilon(1) + lambda(11)*epsilon(2);
      // p log density
      plogDens = -0.5 * (Obs+5) * log(sigmaSq)  -  1.0 / sigmaSq;  -   pow(gamma, 2) / 200 -  pow(y(0) - gamma, 2) / (2*sigmaSq);
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(y(t) - gamma - phi*(y(t-1)-gamma), 2) / (2*sigmaSq);
      }
      // q log density
      qlogDens = lambda(0)  +  lambda(3) * epsilon(0);
      for(int i = 0; i < 3; ++i){
        qlogDens += log(fabs(lambda((i+1)*dim+i)))  -  0.5 * pow(epsilon(i), 2);
      }
    // Basic Dynamic Linear model
    } else if(dim==Obs+6){
      Obs += 1; // In this case there is no y0, so Obs should be the length of y.
      // Transform variables
      // lambda 0:(dim-1) corresponds to the mean vector
      // U(i, j) is pushed to lambda(dim*i + j) when the parameters are converted to a vector.
      // Slightly different formula for this element in the loop as i and j (and lambda indices) start at 0, not 1.
      Matrix<T, Dynamic, 1> transformed(dim);
      for(int i = 0; i < dim; ++i){
        transformed(i) = lambda(i);
        if(meanfield & (i > 3)){ 
          // If meanfield is true, the latent states have a diagonal covariance.
          // We still allow dependencies in the theta part of the matrix, which happens if i <= 3.
          // So q(theta, x_{0:T}) = q(theta) \prod_{t=0}^T q(x_t) with a four dimensional theta
          transformed(i) += lambda((i+1)*dim+i) * epsilon(i);
        } else {
          // dependencies allows for theta (always) and x (if meanfield is false)
          for(int j = 0; j <= i; ++j){
            transformed(i) += lambda((i+1)*dim+j) * epsilon(j);
          }
        }
      }
      T sigmaSqY =  exp(transformed(0));
      T sigmaSqX =  exp(transformed(1));
      phi = transformed(2);
      gamma = transformed(3);
      Matrix<T, Dynamic, 1> x = transformed.bottomRows(Obs+1);
      // p log density
      plogDens = -0.5 * (Obs+4) * log(sigmaSqY)  -  0.5 * (Obs+5) * log(sigmaSqX)  -  1.0 / sigmaSqY  -
        1.0 / sigmaSqX;  -   pow(gamma, 2) / 200 -  pow(x(0) - gamma, 2) / (2*sigmaSqX);
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(x(t) - gamma - phi*(x(t-1)-gamma), 2) / (2*sigmaSqX);
        plogDens -= pow(y(t-1) - x(t), 2) / (2*sigmaSqY);
      }
      // q log density
      // Jacobian has log(sigmaSqX) + log(sigmaSqY) + sum_i log(U_ii)
      qlogDens = transformed(0) + transformed(1);
      for(int i = 0; i < 3; ++i){
        // logq(epsilon) reduces to sum_i -0.5*eps(i)^2
        qlogDens += log(fabs(lambda((i+1)*dim+i)))  -  0.5 * pow(epsilon(i), 2);
      }
      // Stochastic Volatility Model
    } else if(dim==Obs+5) {
      Obs += 1;
      Matrix<T, Dynamic, 1> transformed(dim);
      for(int i = 0; i < dim; ++i){
        transformed(i) = lambda(i);
        if(meanfield & (i > 2)){ 
          transformed(i) += lambda((i+1)*dim+i) * epsilon(i);
        } else {
            for(int j = 0; j <= i; ++j){
            transformed(i) += lambda((i+1)*dim+j) * epsilon(j);
          }
        }
      }
      sigmaSq = exp(transformed(0));
      phi = transformed(1);
      gamma = transformed(2);
      Matrix<T, Dynamic, 1> x = transformed.bottomRows(Obs+1);
      // p log density
      plogDens = -0.5 * (Obs+5) * log(sigmaSq)  -  1.0 / sigmaSq;  -   pow(gamma, 2) / 200 -  pow(x(0) - gamma, 2) / (2*sigmaSq);
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(x(t) - gamma - phi*(x(t-1)-gamma), 2) / (2*sigmaSq);
        plogDens += 0.5 * (y(t-1) - x(t) - exp(y(t-1)-x(t)));
      }
      // q log density
      qlogDens = transformed(0);
      for(int i = 0; i < 3; ++i){
        qlogDens += log(fabs(lambda((i+1)*dim+i)))  -  0.5 * pow(epsilon(i), 2);
      }
    }
    // Each estimate of the ELBO = log(p(y, theta, x)) - log(q(theta, x))
    // These are averaged over epsilon draws in AD_SGA to give \hat{ELBO} \approx E_q (log(p(y, theta, x)) - log(q(theta, x)))
    return plogDens - qlogDens;
  }
};

// [[Rcpp::export]]
Rcpp::List AD_SGA(const Rcpp::NumericMatrix yIn, Rcpp::NumericMatrix lambdaIn, int dim,
                   int M, int maxIter, bool meanfield = false, 
                   double threshold=0.01, double alpha=0.1, double beta1=0.9, double beta2=0.999){
  // Convert to Eigen Format
  Map<MatrixXd> y(Rcpp::as<Map<MatrixXd> >(yIn));
  Map<MatrixXd> lambdaMat(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  int rows = lambdaMat.rows();
  int cols = lambdaMat.cols();
  
  double e = pow(10, -8);
  // Initialise new matrices
  // lambda must be treated as a vector
  // Eigen is smart, this is a new way to look at an existing object, not creating a new object
  Map<MatrixXd> lambda(lambdaMat.data(), rows*cols, 1);
 
  // gradient = store each evaluation of the gradient
  // meanGradient/meanGradientSq = store estimates of the gradient and gradient^2 as averages of each evaluation
  // Mt/Vt = used for adam optimiser
  // LB = store evaluations of the ELBO
  Matrix<double, Dynamic, 1> gradient(rows*cols), meanGradient(rows*cols), meanGradientSq(rows*cols), 
                            Mt(rows*cols), Vt(rows*cols), LB(maxIter);
  
  Mt.fill(0);
  Vt.fill(0);
  LB.fill(0);
  // Loop control
  int iter = 0;
  double diff = threshold + 1;
  double meanLB = 0;
  double omeanLB;
  if(false){
  while(diff > threshold | iter < 100) {
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // Calculate Gradient and ELBO
    // Set estimates of the derivative to zero
    meanGradient.fill(0);
    meanGradientSq.fill(0);
    // stores estimates of the ELBO
    double LBrep;
    for(int m = 0; m < M; ++m){
      // Eigen doesn't seem to have random normal sampling, so sample through Rcpp and convert the result to Eigen's format.
      Rcpp::NumericMatrix eps(dim, 1);
      eps(Rcpp::_, 0) = Rcpp::rnorm(dim);
      Map<MatrixXd> epsilon(Rcpp::as<Map<MatrixXd> >(eps));
      // Push epsilon values into the ELBO
      ELBO f(y, epsilon, meanfield, dim);
      // Reset partial derivatives to zero
      stan::math::set_zero_all_adjoints();
      // Store the estimate of f(lambda) in LBrep and df(lambda) / dlambda in gradient
      stan::math::gradient(f, lambda, LBrep, gradient);
      // meanGradient_i^(m) = mean(gradient_i^(m)), similar for meanGradientSq
      for(int i = 0; i < rows*cols; ++i){
        meanGradient(i) += gradient(i) / M;
        meanGradientSq(i) += pow(gradient(i), 2) / M;
      }
      // LB = mean(LB^(m))
      LB(iter-1) += LBrep / M;
    }
    // Update Parameters
    for(int i = 0; i < rows*cols; ++i){
      // Biased estimates of E(gradient) and E(gradient^2)
      Mt(i) = beta1 * Mt(i)  +  (1 - beta1) * meanGradient(i);
      Vt(i) = beta2 * Vt(i)  +  (1 - beta2) * meanGradientSq(i);
      // Correct Bias
      double MtHat = Mt(i) / (1 - pow(beta1, iter));
      double VtHat = Vt(i) / (1 - pow(beta2, iter));
      // Seems to work better if you skip the update on the first iteration
      if(iter > 1){
        // adam updating rule, smaller alpha for variance parameters
        if(i >= dim){
          lambda(i) += alpha / 5 * MtHat / (sqrt(VtHat) + e);
        } else {
          lambda(i) += alpha * MtHat / (sqrt(VtHat) + e);
        }
      }
    }
    // Due to noise in the ELBO evaluations, we look to see if there has been a change in the average value over the last 5 iterations
    if(iter > 5){
      omeanLB = meanLB;
      meanLB = 0.2 * (LB(iter-1) + LB(iter-2) + LB(iter-3) + LB(iter-4) + LB(iter-5));
      diff = std::fabs(meanLB - omeanLB);
    }
    // Print Results
    if(iter % 50 == 0){
      Rcpp::Rcout << "Iteration: " << iter << std::endl << "ELBO: " << LB(iter-1) << std::endl;
    }
  }
  }
  if(iter <= maxIter){
    // iter goes up by one before checking the maxIter condition, so need to use LB(iter-2)
    Rcpp::Rcout << "Converged after " << iter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  } else {
    Rcpp::Rcout << "Warning, failed to converge after " << maxIter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  }
  // Convert to useful format in R
  MatrixXd Mu = lambda.topRows(dim);
  Map<MatrixXd> U(lambda.bottomRows(dim * dim).data(), dim, dim);
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("U") = U,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}