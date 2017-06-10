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

struct ELBO {
  const Matrix<double, Dynamic, 1> y;
  const Matrix<double, Dynamic, 1> epsilon;
  const bool meanfield;
  ELBO(const MatrixXd& yin, const MatrixXd& epsin, const bool mf) : y(yin), epsilon(epsin), meanfield(mf) {}
  template <typename T>
  T operator()(const Matrix<T, Dynamic, 1>& lambda) 
    const {
    using std::pow; using std::log; using std::exp; using std::sqrt;
    int Obs = y.size() - 1;
    int dim = epsilon.size();
    T sigmaSq;
    T phi;
    T gamma;
    T plogDens;
    T qlogDens;
    // AR 1 model without a mean
    if(dim == 2){
      sigmaSq = exp(lambda(0) + lambda(2)*epsilon(0));
      phi = lambda(1) + lambda(4)*epsilon(0) + lambda(5)*epsilon(1);
      // p log density
      plogDens = -0.5 * (Obs+5) * log(sigmaSq)  -  1.0 / sigmaSq;
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(y(t) - phi*y(t-1), 2) / (2*sigmaSq);
      }
      // q log density
      qlogDens = lambda(0)  +  lambda(2) * epsilon(0)  +  log(sqrt(pow(lambda(2), 2)))  +  log(sqrt(pow(lambda(5), 2)));
      qlogDens -= 0.5 * (pow(epsilon(0), 2) + pow(epsilon(1), 2));
    // AR 1 model with a mean
    } else if(dim==3){
      sigmaSq = exp(lambda(0) + lambda(3)*epsilon(0));
      phi = lambda(1) + lambda(6)*epsilon(0) + lambda(7)*epsilon(1);
      gamma = lambda(2) + lambda(9)*epsilon(0) + lambda(10)*epsilon(1) + lambda(11)*epsilon(2);
      // p log density
      plogDens = -0.5 * (Obs+5) * log(sigmaSq)  -  1.0 / sigmaSq;  -   pow(gamma, 2) / 200;
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(y(t) - gamma - phi*(y(t-1)-gamma), 2) / (2*sigmaSq);
      }
      // q log density
      qlogDens = lambda(0)  +  lambda(3) * epsilon(0);
      for(int i = 0; i < 3; ++i){
        qlogDens += log(sqrt(pow(lambda((i+1)*dim+i), 2)))  -  0.5 * pow(epsilon(i), 2);
      }
    // Basic State-Space model
    } else {
      Obs += 1;
      Matrix<T, Dynamic, 1> transformed(dim);
      for(int i = 0; i < dim; ++i){
        transformed(i) = lambda(i);
        if(meanfield){
          transformed(i) += lambda((i+1)*dim+i) * epsilon(i);
        } else {
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
        1.0 / sigmaSqX;  -   pow(gamma, 2) / 200;
      for(int t = 1; t < Obs+1; ++t){
        plogDens -= pow(x(t) - gamma - phi*(x(t-1)-gamma), 2) / (2*sigmaSqX);
        plogDens -= pow(y(t-1) - x(t), 2) / (2*sigmaSqY);
      }
      // q log density
      qlogDens = lambda(0)  +  lambda(dim) * epsilon(0)  +  lambda(1)  +  lambda(2*dim+1) * epsilon(1);
      if(!meanfield){
        qlogDens += lambda(2*dim) * epsilon(0); 
      }
      for(int i = 0; i < 3; ++i){
        qlogDens += log(sqrt(pow(lambda((i+1)*dim+i), 2)))  -  0.5 * pow(epsilon(i), 2);
      }
    }
    return plogDens - qlogDens;
  }
};

// [[Rcpp::export]]
Rcpp::List AD_SGA(const Rcpp::NumericMatrix yIn, const Rcpp::NumericMatrix lambdaIn, int lrows,
                   int M, int maxIter, bool meanfield = false, 
                   double threshold=0.001, double alpha=0.1, double beta1=0.9, double beta2=0.999){
  // Convert to Eigen Format
  Map<MatrixXd> y(Rcpp::as<Map<MatrixXd> >(yIn));
  Map<MatrixXd> lambdaMat(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  
  double e = pow(10, -8);
  // Initialise new matrices
  int lcols = lrows + 1;
  Map<MatrixXd> lambda(lambdaMat.data(), lcols*lrows, 1);
  Matrix<double, Dynamic, 1> gradient(lrows*lcols), meanGradient(lrows*lcols), meanGradientSq(lrows*lcols), 
                            Mt(lrows*lcols), Vt(lrows*lcols), LB(maxIter);
  Mt.fill(0);
  Vt.fill(0);
  LB.fill(0);
  // Loop control
  int iter = 0;
  double diff = threshold + 1;
  while(diff > threshold | iter < 50) {
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // Calculate Gradient and LB
    meanGradient.fill(0);
    meanGradientSq.fill(0);
    double LBrep;
    for(int m = 0; m < M; ++m){
      Rcpp::NumericMatrix eps(lrows, 1);
      eps(Rcpp::_, 0) = Rcpp::rnorm(lrows);
      Map<MatrixXd> epsilon(Rcpp::as<Map<MatrixXd> >(eps));
      ELBO f(y, epsilon, meanfield);
      stan::math::set_zero_all_adjoints();
      stan::math::gradient(f, lambda, LBrep, gradient);
      for(int i = 0; i < lrows*lcols; ++i){
        meanGradient(i) += gradient(i) / M;
        meanGradientSq(i) += pow(gradient(i), 2) / M;
      }
      LB(iter-1) += LBrep / M;
    }
    // Update Parameters
    for(int i = 0; i < lrows*lcols; ++i){
      Mt(i) = beta1 * Mt(i)  +  (1 - beta1) * meanGradient(i);
      Vt(i) = beta2 * Vt(i)  +  (1 - beta2) * meanGradientSq(i);
      double MtHat = Mt(i) / (1 - pow(beta1, iter));
      double VtHat = Vt(i) / (1 - pow(beta2, iter));
      if(iter > 1){
        lambda(i) += alpha * MtHat / (sqrt(VtHat) + e);
      }
    }
    if(iter > 1){
      diff = sqrt(pow(LB(iter-1) - LB(iter-2), 2));
    }
    // Print Results
    if(iter % 50 == 0){
      Rcpp::Rcout << "Iteration: " << iter << std::endl << "ELBO: " << LB(iter-1) << std::endl;
    }
  }
  // Convert to useful format
  MatrixXd Mu = lambda.topRows(lrows);
  Map<MatrixXd> U(lambda.bottomRows(lrows * lrows).data(), lrows, lrows);
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("U") = U,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}