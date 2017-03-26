// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <list>
using namespace Rcpp;
using namespace arma;
using namespace std;
#define Pi 3.14159265358979

class Distribution{
public:
  virtual void set_values(int, int) {};
  virtual double logdens(double x) {return -99.99;};
  virtual vec sample (int n_samples) {return zeros<vec>(n_samples);};
};

class Normal: public Distribution{
  double mean, variance;
  public:
    Normal(double, double);
    void set_values(double m, double v){
      mean = m;
      variance = v;
    }
    double logdens (double x) {
      return -0.5 * log(2*3.141593*variance) - pow(x - mean, 2) / (2 * variance);}
    vec sample (int n_samples){
      return  mean + sqrt(variance) * randn<vec>(n_samples);
    }
    double MeanDeriv(double x){
      return (x - mean) / variance;
    }
    double VarDeriv(double x){
      return (pow(x - mean, 2)/(2*pow(variance, 2)) - 1 / (2*variance));
    }
    double StdDeriv(double x){
      return pow(x - mean, 2)/(pow(variance, 1.5)) - 1 / sqrt(variance); 
    }
    double XDeriv(double x){
      return (mean - x) / variance;
    }
};

Normal::Normal(double m, double v){
  mean = m;
  variance = v;
}

class InverseGamma: public Distribution{
    double shape, scale;
  public:
    InverseGamma(double, double);
    void set_values(double sh, double sc){
      shape = sh;
      scale = sc;
    }
    double logdens (double x) {
      double cons = shape * log(scale) - log(tgamma(shape));
      double kernel = -(shape + 1) * log(x) - scale / x;
      return cons + kernel;
    }
    vec sample (int n_samples){
      return 1.0 / randg<vec>(n_samples, distr_param(shape, 1/scale));
    }
};

InverseGamma::InverseGamma(double sh, double sc){
  shape = sh;
  scale = sc;
}



