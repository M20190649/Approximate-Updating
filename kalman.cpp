// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
//#include <algorithm>
//#include <queue>
//#include <vector>
//#include <functional>   
#include <cmath>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
rowvec KalmanUpdate(double al, double pl, double d, double R, double Q, double c, double Z, double H, double y){
  double a1 = d + R*al;
  double P1 = R*pl*R + Q;
  double v = y - c - Z*a1;
  double F = Z*P1*Z+H;
  double M = P1*Z;
  double a2 = a1 + M*v/F;
  double P2 = P1 - M*M/F;
  rowvec out(2);
  out(0) = a2;
  out(1) = P2;
  return out;
}

// [[Rcpp::export]]
mat KalmanLoop(double a0, double p0, double d, double R, double Q, double c, double Z, double H, vec y) {
  int T = y.size();
  mat out(T, 2);
  out.row(0)= KalmanUpdate(a0, p0, d, R, Q, c, Z, H, y[0]);
  for(int i = 1; i < T; ++i) {
    out.row(i) = KalmanUpdate(out(i-1, 0), out(i-1, 1), d, R, Q, c, Z, H, y[i]);
  }
  return out;
}

// [[Rcpp::export]]
double backwardssampler(double adraw, double at, double pt, double d, double R, double Q) {
  double Vstar = adraw - d - R*at;
  double Fstar = R*pt*R + Q;
  double Mstar = pt*R;
  double Astar = at + Mstar*Vstar/Fstar;
  double Pstar = pt - Mstar*Mstar/Fstar;
  double out = Astar + sqrt(Pstar)*randn<vec>(1)[0];
  return out;
}

// [[Rcpp::export]]
vec backwardsloop(mat kalman, double d, double R, double Q) {
  int T = kalman.n_rows;
  vec out(T);
  out(T-1) = kalman(T-1, 0) + sqrt(kalman(T-1, 1))*randn<vec>(1)[0];
  for(int i = T-2; i > -1; --i) {
    out(i) = backwardssampler(out(i+1), kalman(i, 0), kalman(i, 1), d, R, Q);
  }
  return out;
}