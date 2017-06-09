// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp> // Making use of many of the distributions in Boost
#include "likelihood.h"
using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost::math;



// Boost does not provide a location-scale t distribution, so it is implemented manually
struct locScaleT {
  double loc, scale, df;
  locScaleT (double l, double s, double d){
    loc = l;
    scale = s;
    df = d;
  }
};

// Truncnorm is also implemented manually
struct truncNorm {
  double mean, var;
  truncNorm (double m, double v){
    mean = m;
    var = v;
  }
};

// As is truncated loc/scale T
struct truncT {
  double loc, scale, df;
  truncT (double l, double s, double d){
    loc = l;
    scale = s;
    df = d;
  }
};

// Overloading boost's pdf and quantile function for the new distributions
double pdf(locScaleT dist, double x){
  boost::math::students_t_distribution<> aux(dist.df);
  return 1.0/dist.scale * pdf(aux, (x-dist.loc)/dist.scale);
}

double quantile(locScaleT dist, double p){
  boost::math::students_t_distribution<> auxDist(dist.df);
  return quantile(auxDist, p)*dist.scale + dist.loc;
}

double pdf(truncNorm dist, double x){
  boost::math::normal_distribution<> aux(dist.mean, sqrt(dist.var));
  double area = cdf(aux, 1); 
  // While we should measure the area from -1 to 1, the probability mass below -1 is extremely close to zero
  return pdf(aux, x) / area;
}

double quantile(truncNorm dist, double p){
  boost::math::normal_distribution<> aux(dist.mean, sqrt(dist.var));
  double area = cdf(aux, 1);
  return quantile(dist, p*area);
}

double pdf(truncT dist, double x){
  boost::math::students_t_distribution<> aux(dist.df);
  double area = cdf(aux, (1.0-dist.loc)/dist.scale);
  return pdf(aux, (x-dist.loc)/dist.scale) / (dist.scale*area);
}

double quantile(truncT dist, double p){
  boost::math::students_t_distribution<> aux(dist.df);
  double area = cdf(aux, (1.0-dist.loc)/dist.scale);
  return quantile(aux, area*p)*dist.scale + dist.loc;
}

double pLogDens(vec &y, vec x, double sigmaSqY, double sigmaSqX, double phi, double gamma){
  int TS = y.n_elem;
  double prior = - pow(gamma, 2) / 200  -  2 * log(sigmaSqY)  -  2 * log(sigmaSqX)  -  1.0/sigmaSqY  -  1.0/sigmaSqX  +
    19.0 * log(1 + phi)  +  0.5 * log(1 - phi);
  double states = 0.5*log(1-pow(phi, 2))  -  0.5*log(sigmaSqX)  -  (1-pow(phi, 2)) * pow(x[0]-gamma, 2) / (2*sigmaSqX);
  double data = 0;
  for(int t = 1; t < TS+1; ++t){
    states -= pow(x[t]-gamma-phi*(x[t-1]-gamma), 2) / (2*sigmaSqX);
    data -= pow(y[t-1]-x[t], 2) / (2*sigmaSqY);
  }
  return data  +  states  +  prior;
}

// Transform a vector of uniforms to a vector of realisations from a single variable
vec SingleMarginalTransform (vec & unifs, vec & thetaDist, double & xDist, mat & thetaParams, mat & xParams, int T, int j){
  vec output;
  int n = unifs.n_elem;
  if(j <= 1){
    if(thetaDist[j]==1){
      boost::math::weibull_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    } else if(thetaDist[j]==2){
     boost::math::gamma_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    } else if(thetaDist[j]==3){
     boost::math::inverse_gamma_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    } else {
     boost::math::lognormal_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    }
  } else if (j == 2){
    if(thetaDist[j]==1){
      truncNorm dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    } else {
      truncT dist(thetaParams(0, j), thetaParams(1, j), thetaParams(2, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    }
  } else if(j == 3) {
      if(thetaDist[j]==1){
        boost::math::normal_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
        for(int i = 0; i < n; ++i){
          output[i] = quantile(dist, unifs[i]);
        }
    } else {
      locScaleT dist(thetaParams(0, j), thetaParams(1, j), thetaParams(2, j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    }
  } else {
    if(xDist == 2){
     boost::math::normal_distribution<> dist(xParams(0, T+5-j), sqrt(xParams(1, T+5-j)));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    } else {
      locScaleT dist(xParams(0, T+5-j), xParams(1, T+5-j), xParams(2, T+5-j));
      for(int i = 0; i < n; ++i){
        output[i] = quantile(dist, unifs[i]);
      }
    }
  }

  return output;
}

// Transform a matrix of uniforms to a matrix of realisations of theta & xt
mat JointMarginalTransform (mat & unifs, vec & thetaDist, double & xDist, mat & thetaParams, mat & xParams, int T){
  int n = unifs.n_rows;
  int p = unifs.n_cols;
  mat output(n, p);
  
  for(int j = 0; j < 2; ++j){
    if(thetaDist[j]==1){
      boost::math::weibull_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output(i, j) = quantile(dist, unifs(i, j));
      }
    } else if(thetaDist[j]==2){
      boost::math::gamma_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output(i, j) = quantile(dist, unifs(i, j));
      }
    } else if(thetaDist[j]==3){
      boost::math::inverse_gamma_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output(i, j) = quantile(dist, unifs(i, j));
      }
    } else {
      boost::math::lognormal_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      for(int i = 0; i < n; ++i){
        output(i, j) = quantile(dist, unifs(i, j));
      }
    }
  }
  int j = 2;
  if(thetaDist[j]==1){
    truncNorm dist(thetaParams(0, j), thetaParams(1, j));
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  } else {
    truncT dist(thetaParams(0, j), thetaParams(1, j), thetaParams(2, j));
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  }
  j = 3;
  if(thetaDist[j]==1){
    boost::math::normal_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  } else {
    locScaleT dist(thetaParams(0, j), thetaParams(1, j), thetaParams(2, j));
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  }
  for(int j = 4; j < p; ++j){
    if(xDist == 2){
      boost::math::normal_distribution<> dist(xParams(0, T+5-j), sqrt(xParams(1, T+5-j)));
      for(int i = 0; i < n; ++i){
        output(i, j) = quantile(dist, unifs(i, j));
      }
    } else {
      locScaleT dist(xParams(0, T+5-j), xParams(1, T+5-j), xParams(2, T+5-j));
      for(int i = 0; i < n; ++i){
        output(i, j) = quantile(dist, unifs(i, j));
      }
    }
  }
  return output;
}

// Evaluate q(x, theta) - Requires BiCopPdf to be implemented in c++
double QLogDens (vec & unifsdep, vec & sims, vec & thetaDist, double & xDist, mat & thetaParams,
                 mat & xParams, mat & VineMatrix, mat & VineFamily, mat & VinePar1, mat & VinePar2, int T){
  double margins = 0;
  int p = unifsdep.n_elem;
  
  for(int j = 0; j < 2; ++j){
    if(thetaDist[j]==1){
      boost::math::weibull_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      margins += log(pdf(dist, sims[j]));
    } else if(thetaDist[j]==2){
      boost::math::gamma_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      margins += log(pdf(dist, sims[j]));
    } else if(thetaDist[j]==3){
      boost::math::inverse_gamma_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      margins += log(pdf(dist, sims[j]));
    } else {
      boost::math::lognormal_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
      margins += log(pdf(dist, sims[j]));
    }
  }
  int j = 2;
  if(thetaDist[j]==1){
    truncNorm dist(thetaParams(0, j), thetaParams(1, j));
    margins += log(pdf(dist, sims[j]));
  } else {
    truncT dist(thetaParams(0, j), thetaParams(1, j), thetaParams(2, j));
    margins += log(pdf(dist, sims[j]));
  }
  j = 3;
  if(thetaDist[j]==1){
    boost::math::normal_distribution<> dist(thetaParams(0, j), thetaParams(1, j));
    margins += log(pdf(dist, sims[j]));
  } else {
    locScaleT dist(thetaParams(0, j), thetaParams(1, j), thetaParams(2, j));
    margins += log(pdf(dist, sims[j]));
  }
  
  for(int j = 4; j < p; ++j){
    if(xDist == 2){
      boost::math::normal_distribution<> dist(xParams(0, T+5-j), sqrt(xParams(1, T+5-j)));
      margins += log(pdf(dist, sims[j]));
    } else {
      locScaleT dist(xParams(0, T+5-j), xParams(1, T+5-j), xParams(2, T+5-j));
      margins += log(pdf(dist, sims[j]));
    }
  }
  double copulas = 0;
  for(int i = T; i < T+5; ++i){
    for(int j = 0; j < T+5; ++j){
      if(VineFamily(i, j) != 0){
        int u1 = VineMatrix(j, j) -1;
        int u2 = VineMatrix(i, j) -1;
        double loglik = 0;
        int n = 1;
        LL(&VineFamily(i, j), &n, &unifsdep[u1], &unifsdep[u2], &VinePar1(i, j), &VinePar2(i, j), &loglik);
        copulas += loglik;
      }
    }
  }
  return margins + copulas;
}

// Evaluate ELBO - Requires QLogDens
double ELBO (mat & unifsdep, mat & sims, vec & y, vec & thetaDist, double & xDist, mat & thetaParams, mat & xParams,
                     mat & VineMatrix, mat & VineFamily, mat & VinePar, mat & VinePar2, int & N){
  double elbo = 0;
  int T = y.n_elem;
  for(int i = 0; i < N; ++i){
    vec simscol = sims.row(i).t();
    vec unifsdepcol = unifsdep.row(i).t();
    vec x = simscol.tail(T+1);
    double sigmaSqX = exp(simscol[0]);
    double sigmaSqY = exp(simscol[1]);
    elbo += pLogDens(y, x, sigmaSqX, sigmaSqY, simscol[2], simscol[3]) -
      QLogDens(unifsdepcol, simscol, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, T);
  }
  return elbo/N;
}

// Simulate from the Vine - Requires H & HInv Functions
// [[Rcpp::export]]
mat VineSim (mat & unifs, mat & VineMatrix, mat & VineFamily, mat & VinePar, mat & VinePar2, int T){
  int n = T + 5;
  int nSamples = unifs.n_rows;
  mat m = VineMatrix;
  // Re-arrange Vine Matrix so diagonal is from n to 1.
  vec mTheta = {m(T+1, T+1), m(T+2, T+2), m(T+3, T+3), m(T+4, T+4)};
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; ++j){
      if(m(T+i, T+j) == mTheta(0)){
        m(T+i, T+j) = 4;
      } else if(m(T+i, T+j) == mTheta(1)){
        m(T+i, T+j) = 3;
      } else if(m(T+i, T+j) == mTheta(2)){
        m(T+i, T+j) = 2;
      } else if(m(T+i, T+j) == mTheta(3)){
        m(T+i, T+j) = 1;
      }
    }
  }
 
  mat mmax(n, n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      mmax(i, j) = m.row(j).tail(n-i).max();
    }
  }
  
  cube vd(nSamples, n, n, fill::zeros), vind(nSamples, n, n, fill::zeros);
  vd.slice(n-1) = unifs;
  cube z1(nSamples, n, n, fill::zeros), z2(nSamples, n, n, fill::zeros);
  mat x(nSamples, n, fill::ones);
  x.col(0) = unifs.col(0);
  int nVine = 1;
  
  // Vine Algorithm looping through number of samples to take
  for(int j = 0; j < nSamples; ++j){
    for(int k = n-2; k >= 0; --k){
      for(int i = max(k+1, n-5); i < n; ++i){
        if(mmax(i, k) == m(i, k)){
          z2(j, i, k) = vd(j, i, n-mmax(i, k)-1);
        } else {
          z2(j, i, k) = vind(j, i, n-mmax(i, k)-1);
        }
        Hinv(&VineFamily(i, k), &nVine, &vd(j, n-1, k), &z2(j, i, k), &VinePar(i, k), &VinePar2(i, k), &vd(j, n-1, k));
      }
      
      x(j, n-k-1) = vd(j, n-1, k);
      for(int i = n-1; i >= max(k+1, n-6); --i){
        z1(j, i, k) = vd(j, i, k);
        Hfunc(&VineFamily(i, k), &nVine, &z1(j, i, k), &z2(j, i, k), &VinePar(i, k), &VinePar2(i, k), &vd(j, i-1, k));
        vind(j, i-1, k) = vd(j, i-1, k);
      }
    }
  }
  
  // Re-arrange rows to match original Vine Matrix input.
  mat output(nSamples, n);
  output.submat(0, 4, nSamples-1, n-1) = x.submat(0, 4, nSamples-1, n-1);
  output.col(mTheta[3]) = x.col(0);
  output.col(mTheta[2]) = x.col(1);
  output.col(mTheta[1]) = x.col(2);
  output.col(mTheta[0]) = x.col(3);
  return output;
}

// dlog(p(y, x, theta)) / dx and dtheta
vec pDeriv (vec &y, vec &x, double &sigmaSqY, double &sigmaSqX, double &phi, double &gamma, int &T, int &S, bool &xderiv){
  vec derivs(T+S+5, fill::zeros);
  derivs[0] = -0.5 * (T+S+4) / sigmaSqY  +  pow(sigmaSqY, -2);
  derivs[1] = -0.5 * (T+S+5) / sigmaSqX  +  pow(sigmaSqX, -2)  + 
    (1-pow(phi, 2)) * pow(x[0]-gamma, 2) / (2*pow(sigmaSqX, 2));
  derivs[2] = -phi / (1-pow(phi, 2))  +   phi / sigmaSqX * pow(x[0]-gamma, 2)  +
    19.0 * phi / (1 + phi)  -  0.5 * phi / (1 - phi);
  derivs[3] = -gamma/100  +  (1-pow(phi, 2)) * (x[0]-gamma) / sigmaSqX;
  for(int t = 1; t < T+S+1; ++t){
    derivs[0] += pow(y[t-1]-x[t], 2) / (2*pow(sigmaSqY, 2));
    derivs[1] += pow(x[t]-gamma-phi*(x[t-1]-gamma), 2) / (2*pow(sigmaSqX, 2));
    derivs[2] += (x[t-1]-gamma) * (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX;
    derivs[3] += (1-phi) * (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX;
  }
  if(!xderiv){
    return derivs;
  }
  if(S == 0){
    derivs[4] = (phi*x[1]-x[0]+gamma*(1-phi)) / sigmaSqX;
    for(int i = 5; i < T+4; ++i){
      int t = i - 4;
      derivs[i] = (y[t-1]-x[t]) / sigmaSqY  -  (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX  +
        phi * (x[t+1]-gamma-phi*(x[t]-gamma)) / sigmaSqX;
    }
  } else {
    for(int i = T+5; i < T+S+4; ++i){
      int t = i - 4;
      derivs[i] = (y[t-1]-x[t]) / sigmaSqY  -  (x[t]-gamma-phi*(x[t-1]-gamma)) / sigmaSqX +
        phi * (x[t+1]-gamma-phi*(x[t]-gamma)) / sigmaSqX;
    }
  }
  derivs[T+S+4] = (y[T+S-1]-x[T+S]) / sigmaSqY  -  (x[T+S]-gamma-phi*(x[T+S-1]-gamma)) / sigmaSqX;
  return derivs;
}

// Take Partial Derivatives - Requires all of the above. 
// Need to come up with non-numeric approach because eta derivatives all involve re-simulation from the vine.
// Unless the fully c++ vine sim is super fast.
double Partial(mat & unifs, mat & unifsdep, mat & sims, vec & y, int T, int S, vec & thetaDist, double & xDist, mat & thetaParams,
               mat &xParams, mat & VineMatrix, mat & VineFamily, mat & VinePar, mat & VinePar2, int M, int what, int i, int j, bool par1 = true){
  double h = 0.0000001;
  double deriv = 0;
  if(what == 1){
    mat thetaParams2 = thetaParams;
    thetaParams2(i, j) += h;
    mat newsim = sims;
    vec unidep = unifsdep.row(j).t();
    newsim.col(j) = SingleMarginalTransform(unidep, thetaDist, xDist, thetaParams2, xParams, T+S, j);
    for(int i = 0; i < M; ++i){
      vec newsimcol = newsim.row(i).t();
      vec simcol = sims.row(i).t();
      vec unidepcol = unifsdep.row(i).t();
      deriv += (pLogDens(y, newsimcol.tail(T+1), exp(newsimcol[0]), exp(newsimcol[1]), newsimcol[2], newsimcol[3]) -
        pLogDens(y, simcol.tail(T+1), exp(simcol[0]), exp(simcol[1]), simcol[2], simcol[3]) -
        QLogDens(unidepcol, newsimcol, thetaDist, xDist, thetaParams2, xParams, VineMatrix, VineFamily, VinePar, VinePar2, T+S) +
        QLogDens(unidepcol, simcol, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, T+S))/h;
    }
  } else if (what == 2){
    mat xParams2 = xParams;
    xParams(i, j) += h;
    mat newsim = sims;
    vec unifsdepcol = unifsdep.row(j+5+T).t();
    newsim.col(j+5+T) = SingleMarginalTransform(unifsdepcol, thetaDist, xDist, thetaParams, xParams2, T+S, j+5+T);
    for(int i = 0; i < M; ++i){
      vec newsimcol = newsim.row(i).t();
      vec simcol = sims.row(i).t();
      vec unidepcol = unifsdep.row(i).t();
      deriv += (pLogDens(y, newsimcol.tail(T+1), exp(newsimcol[0]), exp(newsimcol[1]), newsimcol[2], newsimcol[3]) -
        pLogDens(y, simcol.tail(T+1), exp(simcol[0]), exp(simcol[1]), simcol[2], simcol[3]) -
        QLogDens(unidepcol, newsimcol, thetaDist, xDist, thetaParams, xParams2, VineMatrix, VineFamily, VinePar, VinePar2, T+S) +
        QLogDens(unidepcol, simcol, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, T+S))/h;
    }
  } else {
    mat Vine2Par = VinePar;
    mat Vine2Par2 = VinePar2;
    if (what == 3) {
      if (par1){
        Vine2Par(T+S+2+i, T+S+1+j) += h;
      } else {
        Vine2Par2(T+S+2+i, T+S+1+j) += h;
      }
    } else if(what == 4) {
      if (par1){
        Vine2Par(T+S+1+i, j) += h;
      } else {
        Vine2Par2(T+S+1+i, j) += h;
      }
    } else {
      if (par1){
        Vine2Par(T+S+1, j) += h;
      } else {
        Vine2Par2(T+S+1, j) += h;
      }
    }
    mat unifsSub = unifs.submat(0, 0, M-1, T+S+5);
    mat newdep = VineSim(unifsSub, VineMatrix, VineFamily, VinePar, VinePar2, T+S);
    mat newsim = JointMarginalTransform(newdep, thetaDist, xDist, thetaParams, xParams, T+S);
    for(int i = 0; i < M; ++i){
      vec newsimcol = newsim.row(i).t();
      vec simcol = sims.row(i).t();
      vec unidepcol = unifsdep.row(i).t();
      deriv += (pLogDens(y, newsimcol.tail(T+1), exp(newsimcol[0]), exp(newsimcol[1]), newsimcol[2], newsimcol[3]) -
        pLogDens(y, simcol.tail(T+1), exp(simcol[0]), exp(simcol[1]), simcol[2], simcol[3]) -
        QLogDens(unidepcol, newsimcol, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, Vine2Par, Vine2Par2, T+S) +
        QLogDens(unidepcol, simcol, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, T+S))/h;
    }
  }
  return deriv / M;
}

//[[Rcpp::export]]
Rcpp::List CopulaSGA (vec y, int S, vec thetaDist, double xDist, mat thetaParams, mat xParams1, mat VineMatrix, mat VineFamily, mat VinePar, 
      mat VinePar2, int M, int maxIter, double threshold=0.01, double alpha=0.01, double beta1=0.9, double beta2=0.999, double e = 0.0000001){
  // initialisation {
  int T = y.n_elem - S;
  int dimx = 2;
  if(xDist == 2) {
    dimx = 3;
  }
  mat xParams(dimx, T+S+1, fill::randu);
  xParams.submat(0, S, dimx-1, T+S) = xParams1;

  mat MtLamT(3, 4, fill::zeros), VtLamT(3, 4, fill::zeros),
      MtLamX(xDist+1, S, fill::zeros), VtLamX(xDist+1, S, fill::zeros),
      MtEtaT1(4, 4, fill::zeros), VtEtaT1(4, 4, fill::zeros), MtEtaT2(4, 4, fill::zeros), VtEtaT2(4, 4, fill::zeros),
      MtEtaX1(4, S, fill::zeros), VtEtaX1(4, S, fill::zeros),  MtEtaX2(4, S, fill::zeros), VtEtaX2(4, S, fill::zeros);
  vec MtEtaXX1(S, fill::zeros), VtEtaXX1(S, fill::zeros), MtEtaXX2(S, fill::zeros), VtEtaXX2(S, fill::zeros);
  
  mat FullVineMatrix(T+S+5, T+S+5, fill::zeros), FullVineFamily(T+S+5, T+S+5, fill::zeros),
      FullVinePar(T+S+5, T+S+5, fill::zeros), FullVinePar2(T+S+5, T+S+5, fill::zeros);
  for(int i = 0; i < T+S+1; ++i){
    FullVineMatrix(i, i) = 5 + i;
  }
  for(int i = 1; i < T+S+1; ++i){
    for(int j = 0; j < i; ++j){
      FullVineMatrix(i, j) = i - j + 4;
    }
  }
  
  FullVineMatrix.submat(S, S, T+S+4, T+S+4) = VineMatrix;
  FullVineFamily.submat(S, S, T+S+4, T+S+4) = VineFamily;
  FullVinePar.submat(S, S, T+S+4, T+S+4) = VinePar;
  FullVinePar2.submat(S, S, T+S+4, T+S+4) = VinePar2;
  for(int i = 0; i < S; ++i){
    for(int t = T+S+1; t < T+S+5; ++t){
      FullVineMatrix(t, i) = FullVineMatrix(t, S);
    }
    for(int t = T+S; t < T+S+5; ++t){
      FullVineFamily(t, i) = FullVineFamily(t, S);
      FullVinePar(t, i) = FullVinePar(t, S);
      FullVinePar2(t, i) = FullVinePar2(t, S);
    }
  }
  VineMatrix = FullVineMatrix;
  VineFamily = FullVineFamily;
  VinePar = FullVinePar;
  VinePar2 = FullVinePar2;
  //}
  
  // loop control {
  int iter = 0;
  mat unifs = randu<mat>(max(25, M), T+S+5);
  mat unifsdep = VineSim(unifs, VineMatrix, VineFamily, VinePar, VinePar2, T+S);
  mat sims = JointMarginalTransform(unifsdep, thetaDist, xDist, thetaParams, xParams, T+S);
  vec LB(maxIter + 1);
  int LBrep = 25;
  LB(iter) = ELBO(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, LBrep);
  double meanLB = 1;
  double meanLBold;
  double diff = threshold + 1;
  
  while(diff > threshold){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // }
    
    // derivatives {
    mat PLamT(3, 4, fill::zeros);
    mat PLamX(xDist+1, S, fill::zeros);
    mat PEtaT1, PEtaT2(4, 4, fill::zeros);
    mat PEtaX1, PEtaX2(4, S, fill::zeros);
    vec PEtaXX1, PEtaXX2(S, fill::zeros);
    
    for(int j = 0; j < 4; ++j){
      int u = 2;
      if(j>1){
        u = thetaDist[j] + 1;
      }
      for(int i = 0; i < u; ++i){
        PLamT(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
              xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 1, i, j);
      }
    }
    for(int j = 0; j < S; ++j){
      for(int i = 0; i <= xDist; ++i){
        PLamX(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
              xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 2, i, j);
      }
    }
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j <= i; ++j){
        if(VinePar(T+S+2+i, T+S+2+j) != 0){
          PEtaT1(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 3, i, j);
        }
        if(VinePar2(T+S+2+i, T+S+2+j) != 0){
          PEtaT2(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 3, i, j, false);
        }
      }
    }
    for(int j = 0; j < S; ++j){
      for(int i = 0; i < 4; ++i){
        if(VinePar(T+S+1+i, j) != 0){
          PEtaX1(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 4, i, j);
        }
        if(VinePar2(T+S+1+i, j) != 0){
          PEtaX2(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 4, i, j, false);
        }
      }
    }
    for(int j = 0; j < S; ++j){
      if(VinePar(T+S, j) != 0){
        PEtaXX1(j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
               xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 5, 0, j);
      }
      if(VinePar2(T+S, j) != 0){
        PEtaXX2(j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
               xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 5, 0, j, false);
      }
    }
    // }
    
    // updating parameters {
    MtLamT = beta1*MtLamT + (1-beta1)*PLamT;
    MtLamX = beta1*MtLamX + (1-beta1)*PLamX;
    MtEtaT1 = beta1*MtEtaT1 + (1-beta1)*PEtaT1;
    MtEtaT1 = beta1*MtEtaT2 + (1-beta1)*PEtaT2;
    MtEtaX1 = beta1*MtEtaX1 + (1-beta1)*PEtaX1;
    MtEtaX2 = beta1*MtEtaX2 + (1-beta1)*PEtaX2;
    MtEtaXX1 = beta1*MtEtaXX1 + (1-beta1)*PEtaXX1;
    MtEtaXX2 = beta1*MtEtaXX2 + (1-beta1)*PEtaXX2;
      
    VtLamT = beta2*VtLamT + (1-beta2)*pow(PLamT, 2);
    VtLamX = beta2*VtLamX + (1-beta2)*pow(PLamX, 2);
    VtEtaT1 = beta2*VtEtaT1 + (1-beta2)*pow(PEtaT1, 2);
    VtEtaT1 = beta2*VtEtaT2 + (1-beta2)*pow(PEtaT2, 2);
    VtEtaX1 = beta2*VtEtaX1 + (1-beta2)*pow(PEtaX1, 2);
    VtEtaX2 = beta2*VtEtaX2 + (1-beta2)*pow(PEtaX2, 2);
    VtEtaXX1 = beta2*VtEtaXX1 + (1-beta2)*pow(PEtaXX1, 2);
    VtEtaXX2 = beta2*VtEtaXX2 + (1-beta2)*pow(PEtaXX2, 2);
    
    thetaParams += alpha * (MtLamT/(1-pow(beta1, iter))) / (sqrt(VtLamT/(1-pow(beta2, iter)))+e);
    xParams.submat(0, 0, xDist+1, S-1) += alpha * (MtLamX/(1-pow(beta1, iter))) / (sqrt(VtLamX/(1-pow(beta2, iter)))+e);
    VinePar.submat(T+1, T+1, T+4, T+4) += alpha * (MtEtaT1/(1-pow(beta1, iter))) / (sqrt(VtEtaT1/(1-pow(beta2, iter)))+e);
    VinePar2.submat(T+1, T+1, T+4, T+4) += alpha * (MtEtaT2/(1-pow(beta1, iter))) / (sqrt(VtEtaT2/(1-pow(beta2, iter)))+e);
    VinePar.submat(T+1, 1, T+4, S-1) += alpha * (MtEtaX1/(1-pow(beta1, iter))) / (sqrt(VtEtaX1/(1-pow(beta2, iter)))+e);
    VinePar2.submat(T+1, 1, T+4, S-1) += alpha * (MtEtaX2/(1-pow(beta1, iter))) / (sqrt(VtEtaX2/(1-pow(beta2, iter)))+e);
    VinePar.submat(T, 0, T, S-1) += alpha * (MtEtaXX1/(1-pow(beta1, iter))) / (sqrt(VtEtaXX1/(1-pow(beta2, iter)))+e);
    VinePar2.submat(T, 0, T, S-1) += alpha * (MtEtaXX2/(1-pow(beta1, iter))) / (sqrt(VtEtaXX2/(1-pow(beta2, iter)))+e);
    
    // }
    
    // loop control {
    unifs = randu<mat>(max(25, M), T+5);
    unifsdep = VineSim(unifs, VineMatrix, VineFamily, VinePar, VinePar2, T);
    sims = JointMarginalTransform(unifsdep, thetaDist, xDist, thetaParams, xParams, T);
    LB(iter) = ELBO(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, LBrep);
    if(iter % 5 == 0){
      meanLBold = meanLB;
      meanLB = 0.2 * (LB[iter]+LB[iter-1]+LB[iter-2]+LB[iter-3]+LB[iter-4]);
      diff = abs(meanLB - meanLBold);
    }
    if(iter % 100 == 0 & diff > threshold){
      Rcpp::Rcout << iter << std::endl;
    }
  } // close while loop
  // }
  
  
  // returning output {
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  } else {
    LB = LB.head(iter);
  }
  Rcpp::Rcout << iter << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("Theta") = thetaParams,
                            Rcpp::Named("X") = xParams,
                            Rcpp::Named("Vine") = Rcpp::List::create(
                              Rcpp::Named("Matrix") = VineMatrix,
                              Rcpp::Named("Family") = VineFamily,
                              Rcpp::Named("Par") = VinePar,
                              Rcpp::Named("Par2") = VinePar2));
  // }
}
      
