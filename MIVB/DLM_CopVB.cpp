// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp> // Making use of many of the distributions in Boost
using namespace Rcpp;
using namespace arma;
using namespace std;

// Boost does not provide a location-scale t distribution, so it is implemented manually
struct locScaleT {
  double loc, scale, df;
  locScaleT locScaleT (double, double, double);
};

locScaleT locScaleT (double l, double s, double d){
  loc = l;
  scale = s;
  df = d;
}

// Overloading boost's pdf and quantile function for the new distribution
double pdf(locScaleT dist, double x){
  student_t_distribution auxDist(dist.df);
  return 1.0/dist.scale * pdf(auxDist, (x-dist.loc)/dist.scale);
}

double quantile(locScaleT dist, double p){
  student_t_distribution auxDist(dist.df);
  return quantile(auxDist, p)*dist.scale + dist.loc;
}

// Evaluate p(y, x, theta)
double PLogDens (vec & y, vec & sims){
  inverse_gamma_distribution sigSqY(1, 1);
  inverse_gamma_distribution sigSqX(1, 1);
  normal_distribution gamma(0, 10);
  double dataDens = 0;
  double priorDens = log(pdf(sigSqY, sims[0])) + log(pdf(sigSqX, sims[1])) + log(1/2) + log(pdf(gamma(sims[3])));
  normal_distribution x0(0, sqrt(sims[1] / (1 - pow(sims[2], 2))));
  priorDens += log(pdf(x0, sims[4]));
  for(int t = 5; t < T+5; ++5){
    normal_distribution xt(sims[2]*sims[t-1], sqrt(sims[1]));
    normal_distribution yt(sims[3]+sims[t], sqrt(sims[0]));
    priorDens += log(pdf(xt, sims[t]));
    dataDens += log(pdf(yt, y[t-5]));
  }
  return dataDens + priorDens;
}

// Transform a vector of uniforms to a vector of realisations from a single variable
vec SingleMarginalTransform (vec & unifs, vec & thetaDist, double & xDist, mat & thetaParmas, mat & xParams, int T, int j){
  vec output;
  int n = unifs.n_elem;
  if(j <= 1){
    if(thetaDist[j]==1){
      weibull_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else if(thetaDist[j]==2){
      gamma_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else if(thataDist[j]==3){
      inverse_gamma_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else {
      lognormal_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    }
  } else if (j <= 3){
    if(thetaDist[j]==1){
      normal_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else {
      locScaleT dist(thetaParams[0, j], thetaParams[1, j], thetaParams[2, j]);
    }
  } else {
    if(xDist == 2){
      normal_distribution dist(xParams[0, T+5-j], sqrt(xParams[1, T+5-j]));
    } else {
      locScaleT dist(xParams[0, T+5-j], xParams[1, T+5-j], xParams[2, T+5-j]);
    }
  }
  for(int i = 0; i < n; ++i){
    output[i] = quantile(dist, unifs[i]);
  }
  return output;
}

// Transform a matrix of uniforms to a matrix of realisations of theta & xt
mat JointMarginalTransform (mat & unifs, vec & thetaDist, double & xDist, mat & thetaParmas, mat & xParams, int T){
  int n = unifs.n_rows;
  int p = unifs.n_cols;
  mat output(n, p);
  
  for(int j = 0; j < 2; ++j){
    if(thetaDist[j]==1){
      weibull_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else if(thetaDist[j]==2){
      gamma_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else if(thataDist[j]==3){
      inverse_gamma_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else {
      lognormal_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    }
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  }
  for(int j = 2; j < 4; ++j){
    if(thetaDist[j]==1){
      normal_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else {
      locScaleT dist(thetaParams[0, j], thetaParams[1, j], thetaParams[2, j]);
    }
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  }
  for(int j = 4; j < p; ++j){
    if(xDist == 2){
      normal_distribution dist(xParams[0, T+5-j], sqrt(xParams[1, T+5-j]));
    } else {
      locScaleT dist(xParams[0, T+5-j], xParams[1, T+5-j], xParams[2, T+5-j]);
    }
    for(int i = 0; i < n; ++i){
      output(i, j) = quantile(dist, unifs(i, j));
    }
  }
  return output;
}

// Evaluate q(x, theta) - Requires BiCopPdf to be implemented in c++
double QLogDens (vec & unifsdep, vec & sims, vec & thetaDist, double & xDist, mat & thetaParams,
                 mat & xParams, mat & VineMatrix, mat & VineFamily, mat & VinePar1, mat & VinePar2){
  double margins = 0;
  int p = unifsdep.n_elem;
  
  for(int j = 0; j < 2; ++j){
    if(thetaDist[j]==1){
      weibull_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else if(thetaDist[j]==2){
      gamma_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else if(thataDist[j]==3){
      inverse_gamma_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else {
      lognormal_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    }
    margins += log(pdf(dist, sims[j]));
  }
  for(int j = 2; j < 4; ++j){
    if(thetaDist[j]==1){
      normal_distribution dist(thetaParams[0, j], thetaParams[1, j]);
    } else {
      locScaleT dist(thetaParams[0, j], thetaParams[1, j], thetaParams[2, j]);
    }
    margins += log(pdf(dist, sims[j]));
  }
  for(int j = 4; j < p; ++j){
    if(xDist == 2){
      normal_distribution dist(xParams[0, T+5-j], sqrt(xParams[1, T+5-j]));
    } else {
      locScaleT dist(xParams[0, T+5-j], xParams[1, T+5-j], xParams[2, T+5-j]);
    }
    margins += log(pdf(dist, sims[j]));
  }
  double copulas = 0;
  for(int i = T; i < T+5; ++i){
    for(int j = 0; j < T+5; ++j){
      if(VineFamily(i, j) != 0){
        int u1 = VineMatrix(j, j) -1;
        int u2 = VineMatrix(i, j) -1;
        copulas += log(BiCopPDF(unifsdep[u1], unifsdep[u2], VineFamily(i, j), 
                                VinePar(i, j), VinePar2(i, j)));
      }
    }
  }
  return margins + copulas;
}

// Evaluate ELBO - Requires QLogDens
double ELBO (mat & unifsdep, mat & sims, vec & y, vec & thetaDist, double & xDist, mat & thetaParams, mat & xParams,
                     mat & VineMatrix, mat & VineFamily, mat & VinePar, mat & VinePar2, int N){
  double elbo = 0;
  for(int i = 0; i < N; ++i){
    elbo += PLogDens(y, sims.row(i)) - QLogDens(unifsdep.row(i), sims.row(i), thetaDist, xDist, thetaParams, xParams, 
                     VineMatrix, VineFamily, VinePar, VinePar2);
  }
  return elbo/N;
}

// Simulate from the Vine - Requires H & HInv Functions
mat VineSim function(mat & unifs, mat & VineMatrix, mat & VineFamily, mat & VinePar, mat & VinePar2, int T){
  int n = unifs.n_col;
  int nr = unifs.n_row;
  mat m = VineMatrix;
  vec mTheta = (m.diag).tail(4);
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; ++j){
      switch(m[T+i, T+j]){
        case 0 :
          break;
        case mTheta[0]:
          m[T+i, T+j] = 4;
          break;
        case mTheta[1];
          m[T+i, T+j] = 3;
          break;
        case mTheta[2];
          m[T+i, T+j] = 2;
          break;
        case mTheta[3];
          m[T+i, T+j] = 1;
          break;
      }
    }
  }
 
  mat mmax(n, n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      mmax(i, j) = m.row(j).tail(n-i).max();
    }
  }
  cube vd(nr, n, n, fill::zeros);
  cube vind(nr, n, n, fill::zeros);
  vd.slice(n-1) = unifs;
  cube z1(nr, n, n, fill::zeros);
  cube z2(nr, n, n, fill::zeros);
  mat x(nr, n);
  x.col(0) = unifs.col(0);
  
  for(int k = n-2; k >= 0; --k){
    for(int i = k+1; i < n; ++i){
      if(mmax(i, k) == m(i, k)){
        z2.subcube(0, i, k, nr-1, i, k) = vd.subcube(0, i, n-mmax(i, k), nr-1, i, n-mmax(i, k));
      } else {
        z2.subcube(0, i, k, nr-1, i, k) = vind.subcube(0, i, n-mmax(i, k), nr-1, i, n-mmax(i, k));
      }
      for(int j = 0; j < nr; ++j){
        vd(j, n-1, k) = BiCopHInv(vd(j, n-1, k), z2(j, i, k), VineFamily(i, k), VinePar(i, k), VinePar2(i, k));
      }
    }
    for(int j = 0; j < nr; ++j){
      x(j, n-k-1) = vd(j, n-1, k);
    }
    for(int i = n-1; i >= k+1; --i){
      z1.subcube(0, i, k, nr-1, i, k) = vd.subcube(0, i, k, nr-1, i, k);
      for(int j = 0; j < nr; ++j){
        vd(j, i-1, k) = vind(j, i-1, k) = BiCopH(z1(j, i, k), z2(j, i, k), VineFamily(i, k), VinePar(i, k), VinePar2(i, k));
      }
    }
  }
  mat output(nr, n);
  output.submat(0, 4, nr, n) = x.submat(0, 4, nr, n);
  output.col(mTheta[3]) = x.col(0);
  output.col(mTheta[2]) = x.col(1);
  output.col(mTheta[1]) = x.col(2);
  output.col(mTheta[0]) = x.col(3);
  return output;
}


// Take Partial Derivatives - Requires all of the above. 
// Need to come up with non-numeric approach because eta derivatives all involve re-simulation from the vine.
// Unless the fully c++ vine sim is super fast.
double Partial(mat & unifs, mat & unifsdep, mat & sims, vec & y, int T, int S, vec & thetaDist, double & xDist, mat & thetaParams,
               mat &xParams, mat & VineMatrix, mat & VineFamily, mat & VinePar, mat & VinePar2, int M, char what, int i, int j bool par1 = true){
  double h = 0.0000001;
  double deriv = 0;
  if(what == 'LamT'){
    mat thetaParams2 = thetaParams;
    thetaParams2(i, j) += h;
    vec newsim = sims;
    newsim.col(j) = SingleMarginalTransform(unifsdep.col(j), thetaDist, xDist, thetaParams2, xParams, T, j);
    for(int i = 0; i < M; ++i){
      deriv += (PLogDens(y, newsim.row(i)) - PLogDens(y, sims.row(i)) -
        QLogDens(unifsdep.row(i), newsim.row(i), thetaDist, xDist, thetaParams2, xParams, VineMatrix, VineFamily, VinePar, VinePar2) +
        QLogDens(unifsdep.row(i), sims.row(i), thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2))/h;
    }
  } else if (what == 'LamX'){
    mat xParams2 = xParams;
    xParams(i, j) += h;
    newsim = sims;
    newsim.col(j+5+T-S) = SingleMarginalTransform(unifsdep.col(j+5+T-S), thetaDist, xDist, thetaParams, xParams2, T, j+5+T-S);
    for(int i = 0; i < M; ++i){
      deriv += (PLogDens(y, newsim.row(i)) - PLogDens(y, sims.row(i)) -
        QLogDens(unifsdep.row(i), newsim.row(i), thetaDist, xDist, thetaParams, xParams2, VineMatrix, VineFamily, VinePar, VinePar2) +
        QLogDens(unifsdep.row(i), sims.row(i), thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2))/h;
    }
  } else {
    mat Vine2Par = VinePar;
    mat Vine2Par2 = VinePar2;
    if (what == 'EtaT') {
      if (par1){
        Vine2Par(T+2+i, T+1+j) += h;
      } else {
        Vine2Par2(T+2+i, T+1+j) += h;
      }
    } else if(what == 'EtaX') {
      if (par1){
        Vine2Par(T+1+i, j) += h;
      } else {
        Vine2Par2(T+1+i, j) += h;
      }
    } else {
      if (par1){
        Vine2Par(T+1, j) += h;
      } else {
        Vine2Par2(T+1, j) += h;
      }
    }
    mat newdep = VineSim(unifs.submat(0, 0, M-1, T+5), VineMatrix, VineFamily, VinePar, VinePar2, T);
    mat newsim = JointMarginalTransform(newdep, thetaDist, xDist, thetaParams, xParams);
    for(int i = 0; i < M; ++i){
      deriv += (PLogDens(y, newsim.row(i)) - PLogDens(y, sims.row(i)) -
        QLogDens(newdep.row(i), newsim.row(i), thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, Vine2Par, Vine2Par2) +
        QLogDens(unifsdep.row(i), sims.row(i), thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2))/h;
    }
  }
  return deriv / M;
}

//[[Rcpp::export]]
Rcpp::List CopulaSGA (vec y, int S, vec thetaDist, double xDist, mat thetaParams, mat xParams, mat VineMatrix, mat VineFamily, mat VinePar, 
      mat VinePar2, int M, int maxIter, double threshold=0.01, double alpha=0.01, double beta1=0.9, double beta2=0.999, double e = 0.0000001){
  int T = y.n_elem;
  mat temp = xParams;
  mat xParams(xDist+1, T, fill::randu);
  xParams.submat(0, S, xDist+1, T) = temp;
  delete temp;
  
  mat MtLamT, VtLamT(3, 4, fill::zeros);
  mat MtLamX, VtLamX(xDist+1, S, fill::zeros);
  mat MtEtaT1, VtEtaT1, MtEtaT2, VtEtaT2(4, 4, fill::zeros);
  mat MtEtaX1, VtEtaX1,  MtEtaX2, VtEtaX2(4, S, fill::zeros);
  vec MtEtaXX1, VtEtaXX1, MtEtaXX2, VtEtaXX2(S, fill::zeros
  
  FullVineMatrix(T+5, T+5, fill::zeros);
  FullVineFamily(T+5, T+5, fill::zeros);
  FullVinePar(T+5, T+5, fill::zeros);
  FullVinePar2(T+5, T+5, fill::zeros);
  for(int i = 0; i < T+1; ++i){
    FullVineMatrix(i, i) = 5 + i;
  }
  for(int i = 1; i < T+1){
    for(int j = 0; j < i){
      FullVineMatrix(i, j) = i - j + 4;
    }
  }
  FullVineMatrix.submat(S, S, T+4, T+4) = VineMatrix;
  FullVineFamily.submat(S, S, T+4, T+4) = VineFamily;
  FullVinePar.submat(S, S, T+4, T+4) = VinePar;
  FullVinePar2.submat(S, S, T+4, T+4) = VinePar2;
  for(int i = 0; i < S; ++i){
    for(int t = T+1; t < T+5; ++t){
      FullVineMatrix(t, i) = FullVineMatrix(t, S);
    }
    for(int t = T; t < T+5; ++t){
      FullVineFamily(t, i) = FullVineFamily(t, S);
      FullVinePar(t, i) = FullVinePar(t, S);
      FullVinePar2(t, i) = FullVinePar2(t, S);
    }
  }
  VineMatrix = FullVineMatrix;
  VineFamily = FullVineFamily;
  VinePar = FullVinePar;
  VinePar2 = FullVinePar2;
  delete FullVineMatrix;
  delete FullVineFamily;
  delete FullVinePar;
  delete FullVinePar2;
  
  int iter = 0;
  mat unifs = randu<mat>(max(25, M), T+5);
  mat unifsdep = VineSim(unifs, VineMatrix, VineFamily, VinePar, VinePar2, T);
  mat sims = JointMarginalTransform(unifsdep, thetaDist, xDist, thetaParams, xParams, T);
  vec LB(maxIter + 1);
  LB(iter) = ELBO(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, 25);
  double lastDiff = threshold + 1;
  double meandiff = 0;
  while(lastDiff > 5*threshold & meanDiff > threshold){
    iter += 1;
    if(iter > maxIter){
      break;
    }
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
              xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'LamT', i, j);
      }
    }
    for(int j = 0; j < S; ++j){
      for(int i = 0; i <= xDist; ++i){
        PLamX(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
              xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'LamX', i, j);
      }
    }
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j <= i; ++j){
        if(VinePar(T+2+i, T+2+j) != 0){
          PEtaT1(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'EtaT', i, j);
        }
        if(VinePar2(T+2+i, T+2+j) != 0){
          PEtaT2(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'EtaT', i, j, false);
        }
      }
    }
    for(int j = 0; j < S; ++j){
      for(int i = 0; i < 4; ++i){
        if(VinePar(T+1+i, j) != 0){
          PEtaX1(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'EtaX', i, j);
        }
        if(VinePar2(T+1+i, j) != 0){
          PEtaX2(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
                 xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'EtaX', i, j, false);
        }
      }
    }
    for(int j = 0; j < S; ++j){
      if(VinePar(T, j) != 0){
        PEtaXX1(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
               xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'EtaXX', 0, j);
      }
      if(VinePar2(T, j) != 0){
        PEtaXX2(i, j) = Partial(unifs, unifsdep, sims, y, T, S, thetaDist, xDist, thetaParams,
               xParams, VineMatrix, VineFamily, VinePar, VinePar2, M, 'EtaXX', 0, j, false);
      }
    }
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
    
    unifs = randu<mat>(max(25, M), T+5);
    unifsdep = VineSim(unifs, VineMatrix, VineFamily, VinePar, VinePar2, T);
    sims = JointMarginalTransform(unifsdep, thetaDist, xDist, thetaParams, xParams, T);
    LB(iter) = ELBO(unifsdep, sims, y, thetaDist, xDist, thetaParams, xParams, VineMatrix, VineFamily, VinePar, VinePar2, 25);
    if(iter > 10){
      lastDiff = abs(LB(iter) - LB(iter-1));
      meanDiff = mean(LB.subvec(iter-4, iter) - LB.subvec(iter-5, iter-1));
    }
  }
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  } else {
    LB = LB.head(iter);
  }
  return Rcpp::List::create(Rcpp::Named("Theta") = thetaParams,
                            Rcpp::Named("X") = xParams,
                            Rcpp::Named("Vine") = Rcpp::List::create(
                              Rcpp::Named("Matrix") = VineMatrix,
                              Rcpp::Named("Family") = VineFamily,
                              Rcpp::Named('Par') = VinePar,
                              Rcpp::Named('Par2') = VinePar2));
}
      
