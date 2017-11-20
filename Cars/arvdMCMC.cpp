// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat ARVD_MCMC(mat data, vec hyper, int reps, int lags){
  mat theta (reps, 2 + 2 * lags);
  int T = data.n_rows;
  theta(0, 0) = theta(0, 1) = 0.00001;
  for(int i = 0; i < lags; ++i){
    theta(0, 2 + 2 * i) = theta(0, 3 + 2 * i) = 0.2;
  }
  double vShape = hyper(0) + (T - lags) / 2;
  double dShape = hyper(2) + (T - lags) / 2;
  for(int iter = 1; iter < reps; ++iter){
    // sigSqV
    double vScale = hyper(1);
    for(int t = lags; t < T; ++t){
      double temp = data(t, 0);
      for(int i = 1; i <= lags; ++i){
        temp -= theta(iter-1, 1 + i) * data(t-i, 1); 
      }
      vScale += pow(temp, 2) / 2;
    }
    theta(iter, 0) = 1.0 / randg<vec>(1, distr_param(vShape, 1.0 / vScale))[0];
    
    // sigSqD
    double dScale = hyper(3);
    for(int t = lags; t < T; ++t){
      double temp = data(t, 1);
      for(int i = 1; i <= lags; ++i){
        temp -= theta(iter-1, 1 + lags + i) * data(t-i, 1); 
      }
      dScale += pow(temp, 2) / 2;
    }
    theta(iter, 1) = 1.0 / randg<vec>(1, distr_param(dShape, 1.0 / dScale))[0];
    
    // phi terms
    for(int i = 1; i <= lags; ++i){
      // V
      double meanNumer = theta(iter, 0) * hyper(2 + 2*i);
      double meanDenom = theta(iter, 0);
      for(int t = lags; t < T; ++t){
        meanNumer += hyper(3 + 2*i) * data(t, 0) * data(t-i, 0);
        meanDenom += hyper(3 + 2*i) * pow(data(t-i, 0), 2);
        
        for(int k = 1; k <= lags; ++k){
          if(k != i){
            if(k < i){
              meanNumer -= hyper(3 + 2*i) * theta(iter, 1 + k) * data(t-k, 0) * data(t-i, 0);
            } else {
              meanNumer -= hyper(3 + 2*i) * theta(iter-1, 1 + k) * data(t-k, 0) * data(t-i, 0);
            }
            
          }
        }
      }
      double var = theta(iter, 0) * hyper(3 + 2*i) / meanDenom;
      // Stationary conditions for AR1, AR2
      double draw;
      bool nonStat;
      if(lags == 1){
        nonStat = true;
        while(nonStat){
          draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
          if(draw < 1 & draw > -1){
            nonStat = false;
          }
        }
      } else if(lags == 2){
        if(i == 1){
          int counter = 0;
          nonStat = true;
          while(nonStat){
            draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
            counter += 1;
            if(counter > 20){
              if(draw + theta(iter-1, 3) > 1){
                draw = 0.99 - theta(iter-1, 3);
              } else if(theta(iter-1, 3) - draw < 1){
                draw = theta(iter-1, 3) - 0.99;
              }
            }
            if(draw + theta(iter-1, 3) < 1 & theta(iter-1, 3) - draw < 1){
              nonStat = false;
            }
          }
        } else {
          int counter = 0;
          nonStat = true;
          while(nonStat){
            draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
            counter += 1;
            if(counter > 20){
              if(draw + theta(iter, 2) > 1){
                draw = 0.99 - theta(iter, 2);
              } else if(theta(iter, 2) - draw < 1){
                draw = theta(iter, 2) - 0.99;
              }
              if(draw > 1){
                draw = 0.99;
              } else if(draw < -1){
                draw = -0.99;
              }
            }
            if(draw + theta(iter, 2) < 1 & draw - theta(iter, 2) < 1 & draw < 1 & draw > -1){
              nonStat = false;
            }
          }
        }
      } else {
        draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
      }
      theta(iter, 1 + i) = draw;
      // D
      meanNumer = theta(iter, 1) * hyper(2 + 2*lags + 2*i);
      meanDenom = theta(iter, 1);
      for(int t = lags; t < T; ++t){
        meanNumer += hyper(3 + 2*lags + 2*i) * data(t, 1) * data(t-i, 1);
        meanDenom += hyper(3 + 2*lags + 2*i) * pow(data(t-i, 1), 2);
        for(int k = 1; k <= lags; ++k){
          if(k != i){
            if(k < i){
              meanNumer -= hyper(3 + 2*lags + 2*i) * theta(iter, lags + 1 + k) * data(t-k, 1) * data(t-i, 1);
            } else {
              meanNumer -= hyper(3 + 2*lags + 2*i) * theta(iter-1, lags + 1 + k) * data(t-k, 1) * data(t-i, 1);
            }
          }
        }
      }
      var = theta(iter, 1) * hyper(3 + 2*lags + 2*i) / meanDenom;
      // Stationary conditions for AR1, AR2
      if(lags == 1){
        nonStat = true;
        while(nonStat){
          draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
          if(draw < 1 & draw > -1){
            nonStat = false;
          }
        }
      } else if(lags == 2){
        if(i == 1){
          int counter = 0;
          nonStat = true;
          while(nonStat){
            draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
            counter += 1;
            if(counter > 20){
              if(draw + theta(iter-1, lags + 3) > 1){
                draw = 0.99 - theta(iter-1, lags + 3);
              } else if(theta(iter-1, lags + 3) - draw < 1){
                draw = theta(iter-1, lags + 3) - 0.99;
              }
            }
            if(draw + theta(iter-1, lags+3) < 1 & theta(iter-1, lags+3) - draw < 1){
              nonStat = false;
            }
          }
        } else {
          nonStat = true;
          int counter = 0;
          while(nonStat){
            draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
            counter += 1;
            if(counter > 20){
              if(draw + theta(iter, lags + 2) > 1){
                draw = 0.99 - theta(iter, lags + 2);
              } else if(theta(iter, lags + 2) - draw < 1){
                draw = theta(iter, lags + 2) - 0.99;
              }
              if(draw > 1){
                draw = 0.99;
              } else if(draw < -1){
                draw = -0.99;
              }
            }
            if(draw + theta(iter, lags+2) < 1 & draw - theta(iter, lags+2) < 1 & draw < 1 & draw > -1){
              nonStat = false;
            }
          }
        }
      } else {
        draw = meanNumer / meanDenom  +  pow(var, 0.5) * randn<vec>(1)[0];
      }
      theta(iter, lags + 1 + i)  = draw;
    }
  }
  return theta;
}