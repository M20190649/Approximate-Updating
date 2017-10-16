// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace std;
using namespace arma;

// [[Rcpp::export]]
mat aug(mat L3, mat full){
  int N = L3.n_rows;
  int M = full.n_rows;
  mat out(N, 3);
  rowvec vals(3);
  for(int i = 0; i < N; ++i){
    int proceed = L3(i, 11);
    if(proceed % 10000 == 0){
      vals = {9999, 9999, 9999};
      out.row(i) = vals;
    } else {
      int range = 0;
      bool found = false;
      bool exact = false;
      while(!found){
        range += 10000;
        if(range > M){
          found = true;
          range = M-1;
        } else {
          found = full(range, 0) >= proceed;
          exact = full(range, 0) == proceed;
        }
      }
      double time = L3(i, 12);
      if(exact){
        for(int j = range - 1000; j < range + 1000; ++j){
          if(full(j, 0) == proceed * 1.0 & full(j, 11) == time){
            double relV = L3(i, 10) - full(j, 9);
            double dist = full(j, 6) - L3(i, 7);
            double ttCol = dist / relV;
            vals = {relV, dist, ttCol};
            out.row(i) = vals;
            break;
          }
          vals = {9999, 9999,9999};
          out.row(i) = vals;
        }
      } else {
        for(int j = range - 10000; j < range; ++j){
          if(full(j, 0) == proceed * 1.0 & full(j, 11) == time){
            double relV = L3(i, 10) - full(j, 9);
            double dist = full(j, 6) - L3(i, 7);
            double ttCol = dist / relV;
            vals = {relV, dist, ttCol};
            out.row(i) = vals;
            break;
          }
          vals = {9999, 9999,9999};
          out.row(i) = vals;
        }
      }
    } 
    if(i % 25000 == 0){
      Rcpp::Rcout << i << std::endl;
    }
  }
  return out;
}

// [[Rcpp::export]]
vec timeToChange(mat L3){
  int N = L3.n_rows;
  vec out(N);
  for(int i = 0; i < N - 1; ++i){
    if(L3(i, 2) == 0){
      out(i) = 9999;
    } else {
      for(int j = i + 1; j < N; ++j){
        if(L3(i, 0) != L3(j, 0)){
          out(i) = 9999;
          break;
        }
        if(L3(i, 9) != L3(j, 9)){
          out(i) = j - i;
          break;
        }
      }
    }
  }
  out(N-1) = 9999;
  return out;
}

// [[Rcpp::export]]
vec timeSinceChange(mat L3){
  int N = L3.n_rows;
  vec out(N);
  out(0) = 9999;
  for(int i = 1; i < N - 1; ++i){
    if(L3(i, 2) == 0){
      out(i) = 9999;
    } else {
      for(int j = i - 1; j >=  0; --j){
        if(L3(i, 0) != L3(j, 0)){
          out(i) = 9999;
          break;
        }
        if(L3(i, 9) != L3(j, 9)){
          out(i) = i - j;
          break;
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
vec countLaneChange(mat cars, vec uniqueID){
  int N = uniqueID.n_elem;
  int M = cars.n_rows;
  vec output (N, fill::zeros);
  int IDcounter = 0;
  int numChange = 0;
  for(int i = 1; i < M; ++i){
    if(cars(i, 0) != cars(i-1, 0)){
      output(IDcounter) = numChange;
      IDcounter += 1;
      numChange = 0;
      continue;
    }
    if(cars(i, 16) != cars(i-1, 16)){
      numChange += 1;
    }
  }
  return output;
}





