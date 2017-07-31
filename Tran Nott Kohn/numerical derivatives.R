double IGlogdens (double unif, double par1, double par2){
  boost::math::inverse_gamma_distribution<> dist(exp(par1), exp(par2));
  double x = quantile(dist, unif);
  par1 = exp(par1);
  par2 = exp(par2);
  double density = par1 * log(par2) - lgamma(par1) - (par1+1)*log(x) - par2 / x;
  return density;
}

double Blogdens (double unif, double par1, double par2){
  boost::math::beta_distribution<> dist(exp(par1), exp(par2));
  double x = quantile(dist, unif);
  par1 = exp(par1);
  par2 = exp(par2);
  double density = lgamma(par1 + par2) - lgamma(par1) - lgamma(par2) + (par1-1)*log(x) + (par2-1)*log(1-x);
  return density;
}

double Nlogdens(double unif, double par1, double par2){
  boost::math::normal_distribution<> dist(par1, exp(par2));
  double x = quantile(dist, unif);
  par2 = exp(par2);
  double density = - 0.5 * log(2 * 3.14159) -  log(par2)  - pow(x - par1, 2) / (2 * pow(par2, 2));
  return density;
}

Matrix<double, Dynamic, 1> QDeriv(rowvec unif, MatrixXd lambda){
  double h = 0.00001;
  Matrix<double, Dynamic, 1> output(8);
  output(0) = (IGlogdens(unif(0), lambda(0)+h, lambda(1)) - IGlogdens(unif(0), lambda(0), lambda(1)))/h;
  output(1) = (IGlogdens(unif(0), lambda(0), lambda(1)+h) - IGlogdens(unif(0), lambda(0), lambda(1)))/h;
  output(3) = (Blogdens(unif(1), lambda(2)+h, lambda(3)) - Blogdens(unif(1), lambda(2), lambda(3)))/h;
  output(4) = (Blogdens(unif(1), lambda(2), lambda(3)+h) - Blogdens(unif(1), lambda(2), lambda(3)))/h;
  output(5) = (Nlogdens(unif(2), lambda(4)+h, lambda(5)) - Nlogdens(unif(2), lambda(4), lambda(5)))/h;
  output(6) = (Nlogdens(unif(2), lambda(4), lambda(5)+h) - Nlogdens(unif(2), lambda(4), lambda(5)))/h;
  output(7) = (Nlogdens(unif(3), lambda(6)+h, lambda(7)) - Nlogdens(unif(3), lambda(6), lambda(7)))/h;
  output(8) = (Nlogdens(unif(3), lambda(6), lambda(7)+h) - Nlogdens(unif(3), lambda(6), lambda(7)))/h;
  return output;
}

Matrix<double, Dynamic, 1> FDeriv (rowvec unif, MatrixXd theta, MatrixXd lambda){
  double h = 0.00001;
  Matrix<double, Dynamic, 1> output(8);
  
  boost::math::inverse_gamma_distribution<> sigSqDist1(exp(lambda(0+h)), exp(lambda(1)));
  boost::math::beta_distribution<> phiDist1(exp(lambda(2+h)), exp(lambda(3)));
  boost::math::normal_distribution<> muDist1(lambda(4+h), exp(lambda(5)));
  boost::math::normal_distribution<> xtDist1(lambda(6+h), exp(lambda(7)));
  boost::math::inverse_gamma_distribution<> sigSqDist2(exp(lambda(0)), exp(lambda(1+h)));
  boost::math::beta_distribution<> phiDist2(exp(lambda(2)), exp(lambda(3+h)));
  boost::math::normal_distribution<> muDist2(lambda(4), exp(lambda(5+h)));
  boost::math::normal_distribution<> xtDist2(lambda(6), exp(lambda(7+h)));
  
  output(0) = (quantile(sigSqDist1, unif(0)) - theta(0))/h;
  output(1) = (quantile(sigSqDist2, unif(0)) - theta(0))/h;
  output(2) = (quantile(phiDist1, unif(1)) - theta(1))/h;
  output(3) = (quantile(phiDist2, unif(1)) - theta(1))/h;
  output(4) = (quantile(muDist1, unif(2)) - theta(2))/h;
  output(5) = (quantile(muDist2, unif(2)) - theta(2))/h;
  output(6) = (quantile(xtDist1, unif(3)) - theta(3))/h;
  output(7) = (quantile(xtDist2, unif(3)) - theta(3))/h;
  
  return output;
}


if(false){
  // method 1: find two points in xNew closest to xT
  int lowerIndex = 0;
  int upperIndex = 0;
  for(int i = 1; i < N; ++i){
    if(xNew(lowerIndex) < xNew(i) & xNew(i) < xT){
      lowerIndex = i;
    }
    if(xT < xNew(i) & xNew(i) < xNew(upperIndex)){
      upperIndex = i;
    }
    // linear interpolation between the above points to find log(p(x_T | theta, y_{1:T}))
    xTDens = log(pi(lowerIndex))  +  (xT - xNew(lowerIndex)) * (log(pi(upperIndex)) - log(pi(lowerIndex))) / 
      (xNew(upperIndex) - xNew(lowerIndex));
  }
  
  // method 2: find point closest to xT
  int index = 0;
  for(int i = 1; i < N; ++i){
    if(abs(xNew(i) - xT) < abs(xNew(index) - xT)){
      index = i;
    }
  }
  xTDens = log(pi(index));
}