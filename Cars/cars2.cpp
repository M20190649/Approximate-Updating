// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <boost/math/distributions.hpp> 

using namespace std;
using namespace arma;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// evaluate and automatially differentiate log(p(x_1:T, z_T+1 theta)) wrt theta and z_T+1
struct logP {
  const mat x;
  const int obs;
  const int N;
  const vec epsilon;
  const vec hyperParams;
  logP(const mat& xIn, const int& obsIn, const int& nIn, const vec& epsIn, const vec& hypIn) :
    x(xIn), obs(obsIn), N(nIn), epsilon(epsIn), hyperParams(hypIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambdaV)
    const{
    using std::log; using std::exp; using std::pow; using std::cos; using std::sin;
    
    // Create theta
    Matrix<T, Dynamic, 1> theta(12);
    for(int i = 0; i < 12; ++i){
      theta(i) = lambdaV(i) + lambdaV(12+i) * epsilon(i);
    }
    // Constrained Positive
    T sigSqA = exp(theta(0)), sigSqD = exp(theta(1)), alpha = exp(theta(4)), beta = exp(theta(5)), sigSqX = exp(theta(7)), sigSqY = exp(theta(8));
    // Constrained to (0, 1)
    T lambda = 1.0 / (1 + exp(-theta(2))), phi = 1.0 / (1 + exp(-theta(3))), sT1 = 1.0 / (1 + exp(-theta(11)));
    // Unconstrained
    T gamma = theta(6), aT1 = theta(9), dT1 = theta(10), pi = 3.14159;
    
    // Evaluate log(p(theta))
    T prior = -(hyperParams(0) + 1) * log(sigSqA)  -  hyperParams(1) / sigSqA  - 
      (hyperParams(2) + 1) * log(sigSqD)  -  hyperParams(3) / sigSqD  +  
      (hyperParams(4) - 1) * log(lambda)  +  (hyperParams(5)-1) * log(1-lambda)  +
      (hyperParams(6) - 1) * log(phi)  +  (hyperParams(7)-1) * log(1-phi)  +  
      (hyperParams(8) - 1) * log(alpha)  -  hyperParams(9) * alpha  +  
      (hyperParams(10) - 1) * log(beta)  -  hyperParams(11) * beta  -
      pow(gamma - hyperParams(12), 2) / (2*hyperParams(13))  -
      (hyperParams(14) + 1) * log(sigSqX)  -  hyperParams(15) / sigSqX -
      (hyperParams(16) + 1) * log(sigSqY)  -  hyperParams(17) / sigSqY; 
    // Set up for particle filtering
    
    T xAllDens, xtDens, omegaSum, ZT1Dens;
    Matrix<T, Dynamic, Dynamic> Old(N, 6), New(N, 6), Resample(N, 6), Ahead(N, 6);
    Matrix<T, Dynamic, 1> pweights(N), piSum(N), omega(N), lweights(N), lambdaSum(N);
    Matrix<int, Dynamic, 1> kDraws(N);
    pweights.fill(1.0/N); // initial weights
    
    // Sample z0 from stationary distribution
    for(int k = 0; k < N; ++k){
      Old(k, 0) = sqrt(sigSqA / (1 - pow(phi, 2))) * randn<vec>(1)[0];
      Old(k, 1) = 0;
      Old(k, 2) = pi / 2  +  sqrt(sigSqD) * randn<vec>(1)[0];
      Old(k, 3) = 5.0 + randn<vec>(1)[0];
      Old(k, 4) = 5.0 + randn<vec>(1)[0];
      Old(k, 5) = randn<vec>(1)[0];
    }
    xAllDens = 0;
    
    // Particle Filter loop
    for(int t = 0; t < obs; ++t){
      // Create CDF to use for resampling
      piSum(0) = pweights(0);
      for(int k = 1; k < N; ++k){
        piSum(k) = pweights(k) + piSum(k-1);
      }
      // Resample xt via inverse CDF
      for(int k = 0; k < N; ++k){
        double u = randu<vec>(1)[0];
        for(int i = 0; i < N; ++i){
          if(u < piSum(i)){
            for(int j = 0; j < 6; ++j){
              Resample(k, j) = Old(i, j);
            }
            //Resample.block(1, 6, k, 6*t) = Old.block(1, 6, i, 6*t);
            break;
          }
        }
      }
      
      // Calculate Step Ahead Means and k density
      for(int k = 0; k < N; ++k){
        Ahead(k, 0) = phi * Resample(k, 0);
        Ahead(k, 1) = Resample(k, 1);
        Ahead(k, 2) = pi/2  +  (1 - Ahead(k, 1)) * lambda * (Resample(k, 2) - pi/2)  +
          Ahead(k, 1) * gamma * (Resample(k, 4) - 5);
        Ahead(k, 3) = Resample(k, 3)  +  Ahead(k, 0);
        Ahead(k, 4) = Resample(k, 4)  +   Ahead(k, 3) * cos(Ahead(k, 2));
        Ahead(k, 5) = Resample(k, 5)  +   Ahead(k, 3) * sin(Ahead(k, 2));
        lweights(k) = 1.0 / sqrt(2 * pi * sigSqX) * exp(-pow(x(t, 0) - Ahead(k, 4), 2) / (2*sigSqX)) *
          1.0 / sqrt(2 * pi * sigSqY) * exp(-pow(x(t, 1) - Ahead(k, 5), 2) / (2*sigSqY));
        if(k == 0){
          lambdaSum(k) = lweights(k);
        } else {
          lambdaSum(k) = lweights(k) + lambdaSum(k-1);
        }
      }
      // Normalise
      for(int k = 0; k < N; ++k){
        lambdaSum(k) = lambdaSum(k) / lambdaSum(N-1);
      }
      
      
      // Draw index k using lambda as weights
      for(int k = 0; k < N; ++k){
        double u = randu<vec>(1)[0];
        for(int i = 0; i < N; ++i){
          if(u < lambdaSum(i)){
            kDraws(k) = i;
          }
        }
      }
      
      // Calculate weights
      omegaSum = 0;
      for(int k = 0; k < N; ++k){
        // Step Ahead transitions
        New(k, 0) = phi * Resample(kDraws(k), 0)  +  sqrt(sigSqA) * randn<vec>(1)[0];
        double u = randu<vec>(1)[0];
        if((Resample(kDraws(k), 1) == 0 & u < 1.0 / (1 + 1/(alpha * abs(Resample(kDraws(k), 4) - 5)))) |
           (Resample(kDraws(k), 1) == 1 & u < 1.0 / (1 + 1/(beta * abs(Resample(kDraws(k), 4) - 5))))){
          New(k, 1) = 1;
        } else {
          New(k, 1) = 0;
        }
        New(k, 2) =  pi/2  +  (1 - New(k, 1)) * lambda * (Resample(kDraws(k), 2) - pi/2)  +
          New(k, 1) * gamma * (Resample(kDraws(k), 4) - 5)  +  sqrt(sigSqD) * randn<vec>(1)[0];
        New(k, 3) = Resample(kDraws(k), 3)  +  New(k, 0);
        New(k, 4) = Resample(kDraws(k), 4)  +  New(k, 3) * cos(New(k, 2));
        New(k, 5) = Resample(kDraws(k), 5)  +  New(k, 3) * sin(New(k, 2));
        // Measurement density
        omega(k) =  1.0 / sqrt(2 * pi * sigSqX) * exp(-pow(x(t, 0) - New(k, 4), 2) / (2*sigSqX)) *
          1.0 / sqrt(2 * pi * sigSqY) * exp(-pow(x(t, 1) - New(k, 5), 2) / (2*sigSqY));
        // sum of weights for normalisation
        omegaSum += omega(k);
      }
      
      // Normalise weights
      for(int k = 0; k < N; ++k){
        pweights(k) = omega(k) / omegaSum;
      }
      // log(p(x_1:T | theta)) = sum_t log(p(x_t | theta))
      xtDens = log(omegaSum / N);
      xAllDens += xtDens;
      // reset vOld as we step from t to t+1 (except last iteration)
      if(t < obs - 1){
        Old.block(N, 6, 0, 6*(t+1)) = New.block(N, 6, 0, 6*t);
      }
    } // end of particle filter loop
    // step ahead latent states density
    ZT1Dens = 0; 
    Matrix<T, Dynamic, 1> atDens(N), dtDens(N), prob1(N), stDens(N);
    // p(V_T+1 | theta, x_1:T) approx sum_k pi_k p(V_T+1 | V_T,k, theta)
    for(int k = 0; k < N; ++k){
      atDens(k) = 1.0 / sqrt(2*pi*sigSqA) * exp(-pow(aT1 - phi * New(k, 0), 2) / (2*sigSqA));
      dtDens(k) = 1.0 / sqrt(2*pi*sigSqD) * exp(-pow(dT1- pi/2 - (1-sT1)*lambda*(New(k, 2) -pi/2) + sT1*gamma*(New(k, 4) - 5), 2) / (2*sigSqD));
      prob1(k) = 1.0 / (1 + 1/(pow(alpha, 1-New(k, 1)) * pow(beta, New(k, 1)) * abs(New(k, 4)-5)));
      stDens(k) = pow(prob1(k), sT1) * pow(1-prob1(k), 1-sT1);
      ZT1Dens += pweights(k) * atDens(k) * dtDens(k) * stDens(k);
    }
    // interested in log density
    ZT1Dens = log(ZT1Dens);
    return xAllDens + ZT1Dens + prior;
  }
};


// evaluate and automatially differentiate log det J
struct logDetJ {
  const vec epsilon;
  logDetJ(const vec epsIn) : // constructor
    epsilon(epsIn){}
  template <typename T> // 
  T operator ()(const Matrix<T, Dynamic, 1>& lambda) // derivative is with respect to lambda
    const{
    using std::log; using std::exp; 
    T logdet = 0;
    for(int i = 0; i < 12; ++i){
      if(i == 6 | i == 9 | i == 10){
        logdet += log(lambda(12+i));
      } else if(i == 0 | i == 1 | i == 4 | i == 5 | i == 7 | i == 8){
        logdet += log(lambda(12+i)) + lambda(i) + lambda(12+i)*epsilon(i);
      } else if(i == 2 | i == 3 | i == 11){
        logdet += -log(lambda(12+i)) + lambda(i) + lambda(12+i)*epsilon(i) - 2*log(1 + exp(lambda(i) + lambda(12+i)+epsilon(i)));
      }
    }
    return logdet;
  }
};

// generate sobol points
mat sobol_points(int N, int D) {
  ifstream infile("new-joe-kuo-6.21201" ,ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 
  
  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [N];
  C[0] = 1;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  mat POINTS(N, D, fill::zeros);
  
  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1
  
  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [N];
  X[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS(i, 0) = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  
  // Clean up
  delete [] V;
  delete [] X;
  
  
  // ----- Compute the remaining dimensions -----
  for (unsigned j=1;j<=D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
        V[i] = V[i-s] ^ (V[i-s] >> s); 
        for (unsigned k=1;k<=s-1;k++) 
          V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
      POINTS(i, j) = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
    }
    
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
  
  return POINTS;
}

// randomly shuffle sobol points
mat shuffle(mat sobol){
  int N = sobol.n_rows;
  int D = sobol.n_cols;
  mat output(N, D, fill::zeros);
  // draw a random rule of: switch 1 and 0  /  do not switch for each binary digit.
  vec rule = randu<vec>(16);
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < D; ++j){
      // grab element of the sobol sequence
      double x = sobol(i, j);
      // convert to a binary representation
      uvec binary(16, fill::zeros);
      for(int k = 1; k < 17; ++k){
        if(x > pow(2, -k)){
          binary(k-1) = 1;
          x -= pow(2, -k);
        }
      }
      // apply the transform of tilde(x_k) = x_k + a_k mod 2, where a_k = 1 if rule_k > 0.5, 0 otherwise
      for(int k = 0; k < 16; ++k){
        if(rule(k) > 0.5){
          binary(k) = (binary(k) + 1) % 2;
        }
      }
      // reconstruct base 10 number from binary representation
      for(int k = 0; k < 16; ++k){
        if(binary(k) == 1){
          output(i, j) += pow(2, -(k+1));
        }
      }
    }
  }
  return output;
}

// generate theta for importance sampling.
vec constructTheta(vec epsilon, MatrixXd lambda){
  // Create theta
  vec theta(12);
  for(int i = 0; i < 12; ++i){
    theta(i) = lambda(i) + lambda(12+i) * epsilon(i);
  }
  // Constrained Positive
  theta(0) = exp(theta(0));
  theta(1) = exp(theta(1));
  theta(4) = exp(theta(4));
  theta(5) = exp(theta(5));
  theta(7) = exp(theta(7));
  theta(8) = exp(theta(8));
  // Constrained to (0, 1)
  theta(2) = 1.0 / (1 + exp(-theta(2)));
  theta(3) = 1.0 / (1 + exp(-theta(3)));
  theta(11) = 1.0 / (1 + exp(-theta(11)));
  return theta;
}

// Main VB algorithm
// [[Rcpp::export]]
Rcpp::List VB_Cars2 (mat x, Rcpp::NumericMatrix lambdaIn, vec hyperParams, int S = 1, int N = 100, int maxIter = 5000,
                    int minIter = 1, double alpha = 0.1, double threshold = 0.01, double thresholdIS = 0){
  
  // convert to Eigen format
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  // Initialise various components for the algorithm
  int T = x.n_rows;
  Matrix<double, Dynamic, 1> gradientJ(24), logPeval(S), tempGrad(24);
  MatrixXd gradientP(24, S);
  
  vec meanGradient(24), Mt(24), Vt(24), LB(maxIter+1);
  mat theta(12, S), epsilon(12, S);
  
  LB.fill(0); Mt.fill(0); Vt.fill(0);
  double logQeval, logJeval;
  // Initial SOBOL QMC numbers. first few numbers in each sequence are too similar so we want to work from 101'th to 100+S'th number in each sequence
  mat sobol = sobol_points(S+100, 12);
  // store quasi random "uniform" [0, 1] numbers that result from shuffling the bits of the sobol numbers
  mat unif;
  // Standard normal density to transform [0,1] numbers to standard normal
  boost::math::normal_distribution<> epsDist(0, 1);
  // Loop control
  int iter = 0;
  double diff = threshold + 1;
  double meanLB = 0;
  double omeanLB;
  
  while(diff > threshold | iter < minIter){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // Reset derivatives to zero
    meanGradient.fill(0);
    Vt.fill(0);
    
    // Importance Sample check
    double checkIS = randu<vec>(1)[0];
    // do not IS on first iteration! 
    if(iter == 1 | checkIS > thresholdIS){
      // Did not importance sample: Draw new epsilon/theta, calculate new dlogp / dtheta and other gradients
      
      // Randomly shuffle the sobol numbers for RQMC
      unif = shuffle(sobol);
      // Create estimates of E(dLB/dlam) and E(dLB/dlam^2)
      for(int s = 0; s < S; ++s){
        // transform uniform numbers to standard normal
        logQeval = 0;
        for(int i = 0; i < 12; ++i){
          // Quantile function is unstable for extreme values of unif
          if(unif(s+100, i) > 0.9995){
            epsilon(i, s) = 3.290527;
          } else if (unif(s+100, i) < 0.0005) {
            epsilon(i, s) = -3.290527;
          } else {
            epsilon(i, s) = quantile(epsDist, unif(s+100, i));
          }
          logQeval += - 0.5 * log(6.28319)  -  pow(epsilon(i, s), 2) / 2;
        }
        // create thetas for importance sampler
        theta.col(s) = constructTheta(epsilon.col(s), lambda);
        // Autodiff
        logP p(x, T, N, epsilon.col(s), hyperParams);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(p, lambda, logPeval(s), tempGrad);
        //gradientP.col(s) = tempGrad;  
        
        logDetJ j(epsilon.col(s));
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(j, lambda, logJeval, gradientJ);
        // update estimates of deriv(logp - logq) and (deriv(logp - logq))^2
        for(int i = 0; i < 24; ++i){
          meanGradient(i) += (gradientP(i, s) + gradientJ(i))/ S;
          //meanGradientSq(i) += pow(S * meanGradient(i), 2) / S;
        }
        LB(iter-1) += (logPeval(s) + logJeval - logQeval)/S;
      }
    } else {
      // Importance Sampling: reuse derivatives and simulated values from most recent non IS iteration
      for(int s = 0; s < S; ++s){
        // Transform theta back to the epsilon that would have generated them with current lambda values
        vec  impliedEpsilon(12), impliedTheta(12);
        logQeval = 0;
        for(int i = 0; i < 12; ++i){
          if(i == 0 | i == 1 | i == 2 | i == 4 | i == 5 | i == 7 | i == 8){
            impliedTheta(i) = log(theta(i, s));
          } else if(i == 2 | i == 3 | i == 11){
            impliedTheta(i) = log(theta(i, s) / (1 - theta(i, s)));
          } else {
            impliedTheta(i) = theta(i, s);
          }
          impliedEpsilon(i) = (impliedTheta(i) - lambda(i)) / lambda(12 + i);
          logQeval += log(pdf(epsDist, impliedEpsilon(i)));
        }
        // Take derivatives of J with these new values
        logDetJ j(impliedEpsilon);
        stan::math::set_zero_all_adjoints();
        stan::math::gradient(j, lambda, logJeval, gradientJ);
        // Calculate importance sample weight as q(impEps)/q(eps), q ~ MVN(0, I)
        double weight = 0;
        for(int i = 0; i < 12; ++i){
          weight += pow(epsilon(i, s), 2) - pow(impliedEpsilon(i), 2);
        }
        weight = exp(weight/2);
        // Calculate derivative estimates
        for(int i = 0; i < 24; ++i){
          meanGradient(i) += (gradientP(i, s)  + gradientJ(i)) * weight / S;
        }
        // Update ELBO
        LB(iter-1) += (logPeval(s) + logJeval - logQeval) / S;
      }
    }
    // adagrad updating rule for lambda values
    for(int i = 0; i < 24; ++i){
      Mt(i) += pow(meanGradient(i), 2);
      if(Mt(i) != 0){
        Vt(i) = pow(Mt(i), -0.5);
      }
      if(iter > 1){
        lambda(i) += alpha * Vt(i) * meanGradient(i);
      }
    }
    // check convergence by looking at the difference of the average of past five LB values and the average of the five before that
    if(iter % 5 == 0){
      omeanLB = meanLB;
      meanLB = 0.2 * (LB(iter-1) + LB(iter-2) + LB(iter-3) + LB(iter-4) + LB(iter-5));
      diff = std::fabs(meanLB - omeanLB);
    } 
    // report progress from time to time
    if(iter % 25 == 0){
      Rcpp::Rcout << "Iteration: " << iter << ", ELBO: " << meanLB << std::endl;
    }
  } // End while loop
  if(iter <= maxIter & iter > 2){ // did converge in time
    // iter goes up by one before checking the maxIter condition, so need to use LB(iter-2)
    Rcpp::Rcout << "Converged after " << iter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  } else if(iter < maxIter){ // reached maxIter, so did not converge
    Rcpp::Rcout << "Warning, failed to converge after " << maxIter << " iterations at ELBO = " << LB(iter-2) << std::endl;
  }
  // Grab elements of lambda as the more meaningful approximation mean vector and upper triangular cholesky factor of the variance matrix
  MatrixXd Mu = lambda.topRows(12);
  MatrixXd Sd = lambda.bottomRows(12);
  return Rcpp::List::create(Rcpp::Named("Mu") = Mu,
                            Rcpp::Named("Sd") = Sd,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
}

