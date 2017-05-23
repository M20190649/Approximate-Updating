// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// Y - Distribution object for p(yt | xt, theta) for t= 1, ..., T
// Pdist - Distribution object for p(theta), p(x0) and p(xt | xt-1, theta) for t = 1, ..., T
// Qdist - Distribution object for q(theta, x)
// y - data
// sims - simulated value of theta and x from q
// p = T+5, the dimension of the q distribution
// epsilon - simulated standard normal variables that are transformed to theta and x


// All of the p and q densities are stored as distribution objects, which makes it way easier to adapt the code to new models
struct Distribution{
  // This is the base class for the distributions, and exists so I can have Normal and Inverse Gamma pointers in the same array
  // This function is only here so I can call any of the subclass versions.
  // It gets overwritten by the subclass version anyway.
  virtual double logdens(double x) {return -10000;};
};

// All distribution subclass objects have parameters, a constructor, a set parameter value and a logdensity function.
// The multivariate normal also has extra functions to increase parameter values and to transform standard normal 'noise' variables epsilon to theta and x.
// Only the MVN is used for q distributions so it needs the extra functionality.

struct Normal: Distribution{
  // This just provides methods for the single variate normal distribution
  double mean, variance; // Parameters
  Normal(double, double); // Constructor
  void set_values(double m, double v){ 
    mean = m;
    variance = v;
  }
  double logdens (double x) {
    return -0.5 * log(2*3.141593*variance) - pow(x - mean, 2) / (2 * variance);}
};

Normal::Normal (double m, double v){ // Constructor
  mean = m;
  variance = v;
}

struct MultiNormal: Distribution{
  // This provides methods for the multivariate normal distribution parameterised with the lower triangular cholesky decomposition for variance
  vec mean;
  mat chol;
  MultiNormal(vec, mat);
  void set_values(vec m, mat L){
    mean = m;
    chol = L; 
  }
  void increase_values(vec & m, mat & L){
    mean += m;
    chol += L;
  }
  double logdens(rowvec x){
    double logdet = 0;
    double n = x.n_elem;
    for(int i = 0; i < n; ++i){
      logdet += -log(abs(chol(i, i)));
    }
    double cons = -n / 2 * log(2.0*3.141593);
    mat factor = solve(trimatl(chol), x.t() - mean);
    mat exponent = factor.t() * factor / 2;
    return cons + logdet + exponent(0,0);
  }
  rowvec transform_epsilon(rowvec epsilon){
    int n = epsilon.n_elem;
    rowvec output(n);
    for(int i = 0; i < n; ++i){
      output[i] = mean[i];
      for(int j = 0; j <= i; ++j){
        output[i] += chol(i, j)*epsilon[j];
      }
    }
    return output;
  }
  double chol_diag(int i){
    return chol(i,i);
  }
};

MultiNormal::MultiNormal(vec m, mat ch){
  mean = m;
  chol = ch;
}

struct InverseGamma: Distribution{
  // The distribution for 1/x  where x ~ gamma(shape, rate)
  double shape, scale;
  InverseGamma(double, double);
  void set_values(double sh, double sc){
    shape = sh;
    scale = sc;
  }
  double logdens (double x) {
    // the input x should be log(sigmaSq)
    double constant = shape * log(scale) - log(tgamma(shape));
    double kernel = -(shape + 1) * x - scale / exp(x);
    return constant + kernel;
  }
};

InverseGamma::InverseGamma (double sh, double sc){
  shape = sh;
  scale = sc;
}

struct Uniform: Distribution{
  double min, max;
  Uniform(double, double);
  void set_values(double mi, double ma){
    min = mi;
    max = ma;
  }
  double logdens(double x){
    return - log(max - min);
  }
};

Uniform::Uniform (double mi, double ma){
  min = mi;
  max = ma;
}

// Function to simulate n_samples many draws from the p dimensional multivariate normal distribution for q
cube qSim (MultiNormal* qDist, int n_samples, int p){
  // We want to return the matrix of epsilons (treat each row as a different realisation of a p-dimensional MVN(0, I))
  // We also want the transformed matrix of variables
  // Does not apply the exp transform, the code stores these as log sigma squared.
  // However the maths is written for the posterior/approximation to sigma squared instead of log sigma squared.
  cube output(n_samples, p, 2);
  mat epsilon = randn<mat>(n_samples, p);   
  mat theta (n_samples, p);
  for(int i = 0; i < n_samples; ++i){
    // Each row of theta corresponds to the transform of that row of epsilon
    theta.row(i) = qDist->transform_epsilon(epsilon.row(i));
  }
  output.slice(0) = theta;
  output.slice(1) = epsilon;
  return output;
}

// Function to evaluate the density of log(p(y, x, theta)) for a simulated value of theta & x.
double pdens (Normal* Y[], vec y, Distribution* pDist[], rowvec sims, int T, int p){
  // Evaluates sum(log(p(yt|xt, theta))) + sum(log(p(xt | xt-1, theta))) + sum(log(p(theta))) 
  double priorDens = 0;
  // as the IG logdens function takes log(sigmaSq) as an input do not apply exp transform
  for(int i = 0; i < p; ++i){
    // Prior also includes all of the latent state densities p(xt | xt-1, theta)
    priorDens += pDist[i]->logdens(sims[i]);
  }
  double dataDens = 0;
  for(int t = 0; t < T; ++t){
    // This is the loglikelihood function for yt | xt, theta
    dataDens += Y[t]->logdens(y[t]);
  }
  return priorDens + dataDens;
}

// This function updates the parameters of p and Y using simulated values of theta & x.
// This is because the parameters of these distributions typically depend on other variables.
void updateP (Normal* Y[], Distribution* pDist[], rowvec sims, int T){
  // sims[0] = log(sigmaSqY), sims[1] = log(sigmaSqX), sims[2] = phi, sims[3] = gamma, sims[4],...,sims[T+5] = x_0, ..., x_T
  dynamic_cast<Normal*>(pDist[4])->set_values(0.0, 10.25*exp(sims[1])); //X_0 ~ N(0, sigmaSqX / 1-phi^2) 
  // This should used 1/(1-phi^2) instead of 10.25, but I haven't figured out how to sample phi from a truncated marginal distribution
  // Right now phi can be outside (-1, 1), which makes the variance of x0 negative.
  // The true value of phi is 0.95, which gives x0 a variance of 10.25*sigmaSqX, so use it until I can fix the phi sampler 
  for(int t = 0; t < T; ++t){
    // Compiler doesn't know what subclass pDist is so force it to treat it as Normal.
    dynamic_cast<Normal*>(pDist[t+5])->set_values(sims[2] * sims[t+4], exp(sims[1])); //X_t ~ N(phi*x_t-1, sigmaSqX)
    Y[t]->set_values(sims[3] + sims[t+5], exp(sims[0])); //Y_t ~N(gamma + x_t, sigmaSqY)
  }
} 

// Evaluate the ELBO as E_q(p - q) as a sum of (p-q) from n many simulations from q.
double ELBO (Normal* Y[], vec y, Distribution* pDist[], MultiNormal* qDist, int T, int p, int n = 250){
  // sample theta and x only, don't need epsilon so only keep the theta slice from qSim
  mat sims = qSim(qDist, n, p).slice(0);
  double value = 0;
  for(int i = 0; i < n; ++i){
    // Update p parameters based on that sample of q
    updateP (Y, pDist, sims.row(i), T);
    // Evaluate p - q for that sample
    value += pdens(Y, y, pDist, sims.row(i), T, p) - qDist->logdens(sims.row(i)); 
  }
  // Average over n realisations of q(x, theta)
  return value / n;
}

// Take the derivative of the reparameterised ELBO
// d/dlambda (e_eps (logp(f(eps)) - logq(eps) + log|J|))
// Input i is a marker to take the derivative of the i'th element of {theta, x}.
// Derivative of q is different in diagonal and non-diagonal cases.
rowvec reparamDeriv (vec y, MultiNormal* qDist, rowvec sims, rowvec epsilon, int T, int i, bool meanfield){
  // Aim is to take the derivative of mu_i and the entire i_th row of L as one row vector
  // dpdf = derivative of logp wrt theta
  // dfdm = derivative of theta wrt mu (mean vector of Q)
  // dfdL = derivative of theta wrt L (lower triangle of variance matrix)
  // dqdm = derivative of q wrt mu (easy functional form, can skip the chain rule)
  // dqdL = derivative of q wrt L
  double dpdf = 0; //Initialise to 0 
  double sigmaSqY = exp(sims[0]); // transform to sigmaSq to make p(y, theta) a little bit easier to deal with
  double sigmaSqX = exp(sims[1]);
  double phi = sims(2);
  double gamma = sims(3);
  vec x = sims.tail(T+1);
  // Next is the derivatives of p(y, x, theta) wrt theta or x.
  if(i == 0){ //sigmaSqY
    dpdf = -(0.5*T + 2) / sigmaSqY + 1 / pow(sigmaSqY, 2);
    for(int t = 0; t < T; ++t){
      dpdf += pow(y[t] - gamma - x[t], 2) / (2 * pow(sigmaSqY, 2));
    }
  } else if(i == 1){  //sigmaSqX
    dpdf = -(0.5*T + 2.5) / sigmaSqX + 1 / pow(sigmaSqX, 2) + pow(x[0], 2) * (1 - pow(phi, 2)) / (2 * pow(sigmaSqX, 2));
    for(int t = 1; t < T+1; ++t){
      dpdf += pow(x[t] - phi*x[t-1], 2) / (2 * pow(sigmaSqX, 2));
    }
  } else if(i == 2){ //Phi
    dpdf = - phi / (1 - pow(phi, 2)) + phi * pow(x[0], 2) / sigmaSqX;
    for(int t = 1; t < T+1; ++t){
      dpdf += x[t-1]*(x[t] - phi*x[t-1]) / sigmaSqX;
    }
  } else if(i == 3){ //Gamma
    dpdf = - gamma / 100; // As the prior mean of gamma is 0
    for(int t = 0; t < T; ++t){
      dpdf -= (gamma + x[t] - y[t]) / sigmaSqY;
    }
  } else if(i == 4){ //X0
    dpdf = (x[1]*phi -x[0]) / sigmaSqX;
  } else if(i == T+4){ //XT
    dpdf = (y[T] - gamma - x[T]) / sigmaSqY - (x[T] - phi*x[T-1]) / sigmaSqX;
  } else { //X1 ... XT-1
    int t = i - 5;
    dpdf = (y[t] - gamma - x[t]) / sigmaSqY - (x[t] - phi*x[t-1]) / sigmaSqX +
      phi * (x[t+1] - phi * x[t]) / sigmaSqX;
  }
  
  //theta or x = mu + sum(L*eps) for i > 2. Derivative wrt mu is one.
  double dfdm = 1; 
  // lnq = -( mu1 + mu2 + sum(log(Lii) for i = 1, ..., p) + L11eps1 + L21eps1 + L22eps2) + ln(p(eps)) (L21 = 0 if diagonal approximation)
  // derivative wrt mu_i, i > 2 is zero
  double dqdm = 0;  
  vec dfdL(T+5, fill::zeros);
  if(meanfield){
    dfdL[i] = epsilon[i];
  } else {
    for(int j = 0; j <= i; ++j){ //Each L_ij (j <= i) has a derivative of epsilon_j, zero otherwise
      dfdL[j] = epsilon[j];
    }
  }
 
  vec dqdL(T+5, fill::zeros);
  dqdL[i] = -1/qDist->chol_diag(i); // dq/DLij is -1/Lij for i = j, otherwise = 0 
  if(i < 2){ // special cases for sigma squared variables
    dfdm = exp(sims[i]);
    dfdL[i] *= exp(sims[i]);
    dqdm = -1;
    dqdL[i] -= epsilon[i];
    if(!meanfield & i == 1){ // Derivative with respect to L21 (in non-diagonal case)
      dqdL[0] -= epsilon[0];
      dfdL[0] *= exp(sims[i]);
    }
  }
  rowvec derivs(T+6, fill::zeros);
  derivs[0] = dpdf*dfdm - dqdm;
  derivs.tail(T+5) = dpdf*dfdL.t() - dqdL.t();
  return derivs;
}

void updateQ (MultiNormal* qDist, mat updates, int p){
  // via lambda(t+1) = lambda(t) + Pt dELBO/dlambda
  // updates is a p x p+1 matrix
  // updates contains mu derivatives in first column, L derivatives in the rest
  vec mean = updates.col(0);
  mat chol = updates.submat(0, 1, p-1, p); // all columns except column 1, so from (0, 1) to (p-1, p)
  qDist->increase_values(mean, chol);
}

// Stochastic Gradient Ascent for the DLM.
// Inputs: y = Data up to T+S, S = S (only update past S states, update all if S equals length of y), M = simulations per iteration,
// maxIter = maximum number of iterations, initialM = initial value of q mean, initialL = initial value of lower triangle matrix of q variance,
// threshold = convergence criteria, alpha/beta/beta2 = tuning parameters for algorithm, adam = Use adam if true, use adagrad if false 
// Adam works better. meanfield = Diagonal approximation true/false
// [[Rcpp::export]]
Rcpp::List DLM_SGA(vec y, int S, int M, int maxIter, vec initialM, mat initialL, double threshold=0.01, 
                   double alpha=0.01, double beta1=0.9, double beta2=0.999, bool Adam=true, bool meanfield=false){
  // T is treated as the total length of y, which is T+S in the written report.
  // So we use data up to y_{T-S} then update using y_{T-S+1:T} instead of data up to y_{T} then update to y_{T+S}
  int T = y.n_elem;
  if(S > T){
    // This would imply that the original MCMC is based on data up to y_{T-S < 0}
    Rcpp::Rcout << "Error: S must be equal to or less than the length of y" << std::endl;
    return Rcpp::List::create();
  }
  // Initialise everything we are going to need
  mat e(T+5, T+6);
  e.fill(pow(10, -8));
  mat Mt(T+5, T+6, fill::zeros); // doubles as AdaGrad's Gt
  mat Vt(T+5, T+6, fill::zeros); // doubles as AdaGrad's Pt
  mat MtHat;
  mat VtHat;
  mat partials(T+5, T+6); // Partials will be reset to zero at the start of every iteration
  
  // Set up distribution objects. We will use the pointers to these as arguments for the functions SGA_DLM calls.
  Normal* Y[T];
  Distribution* pDist[T+5]; // pDist contains Normals and Inverse Gammas, so must use the superclass. 
  MultiNormal* qDist;
  pDist[0] = new InverseGamma(1, 1); // Sigma^2_y
  pDist[1] = new InverseGamma(1, 1); // Sigma^2_x
  pDist[2] = new Uniform(-1, 1);  // Phi
  pDist[3] = new Normal(0, 100);  // Gamma
  pDist[4] = new Normal(0, 2);  //X0, inputs to all x distributions don't matter and will be changed depending on simulated data
  for(int i = 0; i < T; ++i){
    Y[i] = new Normal(0, 1); //yt
    pDist[i+5] = new Normal(0, 1); //xt
  }
  qDist = new MultiNormal(initialM, initialL);
  
  // Controls the while loop
  int iter = 0;
  // Vector of ELBO value for each iteration
  vec LB(maxIter+1, fill::zeros);
  // ELBO for initial values
  LB[0] = ELBO(Y, y, pDist, qDist, T, T+5);
  double diff = threshold + 1;
  
  // Repeat until convergence, make sure it doesnt just stop after one - 20 is pretty arbitary but small enough to not matter
  while((diff > threshold) | (iter < 20)){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // Reset partial derivatives
    partials.fill(0);
    // Simulate M times from q
    cube simsCube = qSim(qDist, M, T+5);
    // It is easier to split the cube into two matrices now, as calling a row of a matrix of a cube via Sims.slice(i).row(j) doesn't work. 
    // Subcube views can extract the row but it will be stored as a cube object instead of a row vector.
    // This is a bit redundant memory wise but these are not large objects anyway.                                                            
    mat sims = simsCube.slice(0); 
    mat epsilon = simsCube.slice(1);
    // Derivatives are an average of our M samples
    for(int m = 0; m < M; ++m){ 
      // Take the partial derivative with respect to mu_i and i'th row of L in the same function.
      for(int i = 0; i < 4; ++i){ // Calculate Partial Derivatives wrt global parameters
        partials.row(i) += reparamDeriv(y, qDist, sims.row(m), epsilon.row(m), T, i, meanfield)/M;
      }
      if(S == T){ // If S=T, then we want to estimate all T states AND x0. Below we do T-S+5 to start at X_{T-S}, but when
                  // S=T we also want to estimate x0, so we need T-S+4 which just equals 4.
        for(int i = 4; i < T+5; ++i){
          partials.row(i) += reparamDeriv(y, qDist, sims.row(m), epsilon.row(m), T, i, meanfield)/M;
        } 
      } else { // If S<T, then we want to estimate X_{T-S} to X_T. 
        for(int i = T-S+5; i < T+5; ++i){
          partials.row(i) += reparamDeriv(y, qDist, sims.row(m), epsilon.row(m), T, i, meanfield)/M;
        }
      }
    }
    // ADAM Updating, seems to converge faster and more reliably than AdaGrad
    if(Adam){
      Mt = beta1*Mt + (1-beta1)*partials; // Creates biased estimates of first and second moment
      Vt = beta2*Vt + (1-beta2)*pow(partials, 2);
      MtHat = Mt / (1 - pow(beta1, iter)); // Corrects bias
      VtHat = Vt / (1 - pow(beta2, iter));
      updateQ(qDist, alpha * MtHat / (sqrt(VtHat) + e), T+5); // Size of step depends on second moments and alpha
    } else {
    // AdaGrad Updating is available as it was the original method I used in the confirmation etc. 
    // Only works when S = length of y, adam is better so I never implemented updating only part of the x vector
      Mt += pow(partials, 2); // sum of squared derivatives
      for(int i = 0; i < T+5; ++i){
        if(meanfield){
          Vt(i, 0) = alpha * pow(Mt(i, 0), -0.5); // Size of step depends on squared derivatives
          Vt(i, i+1) = alpha * pow(Mt(i, i+1), -0.5); // Diagonal of L matrix only for meanfield
        } else {
          for(int j = 0; j <= i+1; ++j){
            Vt(i, j) = alpha * pow(Mt(i, j), -0.5); // Each row/column of L for non-meanfield
          }
        }
       }
      updateQ(qDist, Vt % partials, T+5); // Increase values of QDist by stepsize * partial derivatives
    }
    LB[iter] = ELBO(Y, y, pDist, qDist, T, T+5); // Calculate ELBO with new values
    diff = abs(LB[iter] - LB[iter-1]); // Check difference after one extra iteration.
  } // End of while loop
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  }
  // If the ELBO converged, the LB vector has length maxIter+1, extract the part we actually used (+1 to include initial LB).
  // Print some useful information to the console during testing
  //Rcpp::Rcout << "Number of iterations: " << iter << std::endl;
  //Rcpp::Rcout << "Final ELBO: " << LB.tail(1) << std::endl; 
  //Rcpp::Rcout << "Final Change in ELBO: " << diff << std::endl;
  if(meanfield){
    // Not interested in the whole L matrix for the meanfield, just the vector of standard deviations
    vec sd = qDist->chol.diag();
    for(int i = 0; i < T+5; ++i){
      if(sd[i] < 0){
        sd[i] = -sd[i]; // L values can be negative, sd cannot
      }
    }
    return Rcpp::List::create(Rcpp::Named("Mu") = qDist->mean,
                              Rcpp::Named("Sd") = sd,
                              Rcpp::Named("ELBO") = LB,
                              Rcpp::Named("Iter") = iter);
    
  } else {
    return Rcpp::List::create(Rcpp::Named("Mu") = qDist->mean,
                            Rcpp::Named("L") = qDist->chol,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
  }
}

