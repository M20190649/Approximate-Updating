// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
//#include <distributions.h> //figure this out
using namespace Rcpp;
using namespace arma;
using namespace std;


// TODO:
// Find out how to add distribution code as a header to reduce amount of code in one file
// Make phi have a truncated normal marginal so I can fix the x0 prior variance 
// Can repeatedly sample from the MVN until we get a draw in (-1, 1) but I'd like to avoid while loops when the mean of phi can go way outside (-1, 1). Could restrict the mean?
// Find out how to avoid dynamic casting to make it easier to implement new models, but if I'm sticking with DLM's for a while this is not a problem yet

class Distribution{
  // This is the base class for the distributions, and exists so I can have Normal and Inverse Gamma pointers in the same array
public:
  // This function is only here so I can call any of the subclass versions.
  // The code won't compile unless I have a Distribution member function and the subclass member function.
  // There should be a way around this, I can use dynamic cast to tell the compiler I want to explicitly use a normal logdens for example
  // But without the dynamic casting it doesn't know what version of the distribution I want to use, so I need something here
  // It gets overwritten by the subclass version anyway.
  // Mostly I just haven't figured out enough c++ to get this working without dynamic casts
  virtual double logdens(double x) {return -99.99;};
};

class Normal: public Distribution{
  // This just provides methods for the single variate normal distribution
public:
  double mean, variance;
  Normal(double, double);
  void set_values(double m, double v){
    mean = m;
    variance = v;
  }
  void increase_values(double m, double sd){
    mean += m;
    if(sqrt(variance) > -sd){
      variance = pow(sqrt(variance) + sd, 2);
    }
  }
  double logdens (double x) {
    return -0.5 * log(2*3.141593*variance) - pow(x - mean, 2) / (2 * variance);}
  double transform_epsilon(double x){
    return mean + sqrt(variance) * x;
  }
};

Normal::Normal (double m, double v){
  mean = m;
  variance = v;
}

class MultiNormal: public Distribution{
  // This provides methods for the multivariate normal distribution parameterised with the lower triangular cholesky decomposition for variance
public:
  vec mean;
  mat chol;
  MultiNormal(vec, mat);
  void set_values(vec m, mat L){
    mean = m;
    chol = L; 
  }
  void increase_values(vec m, mat L){
    mean += m;
    chol += L;
  }
  double logdens(rowvec x){
    double logdet = 0;
    double N = x.n_elem;
    for(int i = 0; i < N; ++i){
      logdet += -log(abs(chol(i, i)));
    }
    double cons = -N / 2 * log(2.0*3.141593);
    mat sigma = chol * chol.t();
    mat exponent = - ((x - mean.t()) * sigma.i() * (x.t() - mean))/ 2;
    return cons + logdet + exponent(0,0);
  }
  rowvec transform_epsilon(rowvec epsilon){
    int N = epsilon.n_elem;
    rowvec output(N);
    for(int i = 0; i < N; ++i){
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

class InverseGamma: public Distribution{
  // The distribution for 1/x  where x ~ gamma(shape, rate)
  double shape, scale;
public:
  InverseGamma(double, double);
  void set_values(double sh, double sc){
    shape = sh;
    scale = sc;
  }
  void increase_values(double sh, double sc){
    shape += sh;
    scale += sc;
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

cube qSim (MultiNormal* qLatent, int n_samples, int p){
  // We want to return the matrix of epsilons (treat each row as a different realisation of a p-dimensional N(0, I))
  // We also want the transformed matrix of variables
  // Does not apply the exp transform, the algorithm is storing it as log sigma squared. exp is only applied to evaluate log(p(theta, x, y))
  // Maths is written as if we had a distribution for sigma squared instead of log sigma squared.
  cube output(n_samples, p, 2);
  mat epsilon = randn<mat>(n_samples, p);   
  mat theta (n_samples, p);
  for(int i = 0; i < n_samples; ++i){
    // Each row of theta corresponds to the transform of that row of epsilon
    theta.row(i) = qLatent->transform_epsilon(epsilon.row(i));
  }
  output.slice(0) = theta;
  output.slice(1) = epsilon;
  return output;
}

double pdens (Normal* Y[], vec y, Distribution* pLatent[], rowvec latent, int T, int p){
  double priorDens = 0;
  // as the IG logdens function takes log(sigmaSq) as an input we keep the input as latent instead of exp(latent)
  for(int i = 0; i < p; ++i){
    // Prior also includes all of the latent state densities
    priorDens += pLatent[i]->logdens(latent[i]);
  }
  double dataDens = 0;
  for(int t = 0; t < T; ++t){
    dataDens += Y[t]->logdens(y[t]);
  }
  return priorDens + dataDens;
}

void updateP (Normal* Y[], Distribution* pLatent[], rowvec latent, int T){
  // This function updates the parameters of p using a sample of q (when the distributions in p have dependencies)
  // latent[0] = log(sigmaSqY), latent[1] = log(sigmaSqX), latent[2] = phi, latent[3] = mu, latent[4],...,latentT+5] = x_0, ..., x_T
  dynamic_cast<Normal*>(pLatent[4])->set_values(0.0, 4*exp(latent[1])); //X_0 ~ N(0, sigmaSqX / 1-phi^2) 
  // I need an efficient way to sample phi from (-1, 1) as phi outside this range causes negative variance
  // The true value of phi is 0.5, which gives x0 a variance of 4*sigmaSqX, so use it until I can fix the phi sampler 
  for(int t = 0; t < T; ++t){
    // Despite the set_values function in Distribution being overwritten with the version in Normal,
    // unless I explicitly say to use the Normal version (via dynamic cast) this was using the distribution version of set_values
    dynamic_cast<Normal*>(pLatent[t+5])->set_values(latent[2] * latent[t+4], exp(latent[1])); //X_t ~ N(phi*x_t-1, sigmaSqX)
    // Y is a Normal* object not a Distribution* object so it avoids the problem
    Y[t]->set_values(latent[3] + latent[t+5], exp(latent[0])); //Y_t ~N(mu + x_t, sigmaSqY)
  }
} 

double ELBO (Normal* Y[], vec y, Distribution* pLatent[], MultiNormal* qLatent, int T, int p, int n = 250){
  // sample theta - don't need epsilon for this so only keep the theta slice
  mat latentSims = qSim(qLatent, n, p).slice(0);
  double value = 0;
  for(int i = 0; i < n; ++i){
    // Set p parameters based on that sample of q
    updateP (Y, pLatent, latentSims.row(i), T);
    // Evaluate ELBO for that sample
    value += pdens(Y, y, pLatent, latentSims.row(i), T, p) - qLatent->logdens(latentSims.row(i)); 
  }
  // Average over n realisations of theta ~ q
  return value / n;
}

rowvec reparamDeriv (vec y, MultiNormal* qLatent, rowvec latent, rowvec epsilon, int T, int i, bool meanfield = false){
  // Aim is to take the derivative of mu_i and the entire i_th row of L as one row vector
  // We use: dpdf = derivative of p wrt theta
  // dfdm = derivative of theta wrt mu (for the chain rule of the derivative of p wrt mu)
  // dfdL = derivative of theta wrt L
  // dqdm = derivative of q wrt mu (easy functional form, can skip the chain rule)
  // dqdL = derivative of q wrt L
  double dpdf = 0; //Initialise to 0 
  latent[0] = exp(latent[0]); // transform to sigmaSq to make p(y, theta) a little bit easier to deal with
  latent[1] = exp(latent[1]);
  // Next is the derivatives of p(y, theta) wrt theta.
  if(i == 0){ //sigmaSqY
    dpdf = -(T/2 + 2) / latent[i] + 1 / pow(latent[i], 2);
    for(int t = 0; t < T; ++t){
      dpdf += pow(y[t] - latent[3] - latent[t+5], 2) / (2 * pow(latent[i], 2));
    }
  } else if(i == 1){  //sigmaSqX
    dpdf = -(T/2 + 5/2) / latent[i] + 1 / pow(latent[i], 2) + pow(latent[4], 2) / (2 * pow(latent[i], 2));
    for(int t = 5; t < T+5; ++t){
      dpdf += pow(latent[t] - latent[2]*latent[t-1], 2) / (2 * pow(latent[i], 2));
    }
  } else if(i == 2){ //Phi
    for(int t = 5; t < T+5; ++t){
      dpdf += (latent[t]*latent[t-1] - latent[i]*pow(latent[t], 2)) / latent[1];
    }
  } else if(i == 3){ //Mu
    dpdf = - latent[i] / 100;
    for(int t = 0; t < T; ++t){
      dpdf -= (latent[i] + latent[t+5] - y[t]) / latent[0];
    }
  } else if(i == 4){ //X0
    dpdf = (latent[2]*(latent[5]) - (latent[2] + 1)*latent[i]) / latent[1];
  } else if(i == T+5){ //XT
    dpdf = (y[T] - latent[3] - latent[i]) / latent[0] - (latent[i] - latent[2]*latent[i-1]) / latent[1];
  } else { //X1 ... XT-1
    dpdf = (y[i-5] - latent[3] - latent[i]) / latent[0] - (latent[i] - latent[2]*latent[i-1]) / latent[1] +
      latent[2] * (latent[i+1] - latent[2] * latent[i]) / latent[1];
  }
  
  double dfdm = 1; //theta = mu + sum(L*eps) for i > 2. Derivative of 1.
  double dqdm = 0; // dq/dmu is zero except for i = 1, 2
  
  if(meanfield){
    double dfds = epsilon[i];
    double dqds = -1 / qLatent->chol_diag(i);
    if(i < 2){
      dfdm = log(latent[i]);
      dfds *= log(latent[i]);
      dqds -= epsilon[i];
    }
    rowvec derivs(T+6, fill::zeros);
    derivs[0] = dpdf*dfdm - dqdm;
    derivs[1] = dpdf*dfds - dqds;
    return derivs;
  } else {
    vec dfdL(T+5, fill::zeros); 
    for(int j = 0; j <= i; ++j){ //Each L_ij (j <= i) has a derivative of epsilon_j, zero otherwise
      dfdL[j] = epsilon[j];
    }
    // lnq = -( mu1 + mu2 + sum(log(Lii) for i = 1, ..., p) + L11eps1 + L21eps1 + L22eps2) + ln(p(eps))
    vec dqdL(T+5, fill::zeros);
    dqdL[i] = -1/qLatent->chol_diag(i); // dq/DLij is -1/Lij for i = j, otherwise = 0 (except i = 2, j = 1, handled below)
    if(i < 2){
      dqdm = -1; // dq/dmu is -1 for i = 1, 2
      dqdL[0] -= epsilon[0]; // dq/dL11 = -1/L11 - eps1, dq/dL21 = -eps1, so adding eps1 to above
      if(i == 1){ //
        dqdL[1] -= epsilon[1]; // dq/dL22 = -1/L22 - eps2, so adding eps2 to above
      }
      // theta1 = exp(mu1 + L11 eps1), theta2 = exp(mu2 + L21 eps1 + L22 eps2)
      dfdm = latent[i]; // derivative wrt mu is theta
      dfdL *= latent[i]; // derivative wrt Lij is theta_i * epsilon_j for j <= i, dfdL already contains relevant epsilon terms
    }
    double meanDeriv = dpdf*dfdm - dqdm; 
    rowvec LDeriv = dpdf*dfdL.t() - dqdL.t();
    // First element is d(ELBO)/dMu, rest is d(ELBO)/dL_ij for given i and j = 1, ..., T+5
    rowvec derivs(T+6);
    derivs[0] = meanDeriv;
    for(int i = 1; i < T+6; ++i){
      derivs[i] = LDeriv[i-1];
    }
    return derivs;
  }
}

void updateQ (MultiNormal* qLatent, mat updates, int p, bool meanfield){
  // via lambda(t+1) = lambda(t) + Pt dELBO/dlambda
  // updates contains mu derivatives in first column, L derivatives in the rest
  vec mean = updates.col(0);
  mat chol(p, p, fill::zeros);
  if(meanfield){
    chol.diag() = updates.col(1);
  } else {
    chol = updates.submat(0, 1, p-1, p);
  }
  qLatent->increase_values(mean, chol);
}

// [[Rcpp::export]]
Rcpp::List DLM_SGA(vec y, int S, int M, int maxIter, vec initialM, mat initialL, double threshold=0.01, 
                   double alpha=0.01, double beta1=0.9, double beta2=0.999, bool Adam=true, bool meanfield=false){
  // Initialise everything we are going to need
  // T is treated as the total length of y, which is T+S in the written report.
  // So we use data up to y_{T-S} then update using y_{T-S+1:T} instead of data up to y_{T} then update to y_{T+S}
  int T = y.n_elem;
  if(S > T){
    // This would imply that the original MCMC is based on data up to y_{T-S < 0}
    Rcpp::Rcout << "S must be equal to or less than the length of y" << std::endl;
    return Rcpp::List::create();
  }
  mat e(T+5, T+6);
  e.fill(pow(10, -8));
  mat Mt(T+5, T+6, fill::zeros); // doubles as AdaGrad's Gt
  mat Vt(T+5, T+6, fill::zeros); // doubles as AdaGrad's Pt
  mat MtHat;
  mat VtHat;
  mat partials(T+5, T+6); // Partials will be reset to zero at the start of every iteration
  
  // Set up distribution objects. We will use the pointers to these as arguments for the functions SGA_DLM calls.
  Normal* Y[T];
  Distribution* pLatent[T+5]; // pLatent contains Normals and Inverse Gammas, so must use the superclass. This makes dynamic casts required.
  MultiNormal* qLatent;
  pLatent[0] = new InverseGamma(1, 1);
  pLatent[1] = new InverseGamma(1, 1); 
  pLatent[2] = new Normal(0, 0.25); 
  pLatent[3] = new Normal(0, 5); 
  pLatent[4] = new Normal(0, 2);
  for(int i = 0; i < T; ++i){
    Y[i] = new Normal(0, 1);
    pLatent[i+5] = new Normal(0, 1);
  }
  qLatent = new MultiNormal(initialM, initialL);
  
  // Controls the while loop
  int iter = 0;
  // Vector of ELBO value for each iteration
  vec LB(maxIter+1, fill::zeros);
  // ELBO for initial values
  LB[0] = ELBO(Y, y, pLatent, qLatent, T, T+5);
  double diff = threshold + 1;
  
  // Repeat until convergence, make sure it doesnt just stop after one - 20 is pretty arbitary but small enough to not matter
  while((diff > threshold) | (iter < 20)){
    iter += 1;
    if(iter > maxIter){
      break;
    }
    // Reset partial derivatives
    partials.fill(0);
    // Simulate S times from q
    cube Sims = qSim(qLatent, M, T+5);
    // It is easier to split the cube into two matrices now, as calling a row of a matrix of a cube via Sims.slice(i).row(j) doesn't work. 
    // Subcube views can extract the row but it will be stored as a cube object instead of a row vector.
    // This is a bit redundant memory wise but S usually equals 1 so these are not large objects anyway.                                                            
    mat latentSims = Sims.slice(0); 
    mat epsilon = Sims.slice(1);
    // Derivatives are an average of our M samples
    for(int m = 0; m < M; ++m){ 
      updateP(Y, pLatent, latentSims.row(m), T); // Update p distribution parameters based on simulation results
      for(int i = 0; i < 4; ++i){ // Calculate Partial Derivatives wrt global parameters
        partials.row(i) += reparamDeriv(y, qLatent, latentSims.row(m), epsilon.row(m), T, i, meanfield)/M;
      }
      if(S == T){ // If S=T, then we want to estimate all T states AND x0. We can do T-S+5 to start at X_{T-S}, but when
                  // S=T we also want to estimate x0, so we need T-S+4 which just equals 4.
        for(int i = 4; i < T+5; ++i){
          partials.row(i) += reparamDeriv(y, qLatent, latentSims.row(m), epsilon.row(m), T, i, meanfield)/M;
        } 
      } else { // If S<T, then we want to estimate X_{T-S} to X_T. 
        for(int i = T-S+5; i < T+5; ++i){
          partials.row(i) += reparamDeriv(y, qLatent, latentSims.row(m), epsilon.row(m), T, i, meanfield)/M;
        }
      }
    }
    // ADAM Updating, seems to converge faster and more reliably than AdaGrad
    if(Adam){
      Mt = beta1*Mt + (1-beta1)*partials; // Creates biased estimates of first and second moment
      Vt = beta2*Vt + (1-beta2)*pow(partials, 2);
      MtHat = Mt / (1 - pow(beta1, iter)); // Corrects bias
      VtHat = Vt / (1 - pow(beta2, iter));
      updateQ(qLatent, alpha * MtHat / (sqrt(VtHat) + e), T+5, meanfield); // Size of step depends on second moments and alpha
    } else {
    // AdaGrad Updating is available as it was the original method I used in the confirmation etc.
      Mt += pow(partials, 2); // sum of squared derivatives
      for(int i = 0; i < T+5; ++i){
        if(meanfield){
          Vt(i, 0) = alpha * pow(Mt(i, 0), -0.5); // Size of step depends on squared derivatives
          Vt(i, 1) = alpha * pow(Mt(i, 1), -0.5);
        } else {
          for(int j = 0; j <= i+1; ++j){
            Vt(i, j) = alpha * pow(Mt(i, j), -0.5);
          }
        }
       }
      updateQ(qLatent, Vt % partials, T+5, meanfield);
    }
    LB[iter] = ELBO(Y, y, pLatent, qLatent, T, T+5); // Calculate ELBO with new values
    diff = abs(LB[iter] - LB[iter-1]); // Check difference after one extra iteration.
  } // End of while loop
  
  LB = LB.head(iter+1); // The LB vector has length maxIter+1, extract the part we actually used (+1 to include initial LB).
  // Print some useful information to the console
  Rcpp::Rcout << "Number of iterations: " << iter << std::endl;
  Rcpp::Rcout << "Final ELBO: " << LB.tail(1) << std::endl; 
  Rcpp::Rcout << "Final Change in ELBO: " << diff << std::endl;
  if(meanfield){
    // Not interested in the whole L matrix for the meanfield, just the vector of standard deviations
    vec sd(T+5);
    for(int i = 0; i < T+5; ++i){
      sd(i) = abs(qLatent->chol_diag(i));
    }
    return Rcpp::List::create(Rcpp::Named("Mu") = qLatent->mean,
                              Rcpp::Named("Sd") = sd,
                              Rcpp::Named("ELBO") = LB,
                              Rcpp::Named("Iter") = iter);
    
  } else {
    return Rcpp::List::create(Rcpp::Named("Mu") = qLatent->mean,
                            Rcpp::Named("L") = qLatent->chol,
                            Rcpp::Named("ELBO") = LB,
                            Rcpp::Named("Iter") = iter);
  }
}

