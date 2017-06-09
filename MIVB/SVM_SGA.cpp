// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include "Distributions.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// TODO:
// Make phi have a truncated normal marginal so I can fix the x0 prior variance 
// Can repeatedly sample from the MVN until we get a draw in (-1, 1) but I'd like to avoid while loops when the mean of phi can go way outside (-1, 1). Could restrict the mean?
// Find out how to avoid dynamic casting to make it easier to implement new models, but if I'm sticking with DLM's for a while this is not a problem yet

cube qSim (MultiNormal* qLatent, int n_samples, int p){
  // We want to return the matrix of epsilons (treat each row as a different realisation of a p-dimensional N(0, I))
  // We also want the transformed matrix of variables
  // Does not apply the exp transform, the algorithm is T+Sstoring it as log sigma squared. exp is only applied to evaluate log(p(theta, x, y))
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

double pdens (YStardist* Y[], vec y, Distribution* pLatent[], rowvec latent, int T, int p){
  // Evaluates sum(log(p(yt|xt, theta))) + sum(log(p(xt | xt-1, theta))) + sum(log(p(theta))) 
  double priorDens = 0;
  // as the IG logdens function takes log(sigmaSq) as an input this uses latent instead of exp(latent)
  for(int i = 0; i < p; ++i){
    // Prior also includes all of the latent state densities
    priorDens += pLatent[i]->logdens(latent[i]);
  }
  double dataDens = 0;
  for(int t = 0; t < T; ++t){
    // This is the loglikelihood function for yt | xt, theta
    dataDens += Y[t]->logdens(y[t]);
  }
  return priorDens + dataDens;
}

void updateP (YStardist* Y[], Distribution* pLatent[], rowvec latent, int T){
  // This function updates the parameters of p using a sample of q (when the distributions in p have dependencies)
  // latent[0] = log(sigmaSq), latent[1] = gamma, latent[2] = phi, latent[4],...,latentT+4] = a_0, ..., a_T
  dynamic_cast<Normal*>(pLatent[4])->set_values(latent[1] / (1 - latent[2]), (1-pow(0.9, 2))*exp(latent[0])); //a_0 ~ N(gamma / 1-phi, sigmaSq / 1-phi^2) 
  // I need an efficient way to sample phi from (-1, 1) as phi outside this range causes negative variance
  // True phi is 0.9 so substitute this in for now
  for(int t = 0; t < T; ++t){
    // Despite the set_values function in Distribution being overwritten with the version in Normal,
    // unless I explicitly say to use the Normal version (via dynamic cast) this was using the distribution version of set_values (now removed)
    dynamic_cast<Normal*>(pLatent[t+4])->set_values(latent[1] + latent[2] * latent[t+3], exp(latent[0])); //a_t ~ N(gamma + phi*a_t-1, sigmaSq)
    // Y is a not a Distribution* object so it avoids the problem
    Y[t]->set_values(latent[t+4]); //Y_t - alpha_t ~ log(chi^2_1)
  }
} 

double ELBO (YStardist* Y[], vec y, Distribution* pLatent[], MultiNormal* qLatent, int T, int p, int n = 25){
  // Evaluates the ELBO as an expectation over n simulations of theta & x
  // sample theta - don't need epsilon for this so only keep the theta slice
  mat latentSims = qSim(qLatent, n, p).slice(0);
  double value = 0;
  for(int i = 0; i < n; ++i){
    // Set p parameters based on that sample of q
    updateP (Y, pLatent, latentSims.row(i), T);
    // Evaluate ELBO for that sample
    value += pdens(Y, y, pLatent, latentSims.row(i), T, p) - qLatent->logdens(latentSims.row(i)); 
  }
  // Average over n realisations of q(x, theta)
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
  // Next is the derivatives of p(y, theta) wrt theta.
  if(i == 0){ //sigmaSq
    dpdf = -(0.5*T + 2.5) / latent[i] + 1 / pow(latent[i], 2) + pow(latent[4], 2) * (1 - pow(latent[2], 2)) / (2 * pow(latent[i], 2));
    for(int t = 4; t < T+4; ++t){
      dpdf += pow(latent[t] - latent[1] - latent[2]*latent[t-1], 2) / (2 * pow(latent[i], 2));
    }
  } else if(i == 1){  //Gamma
    dpdf = - latent[i] / 100;
    for(int t = 4; t < T+4; ++t){
      dpdf -= (latent[i] + latent[2]*latent[t-1] - latent[t])  / latent[0];
    }
  } else if(i == 2){ //Phi
    dpdf = latent[i] / (1 - pow(latent[i], 2)) - latent[i] * pow(latent[3], 2) / latent[0];
    for(int t = 4; t < T+4; ++t){
      dpdf += (latent[t]*latent[t-1] - latent[i]*pow(latent[t-1], 2) - latent[t-1]*latent[1]) / latent[0];
    }
  } else if(i == 3){ //Alpha0
    dpdf = - latent[i] * (1 - latent[4] + latent[1]) / latent[0];
  } else if(i == T+3){ //AlphaT
    dpdf = (1 + exp(y[T-1] - latent[i])) * 0.5 - (latent[i] - latent[1] - latent[2]*latent[i-1]) / latent[0];
  } else { //A1 ... AT-1
    dpdf = (1 + exp(y[i-4] - latent[i])) * 0.5 - (latent[i] - latent[1] - latent[2]*latent[i-1]) / latent[0] +
      latent[2] * (latent[i+1] - latent[2] * latent[i] - latent[1]) / latent[0];
  }
  
  double dfdm = 1; //theta = mu + sum(L*eps) for i > 2. Derivative of 1.
  double dqdm = 0; // dq/dmu is zero except for i = 1, 2
  vec dfdL(T+4, fill::zeros);
  if(meanfield){
      dfdL[i] = epsilon[i];
  } else {
      for(int j = 0; j <= i; ++j){ //Each L_ij (j <= i) has a derivative of epsilon_j, zero otherwise
        dfdL[j] = epsilon[j];
      }
    }
  // lnq = -( mu1 + mu2 + sum(log(Lii) for i = 1, ..., p) + L11eps1 + L21eps1 + L22eps2) + ln(p(eps))
  vec dqdL(T+4, fill::zeros);
  dqdL[i] = -1.0 / qLatent->chol_diag(i); // dq/DLij is -1/Lij for i = j, otherwise = 0 
  if(i == 0){
    dfdm = latent[i];
    dfdL[i] *= latent[i];
    dqdm = -1;
    dqdL[i] -= epsilon[i];
  }
  rowvec derivs(T+5, fill::zeros);
  derivs[0] = dpdf*dfdm - dqdm;
  derivs.tail(T+4) = dpdf*dfdL.t() - dqdL.t();
  return derivs;
}


void updateQ (MultiNormal* qLatent, mat updates, int p){
  // via lambda(t+1) = lambda(t) + Pt dELBO/dlambda
  // updates contains mu derivatives in first column, L derivatives in the rest
  vec mean = updates.col(0);
  mat chol = updates.submat(0, 1, p-1, p);
  qLatent->increase_values(mean, chol);
}

// [[Rcpp::export]]
Rcpp::List SVM_SGA(vec y, int S, int M, int maxIter, vec initialM, mat initialL, double threshold=0.01, 
                   double alpha=0.01, double beta1=0.9, double beta2=0.999, bool Adam=true, bool meanfield=false){
  // Initialise everything we are going to need
  // T is treated as the total length of y, which is T+S in the written report.
  // So we use data up to y_{T-S} then update using y_{T-S+1:T} instead of data up to y_{T} then update to y_{T+S}
  // Also using M as a stand in for N - the number of simulations per iteration
  int T = y.n_elem;
  if(S > T){
    // This would imply that the original MCMC is based on data up to y_{T-S < 0}
    Rcpp::Rcout << "Error: S must be less than or equal to the length of y" << std::endl;
    return Rcpp::List::create();
  }
  mat e(T+4, T+5);
  e.fill(pow(10, -8));
  mat Mt(T+4, T+5, fill::zeros); // doubles as AdaGrad's Gt
  mat Vt(T+4, T+5, fill::zeros); // doubles as AdaGrad's Pt
  mat MtHat;
  mat VtHat;
  mat partials(T+4, T+5); // Partials will be reset to zero at the start of every iteration
  
  // Set up distribution objects. We will use the pointers to these as arguments for the functions SGA_DLM calls.
  YStardist* Y[T];
  Distribution* pLatent[T+4]; // pLatent contains Normals and Inverse Gammas, so must use the superclass. This makes dynamic casts required.
  MultiNormal* qLatent;
  pLatent[0] = new InverseGamma(1, 1);
  pLatent[1] = new Normal(0, 10); 
  pLatent[2] = new Uniform(-1, 1); 
  pLatent[3] = new Normal(0, 2);
  for(int i = 0; i < T; ++i){
    Y[i] = new YStardist(1);
    pLatent[i+4] = new Normal(0, 1);
  }
  qLatent = new MultiNormal(initialM, initialL);
  
  // Controls the while loop
  int iter = 0;
  // Vector of ELBO value for each iteration
  vec LB(maxIter+1, fill::zeros);
  // ELBO for initial values
  LB[0] = ELBO(Y, y, pLatent, qLatent, T, T+4);
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
    cube Sims = qSim(qLatent, M, T+4);
    // It is easier to split the cube into two matrices now, as calling a row of a matrix of a cube via Sims.slice(i).row(j) doesn't work. 
    // Subcube views can extract the row but it will be stored as a cube object instead of a row vector.
    // This is a bit redundant memory wise but these are not large objects anyway.                                                            
    mat latentSims = Sims.slice(0); 
    mat epsilon = Sims.slice(1);
    // Derivatives are an average of our M samples
    for(int m = 0; m < M; ++m){ 
      updateP(Y, pLatent, latentSims.row(m), T); // Update p distribution parameters based on simulation results
      for(int i = 0; i < 3; ++i){ // Calculate Partial Derivatives wrt global parameters
        partials.row(i) += reparamDeriv(y, qLatent, latentSims.row(m), epsilon.row(m), T, i, meanfield)/M;
      }
      if(S == T){ // If S=T, then we want to estimate all T states AND x0. Below we do T-S+5 to start at X_{T-S}, but when
        // S=T we also want to estimate x0, so we need T-S+4 which just equals 4.
        for(int i = 3; i < T+4; ++i){
          partials.row(i) += reparamDeriv(y, qLatent, latentSims.row(m), epsilon.row(m), T, i, meanfield)/M;
        } 
      } else { // If S<T, then we want to estimate X_{T-S} to X_T. 
        for(int i = T-S+4; i < T+4; ++i){
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
      updateQ(qLatent, alpha * MtHat / (sqrt(VtHat) + e), T+4); // Size of step depends on second moments and alpha
    } else {
      // AdaGrad Updating is available as it was the original method I used in the confirmation etc. 
      // Partial updating not yet implemented
      Mt += pow(partials, 2); // sum of squared derivatives
      for(int i = 0; i < T+4; ++i){
        if(meanfield){
          Vt(i, 0) = alpha * pow(Mt(i, 0), -0.5); // Size of step depends on squared derivatives
          Vt(i, i+1) = alpha * pow(Mt(i, i+1), -0.5);
        } else {
          for(int j = 0; j <= i+1; ++j){
            Vt(i, j) = alpha * pow(Mt(i, j), -0.5);
          }
        }
      }
      updateQ(qLatent, Vt % partials, T+4);
    }
    LB[iter] = ELBO(Y, y, pLatent, qLatent, T, T+4); // Calculate ELBO with new values
    diff = abs(LB[iter] - LB[iter-1]); // Check difference after one extra iteration.
  } // End of while loop
  if(iter <= maxIter){
    LB = LB.head(iter+1); 
  } else {
    LB = LB.head(iter);
  }
  // The LB vector has length maxIter+1, extract the part we actually used (+1 to include initial LB).
  // If the ELBO didn't converge, the loop will add one to iter before checking iter > MaxIter, so don't do +1 in this case
  // Print some useful information to the console
  Rcpp::Rcout << "Number of iterations: " << iter << std::endl;
  Rcpp::Rcout << "Final ELBO: " << LB.tail(1) << std::endl; 
  Rcpp::Rcout << "Final Change in ELBO: " << diff << std::endl;
  if(meanfield){
    // Not interested in the whole L matrix for the meanfield, just the vector of standard deviations
    vec sd = qLatent->chol.diag();
    for(int i = 0; i < T+4; ++i){
      if(sd[i] < 0){
        sd[i] = -sd[i]; // L values can be negative, sd cannot
      }
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
