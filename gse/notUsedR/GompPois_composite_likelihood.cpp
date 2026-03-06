#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// Define the function in C++
// [[Rcpp::export]]
double GompPois_composite_likelihood_cpp(NumericVector Nstar, double theta1, 
                                        double theta2, double b, int nsim,
                                        int d) {
  
  // Initialize variables
  int Times = Nstar.size();
  double curr = 0.0;
  VectorXd mu(2);
  d = d + 1; 
  mu << theta1, theta1; // Set mu vector
  
  // Loop through time points
  for (int t1 = 0; t1 < Times - 1; ++t1) {
      VectorXd v(2);
      v << Times, t1+d; 
      int limit = v.minCoeff();
    for (int t2 = t1 + 1; t2 < limit; ++t2) {
      MatrixXd B(2, 2);
      B << pow(1 + b, 0), pow(1 + b, abs(t2 - t1)),
           pow(1 + b, abs(t2 - t1)), pow(1 + b,0); 
            
      MatrixXd cov_matrix = theta2 * B;
      
      // Cholesky decomposition of the covariance matrix
      Eigen::LLT<MatrixXd> cholSolver(cov_matrix);
      MatrixXd L = cholSolver.matrixL();
      
      // Generate samples from multivariate normal distribution
       MatrixXd Z = MatrixXd::Zero(2, nsim); // Initialize Z with zeros
     // Create a Normal distribution with mean 0 and standard deviation sigma
       NumericVector z_0 = Rcpp::rnorm(nsim, 0, 1);
       NumericVector z_1 = Rcpp::rnorm(nsim, 0, 1);
       Z.row(0) = Map<VectorXd>(z_0.begin(), nsim);
       Z.row(1) = Map<VectorXd>(z_1.begin(), nsim);
     
       MatrixXd N_sample = (mu.replicate(1, nsim) + L * Z).array().exp();
      
      // Calculate log likelihood
      NumericVector logval(nsim);
      for (int i = 0; i < nsim; ++i) {
        logval[i] = R::dpois(Nstar[t1], N_sample(0, i), true) + R::dpois(Nstar[t2], N_sample(1, i),true);
      }
      
      curr -= log(nsim) - log(sum(exp(logval)));
        
    }
  }
  
  return curr;
}