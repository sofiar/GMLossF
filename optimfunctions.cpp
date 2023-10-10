#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// Define the function
double subfun2_bis(Rcpp::NumericVector Ns, int M, double theta1, double theta2, double b) {
  int Times = Ns.size();
  arma::vec v1 = {theta2, theta2*(1+b)};
  arma::vec v2 = {theta2*(1+b), theta2};
  arma::mat cov_matrix(2, 2);
  cov_matrix.col(0)=v1;
  cov_matrix.col(1)=v2;
  
  arma::vec mu = {theta1, theta1};
  
  double aa = 0.0;
  
  for (int t = 0; t < Times - 1; ++t) {
    arma::mat N_sample = exp(arma::mvnrnd(mu,cov_matrix,M));
    //Rcpp::Rcout << N_sample.n_cols << "\n";
    arma::uvec positive_rows = arma::find(arma::all(N_sample > 0, 1));
    N_sample = N_sample.rows(positive_rows);
    
    Rcpp::NumericVector dp1(N_sample.n_cols);
    Rcpp::NumericVector dp2(N_sample.n_cols);
    
    for(int i = 0; i < N_sample.n_cols; i++) {
      dp1[i] = R::dpois(Ns[t], N_sample(0, i), 0);
      dp2[i] = R::dpois(Ns[t+1], N_sample(1, i), 0);
    }
    Rcpp::NumericVector dp = dp1 * dp2;
    dp[dp == 0] = 0.0001;
    //Rcpp::Rcout << dp1 << "\n";
    aa += log(Rcpp::mean(dp));
  }
  
  return aa;
  //return dp;

}


// Register the function with R
RCPP_MODULE(mod) {
  Rcpp::function("subfun2_bis_wrap", &subfun2_bis,
                 Rcpp::List::create(Rcpp::_["NsSEXP"], Rcpp::_["MSEXP"], 
                                    Rcpp::_["theta1SEXP"], Rcpp::_["theta2SEXP"], 
                                    Rcpp::_["bSEXP"]));
}