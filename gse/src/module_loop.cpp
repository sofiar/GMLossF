#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mcmc_loop_cpp(double burn_d, double nsim_d,
                   double thin_d, double verbose_d,
                   NumericVector Nstar, double T_d,
                   NumericVector b_space,
                   double eta1, double eta2,
                   double phi1, double phi2,
                   double theta1_init, double theta2_init,
                   double b_init, NumericVector Z_init,
                   Function update_Z_fn, Function update_pars_fn) {

  R_xlen_t burn    = (R_xlen_t)burn_d;
  R_xlen_t nsim    = (R_xlen_t)nsim_d;
  R_xlen_t thin    = (R_xlen_t)thin_d;
  R_xlen_t verbose = (R_xlen_t)verbose_d;
  int T = (int)T_d;

  double theta1 = theta1_init;
  double theta2 = theta2_init;
  double b = b_init;
  NumericVector Z = clone(Z_init);

  NumericVector keep_b(nsim);
  NumericVector keep_theta1(nsim);
  NumericVector keep_theta2(nsim);

  Function sys_time("Sys.time");
  Function paste0("paste0");

  for (R_xlen_t insim = 1 - burn; insim <= nsim; insim++) {

    R_xlen_t ithin = 0;

    while (ithin < thin) {

      Z = as<NumericVector>(
        update_Z_fn(Nstar, theta1, theta2, b, Z, 1)
      );

      NumericVector params = as<NumericVector>(
        update_pars_fn(T, b_space, eta1, eta2, phi1, phi2, Z)
      );
      b      = params[0];
      theta1 = params[1];
      theta2 = params[2];

      if (insim > 0) ithin++; else ithin = thin;
    }

    if (insim > 0) {
      keep_b[insim - 1]      = b;
      keep_theta1[insim - 1] = theta1;
      keep_theta2[insim - 1] = theta2;

      if ((insim % verbose) == 0) {
        Rcpp::Rcout << "iteration " << insim << " of " << nsim
                    << " completed at time " << as<std::string>(paste0(sys_time()))
                    << std::endl;
      }
    } else {
      if (((insim + burn) % verbose) == 0) {
        Rcpp::Rcout << "burn-in " << (insim + burn) << " of " << burn
                    << " completed at time " << as<std::string>(paste0(sys_time()))
                    << std::endl;
      }
    }

    Rcpp::checkUserInterrupt();
  }

  R_xlen_t nZ = Z.size();
  NumericVector out(3 * nsim + 3 + nZ);
  R_xlen_t pos = 0;
  for (R_xlen_t i = 0; i < nsim; i++) out[pos++] = keep_b[i];
  for (R_xlen_t i = 0; i < nsim; i++) out[pos++] = keep_theta1[i];
  for (R_xlen_t i = 0; i < nsim; i++) out[pos++] = keep_theta2[i];
  out[pos++] = b;
  out[pos++] = theta1;
  out[pos++] = theta2;
  for (R_xlen_t i = 0; i < nZ; i++) out[pos++] = Z[i];
  return out;
}
