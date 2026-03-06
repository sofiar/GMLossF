#include <Rcpp.h>
#include <Rmath.h>  // unif_rand(), norm_rand(), rgamma()
#include <cmath>
#include <algorithm>
#include <vector>
#include "LBFGSB.h"
#include <Eigen/Dense>

using namespace Rcpp;

static inline double logt_kernel_b_single(
    const double* w, int n,
    double b, double eta2, double phi1, double phi2, int T
){
  double sumS_wsq = 0.0; for (int i = 1; i < n - 1; ++i) sumS_wsq += w[i] * w[i];
  double sumH_Lww = 0.0; for (int i = 0; i < n - 1;   ++i) sumH_Lww += w[i] * w[i + 1];
  double sumT_wsq = sumS_wsq + w[0]*w[0] + w[n - 1]*w[n - 1];

  double sumS_w = 0.0;   for (int i = 1; i < n - 1; ++i) sumS_w += w[i];
  double sumT_w = sumS_w + w[0] + w[n - 1];

  const double sumSsq_w = sumS_w * sumS_w;
  const double sumST_w  = sumS_w * sumT_w;
  const double sumTsq_w = sumT_w * sumT_w;

  const double Td  = static_cast<double>(T);
  const double r   = 1.0 + b;
  const double rsq = r * r;
  const double one_minus_rsq = 1.0 - rsq;

  if (one_minus_rsq <= 0.0) return R_NegInf;

  const double quad_eq =
      rsq*(eta2*(T - 2) - 1.0) - 2.0*r*eta2*(T - 1) + eta2*Td + 1.0;

  const double term1 = (1.0 - 0.5*Td) * std::log(one_minus_rsq);
  const double term2 = -0.5 * std::log(quad_eq);

  const double numer_quad   = rsq*sumS_wsq - 2.0*r*sumH_Lww + sumT_wsq;
  const double one_plus_rsq = (1.0 + r) * (1.0 + r);

  const double numer_eta = eta2 * one_minus_rsq *
      (rsq*sumSsq_w - 2.0*r*sumST_w + sumTsq_w);

  const double inner =
      numer_quad/(2.0*phi2*one_minus_rsq)
    - numer_eta /(2.0*phi2*quad_eq*one_plus_rsq);

  const double term3 = -(phi1 + 0.5*Td) * std::log1p(inner);

  return term1 + term2 + term3;
}

 // [[Rcpp::export]]
NumericVector update_pars_cpp(
    int T,
    NumericVector b_space,
    double eta1, double eta2,
    double phi1, double phi2,
    NumericVector Z
){
  const int n  = Z.size();
  const int nb = b_space.size();

  // -------------------------------------------------------------------------
  // (1) Build working vector w = Z - eta1
  // -------------------------------------------------------------------------
  std::vector<double> w(n);
  for (int i = 0; i < n; ++i) w[i] = Z[i] - eta1;

  // -------------------------------------------------------------------------
  // (2) Compute logM over the discrete grid b_space (grid search for proposal)
  //     This finds the maximizer of the log-kernel over the coarse grid.
  // -------------------------------------------------------------------------
  double logM = R_NegInf;
  int best_idx = 0;
  for (int j = 0; j < nb; ++j) {
    const double val = logt_kernel_b_single(w.data(), n, b_space[j], eta2, phi1, phi2, T);
    if (val > logM) { logM = val; best_idx = j; }
  }

// -------------------------------------------------------------------------
// (3) Local refinement around the best grid point using L-BFGS-B (LBFGSpp)
//     This tightens the envelope for the accept–reject sampler.
// -------------------------------------------------------------------------
{
    using namespace LBFGSpp;

    const int dim = 1;

    // Initial guess from the discrete grid
    Eigen::VectorXd x(dim);
    x[0] = b_space[best_idx];

    // Bounds for b
    Eigen::VectorXd lb(dim), ub(dim);
    lb[0] = -1.999;
    ub[0] = -0.001;

    // Define the functor for LBFGSB: LBFGSpp performs MINIMIZATION,
    // so minimize -log target.
    struct Functor {
        const double* w;
        int n, T;
        double eta2, phi1, phi2;

        double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) const {
            const double b = x[0];

            // objective f(b) = - log likelihood kernel
            const double f0 = -logt_kernel_b_single(w, n, b, eta2, phi1, phi2, T);

            // numerical gradient
            const double eps = 1e-6;
            const double b1 = std::min(-0.001, b + eps);
            const double b2 = std::max(-1.999, b - eps);

            const double f1 = -logt_kernel_b_single(w, n, b1, eta2, phi1, phi2, T);
            const double f2 = -logt_kernel_b_single(w, n, b2, eta2, phi1, phi2, T);

            grad.resize(1);
            grad[0] = (f1 - f2) / (2.0 * eps);

            return f0;
        }
    };

    Functor fun = { w.data(), n, T, eta2, phi1, phi2 };

    // Solver parameters
    LBFGSBParam<double> param;
    param.max_iterations = 50;
    param.epsilon = 1e-8;

    LBFGSBSolver<double> solver(param);

    double fx;  // minimized function value

    // Run L-BFGS-B
    solver.minimize(fun, x, fx, lb, ub);

    // Extract optimized b and the LOG target maximum
    const double b_opt     = x[0];
    const double logf_opt  = -fx;

    // Update envelope constant used in accept–reject
    logM = logf_opt;
}

  // -------------------------------------------------------------------------
  // (4) Accept–Reject sampling step for b using envelope defined via logM
  //     Proposal is Uniform(-2, 0); accept with probability exp(lr)
  // -------------------------------------------------------------------------
  double b = b_space[best_idx];
  const int batch = 1;

  for (int it = 0; it < 1000000000; ++it) {
    for (int k = 0; k < batch; ++k) {

      // propose x ~ Uniform(-2,0)
      const double x  = -2.0 + unif_rand()*2.0;

      // log acceptance ratio
      const double lr = logt_kernel_b_single(w.data(), n, x, eta2, phi1, phi2, T) - logM;

      // compare with log(U)
      const double lu = std::log(unif_rand());
      if (lu <= lr) {
        b = x;
        goto B_DONE;
      }
    }
  }
B_DONE:

  // -------------------------------------------------------------------------
  // (5) Update theta2 (Inverse-Gamma full conditional)
  // -------------------------------------------------------------------------
  const double r = 1.0 + b;
  const double rsq = r*r;
  const double one_minus_rsq = 1.0 - rsq;

  double sumS_zsq = 0.0; for (int i = 1; i < n-1; ++i) sumS_zsq += w[i]*w[i];
  double sumH_Lzz = 0.0; for (int i = 0; i < n-1; ++i) sumH_Lzz += w[i]*w[i+1];
  double sumT_zsq = 0.0; for (int i = 0; i < n;   ++i) sumT_zsq += w[i]*w[i];

  double sumS_z   = 0.0; for (int i = 1; i < n-1; ++i) sumS_z   += w[i];
  double sumT_z   = 0.0; for (int i = 0; i < n;   ++i) sumT_z   += w[i];

  const double one_p_eta2_suminvB =
      (rsq*(eta2*(T - 2) - 1.0)
     - 2.0*r*eta2*(T - 1)
     + eta2*T + 1.0) / one_minus_rsq;

  const double denom = (2.0 + b) * (2.0 + b);

  const double quad_form_z =
      (rsq*sumS_zsq - 2.0*r*sumH_Lzz + sumT_zsq)/one_minus_rsq
    - (eta2/one_p_eta2_suminvB) *
      (rsq*sumS_z*sumS_z - 2.0*r*sumS_z*sumT_z + sumT_z*sumT_z) / denom;

  const double shape = phi1 + 0.5 * T;
  const double rate  = phi2 + 0.5 * quad_form_z;

  // theta2 update via inverse-gamma
  const double theta2 = 1.0 / R::rgamma(shape, 1.0 / rate);

  // -------------------------------------------------------------------------
  // (6) Update theta1 (Normal full conditional)
  // -------------------------------------------------------------------------
  double sumT_Z = 0.0; for (int i = 0; i < n; ++i)     sumT_Z += Z[i];
  double sumS_Z = 0.0; for (int i = 1; i < n - 1; ++i) sumS_Z += Z[i];

  const double rowsumsInvB_Z = (sumT_Z - r*sumS_Z) / (2.0 + b);
  const double mean_t1 = (eta2*rowsumsInvB_Z + eta1) / one_p_eta2_suminvB;
  const double sd_t1   = std::sqrt(eta2 * theta2 / one_p_eta2_suminvB);

  // theta1 update
  const double theta1 = mean_t1 + sd_t1 * norm_rand();

  // -------------------------------------------------------------------------
  // (7) Return results in the format expected by mcmc_loop_cpp
  // -------------------------------------------------------------------------
  return NumericVector::create(
    b, theta1, theta2
  );
}