#include <Rcpp.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

// Solve w = W(exp(x)) without forming exp(x) explicitly.
// w satisfies: w + log(w) = x, w > 0.
static inline double lambertW_expArg_scalar(const double x) {
  if (!R_finite(x)) return NA_REAL;

  // Initial guess
  double w;
  if (x < 1.0) {
    w = std::exp(x); // small z regime: W(z) ~ z
    if (w < 1e-300) w = 1e-300;
  } else {
    w = x - std::log(x); // asymptotic
    if (w <= 0.0) w = 1e-6;
  }

  // Newton iterations on f(w) = w + log(w) - x
  for (int it = 0; it < 30; ++it) {
    const double f = w + std::log(w) - x;
    const double fp = 1.0 + 1.0 / w;
    double w_new = w - f / fp;

    if (!R_finite(w_new) || w_new <= 0.0) {
      w_new = 0.5 * w;
      if (w_new <= 0.0) w_new = 1e-6;
    }

    if (std::abs(w_new - w) <= 1e-12 * (1.0 + std::abs(w))) {
      w = w_new;
      break;
    }
    w = w_new;
  }

  return w;
}

static inline double safe_exp_neg(const double x) {
  // returns exp(x) with overflow guard; if overflow, return +inf
  //if (x > 709.0) return std::numeric_limits<double>::infinity();
  //if (x < -745.0) return 0.0; // underflow to 0
  return std::exp(x);
}

static inline double log_pi_accept_ratio(
  const double x,
  const double Nstar,
  const double mu,
  const double xi,
  const double tauSq,
  const double logC
) {
  // log ratio used in the accept-reject check:
  // -exp(x) + x*Nstar - 0.5/tauSq * ( (x-mu)^2 - (x-xi)^2 ) - logC
  const double ex = safe_exp_neg(x);
  if (!R_finite(ex)) return R_NegInf;
  const double dxm = x - mu;
  const double dxxi = x - xi;
  return -ex + x * Nstar - 0.5 / tauSq * (dxm * dxm - dxxi * dxxi) - logC;
}

static inline double compute_xi(const double Nstar, const double mu, const double tauSq) {
  // xi = Nstar*tauSq + mu - W(exp(log(tauSq) + Nstar*tauSq + mu))
  const double x = std::log(tauSq) + Nstar * tauSq + mu;
  const double w = lambertW_expArg_scalar(x);
  return Nstar * tauSq + mu - w;
}

static inline double compute_logC(const double xi, const double Nstar, const double mu, const double tauSq) {
  // logC = -exp(xi) + xi*Nstar - 0.5/tauSq * (xi-mu)^2
  const double exi = safe_exp_neg(xi);
  if (!R_finite(exi) || exi == std::numeric_limits<double>::infinity()) return R_NegInf;
  const double d = xi - mu;
  return -exi + xi * Nstar - 0.5 / tauSq * d * d;
}

static inline double sample_ar_one(
  const double Nstar,
  const double mu,
  const double tauSq,
  const int ntrialAR
) {
  //if (!(tauSq > 0.0) || !R_finite(tauSq)) stop("tauSq must be positive and finite");

  const double xi = compute_xi(Nstar, mu, tauSq);
  const double logC = compute_logC(xi, Nstar, mu, tauSq);
  /*if (!R_finite(xi) || !R_finite(logC) || logC == R_NegInf) {
    stop("Numerical issue in xi/logC computation; check parameters (b, theta1, theta2) and data.");
  }*/

  const double s = std::sqrt(tauSq);

  for (;;) {
    for (int j = 0; j < ntrialAR; ++j) {
      const double x = xi + s * R::rnorm(0.0, 1.0);
      if (!R_finite(x)) continue;
      if (x > 709.0) continue; // would overflow exp(x) in acceptance ratio => reject

      const double lr = log_pi_accept_ratio(x, Nstar, mu, xi, tauSq, logC);
      const double lu = std::log(R::runif(0.0, 1.0));
      if (lu <= lr) return x;
    }
  }
}

// [[Rcpp::export]]
NumericVector update_Z_cpp(
  NumericVector Nstar,
  const double theta1,
  const double theta2,
  const double b,
  NumericVector Z,
  const int ntrialAR = 1
) {
  RNGScope scope;

  const int T = Nstar.size();
  //if (Z.size() != T) stop("Nstar and Z must have the same length");
  //if (T < 2) stop("T must be at least 2");
  //if (ntrialAR < 1) stop("ntrialAR must be >= 1");

  const double a = -b * theta1;
  const double sigmaSq = -b * (2.0 + b) * theta2;
  //if (!(sigmaSq > 0.0) || !R_finite(sigmaSq)) stop("sigmaSq must be positive and finite; check b and theta2");

  // Build odd/even indices (0-based)
  std::vector<int> oddPos;
  std::vector<int> evenPos;
  oddPos.reserve((T + 1) / 2);
  evenPos.reserve(T / 2);
  for (int i = 0; i < T; ++i) {
    if (i % 2 == 0) oddPos.push_back(i);  // 1,3,5,... in 1-based
    else evenPos.push_back(i);            // 2,4,6,...
  }

  const int len_odd = (int)oddPos.size();
  const int len_even = (int)evenPos.size();
  const bool isT_even = (T % 2 == 0);

  NumericVector Z_odd(len_odd), Z_even(len_even);
  NumericVector N_odd(len_odd), N_even(len_even);
  for (int i = 0; i < len_odd; ++i) {
    Z_odd[i] = Z[oddPos[i]];
    N_odd[i] = Nstar[oddPos[i]];
  }
  for (int i = 0; i < len_even; ++i) {
    Z_even[i] = Z[evenPos[i]];
    N_even[i] = Nstar[evenPos[i]];
  }

  // -------------------------
  // Update Z at odd positions
  // -------------------------
  NumericVector mu_odd(len_odd);
  NumericVector tauSq_odd(len_odd);

  // First odd element (time 1) depends on first even (time 2)
  //if (len_even < 1) stop("len_even must be >= 1 (T>=2 required)");
  const double onepb = 1.0 + b;
  const double denom1 = sigmaSq + (onepb * onepb) * theta2;
  mu_odd[0] = (theta1 * sigmaSq + (1.0 + b) * (Z_even[0] - a) * theta2) / denom1;
  tauSq_odd[0] = sigmaSq * theta2 / denom1;

  if (len_odd > 1) {
    if (isT_even) {
      // indices 1..len_odd-1 use paired even neighbors
      for (int j = 1; j < len_odd; ++j) {
        // corresponds to Z_even[j-1] and Z_even[j]
        const double num = a + (1.0 + b) * (Z_even[j - 1] + Z_even[j] - a);
        mu_odd[j] = num / (1.0 + (onepb * onepb));
        tauSq_odd[j] = sigmaSq / (1.0 + (onepb * onepb));
      }
    } else {
      // odd has one extra element; last depends only on last even
      for (int j = 1; j < len_odd - 1; ++j) {
        const double num = a + (1.0 + b) * (Z_even[j - 1] + Z_even[j] - a);
        mu_odd[j] = num / (1.0 + (onepb * onepb));
        tauSq_odd[j] = sigmaSq / (1.0 + (onepb * onepb));
      }
      mu_odd[len_odd - 1] = a + (1.0 + b) * Z_even[len_even - 1];
      tauSq_odd[len_odd - 1] = sigmaSq;
    }
  }

  for (int i = 0; i < len_odd; ++i) {
    Z_odd[i] = sample_ar_one(N_odd[i], mu_odd[i], tauSq_odd[i], ntrialAR);
  }

  // --------------------------
  // Update Z at even positions
  // --------------------------
  NumericVector mu_even(len_even);
  NumericVector tauSq_even(len_even);

  if (isT_even) {
    // length len_even: first len_even-1 use pairs in Z_odd; last uses last odd
    for (int j = 0; j < len_even - 1; ++j) {
      const double num = a + (1.0 + b) * (Z_odd[j] + Z_odd[j + 1] - a);
      mu_even[j] = num / (1.0 + (onepb * onepb));
      tauSq_even[j] = sigmaSq / (1.0 + (onepb * onepb));
    }
    mu_even[len_even - 1] = a + (1.0 + b) * Z_odd[len_odd - 1];
    tauSq_even[len_even - 1] = sigmaSq;
  } else {
    // all even indices use paired odds
    for (int j = 0; j < len_even; ++j) {
      const double num = a + (1.0 + b) * (Z_odd[j] + Z_odd[j + 1] - a);
      mu_even[j] = num / (1.0 + (onepb * onepb));
      tauSq_even[j] = sigmaSq / (1.0 + (onepb * onepb));
    }
  }

  for (int i = 0; i < len_even; ++i) {
    Z_even[i] = sample_ar_one(N_even[i], mu_even[i], tauSq_even[i], ntrialAR);
  }

  // Merge back into Z
  for (int i = 0; i < len_odd; ++i) Z[oddPos[i]] = Z_odd[i];
  for (int i = 0; i < len_even; ++i) Z[evenPos[i]] = Z_even[i];

  return Z;
}

