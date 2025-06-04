################################################################################
############################ Extra functions ###################################
################################################################################

################################################################################
lambertW_expArg = function(x) {

  ctrl = x <= 709
  res = double(length = length(x))

  if (sum(ctrl) > 0) res[ctrl] = lamW::lambertW0(exp(x[ctrl]))
  if (sum(!ctrl) > 0) {
    res[!ctrl] = x[!ctrl] - log(x[!ctrl]) + log(x[!ctrl]) / x[!ctrl] +
      (log(x[!ctrl]) * (log(x[!ctrl]) - 2)) / (2 * x[!ctrl]^2) +
      log(x[!ctrl]) * (2 * log(x[!ctrl])^2 - 9 * log(x[!ctrl]) + 6) /
        (6 * x[!ctrl]^3) + log(x[!ctrl]) *
      (
       -12 + 36 * log(x[!ctrl]) - 22 * log(x[!ctrl])^2 +
         3 * log(x[!ctrl])^3) / (12 * x[!ctrl]^4
    )
  }

  return(res)

}
################################################################################



################################################################################
# GompPois_rng = function(T, theta1, theta2, b) {

#   # This function simulates Gompertz model by transformed parameters

#   # T is the length of the chain

#   res = double(T)
#   a = -b * theta1
#   sigmaSq = -b * (2 + b) * theta2

#   res[1] = theta1 + sqrt(theta2) * rnorm(1)
#   for (t in 2:T) {
#     res[t] = a + (1 + b) * res[t - 1] + sqrt(sigmaSq) * rnorm(1)
#   }

#   res = rpois(T, lambda = exp(res))

#   return(res)

# }
################################################################################



################################################################################
# Simulate_Gompertz = function(TT, a, b, sigma2) {

#   # This function simulates Gompertz model by the original parameters

#   theta1 = -a / b
#   theta2 = -sigma2 / (b * (b + 2))

#   Nt = numeric(TT)
#   Nt[1] = exp(rnorm(1, theta1, sqrt(theta2)))
#   Es = rnorm(TT, 0, sqrt(sigma2))

#   for (t in 2:TT) {
#     Nt[t] = Nt[t - 1] * exp(a + b * log(Nt[t - 1]) + Es[(t - 1)])
#   }

#   Nt.obs = numeric(TT)

#   for (t in 1:TT) {
#     Nt.obs[t] = rpois(1, lambda = Nt[t])
#   }

#   return(Nt.obs)

# }
################################################################################



################################################################################
# GompLognorm_rng = function(T, theta1, theta2, b, tauSq) {

#   ### initialize result object
#   res = double(T)

#   ### get the intercept
#   a = -b * theta1

#   ### get the conditional variance
#   sigmaSq = -b * (2 + b) * theta2

#   ### sample the first draw from the stationary distribution
#   res[1] = theta1 + sqrt(theta2) * rnorm(1)

#   ### sample the remaining draws from the conditional distribution
#   for (t in 2:T) {
#     res[t] = a + (1 + b) * res[t - 1] + sqrt(sigmaSq) * rnorm(1)
#   }

#   ## sample observed variables from Poisson
#   res = rpois(T, lambda = exp(res))

#   return(res)

# }
################################################################################



################################################################################
# GompPois_MoM = function(Nstar) {

#   # Gompertz Model with Poisson Sampling Error: Method of Moments
#   #
#   # the parameters of the model are chosen in order to match the empirical
#   #  mean, variance and autocovariance at lag 1.

#   t1 = mean(Nstar)
#   t2prime = var(Nstar) - t1
#   t3 = acf(Nstar, lag.max = 1, plot = FALSE, type = "covariance")$acf[2]

#   theta2hat = log(1 + t2prime / t1^2)
#   theta1hat = log(t1) - 0.5 * theta2hat
#   theta3hat = log(1 + t3 * exp(-2 * theta1hat - theta2hat))

#   b =  theta3hat / theta2hat - 1
#   a = -theta1hat * b
#   sigmaSq = -theta2hat * b * (2 + b)

#   return(list(
#     theta1 = theta1hat,
#     theta2 = theta2hat,
#     a = a,
#     b = b,
#     sigmaSq = sigmaSq
#   ))

# }
################################################################################



################################################################################
# old version of approx.CL2
GompPois_composite_likelihood = function(
  Nstar, theta1, theta2, b, nsim, d = NULL
) {

  # Computes the approximate composite likelihood considering all possible
  # combinations of points with d distance

  # d < length(Nstar)

  curr = 0
  Times = length(Nstar)
  mu = rep(theta1, 2)

  if (is.null(d)) d = Times
  if (d > Times) cat("d must be less than the length of the time series")

  for (t1 in 1:(Times - 1)) {

    limit = min(t1 + d, Times)

    for (t2 in (t1 + 1):limit) {
      B = (1 + b)^abs(
        outer(c(1, abs(t2 - t1 + 1)), c(1, abs(t2 - t1 + 1)), "-")
      )
      cov_matrix = theta2 * B
      N.sample = exp(MASS::mvrnorm(nsim, mu = mu, Sigma = cov_matrix))
      logval = colSums(dpois(Nstar[c(t1, t2)], t(N.sample), log = TRUE))
      curr = curr - log(nsim) + log(sum(exp(logval)))
    }

  }

  return(curr)

}
################################################################################



################################################################################
log_targetPoissonGauss = function(x, n, tauSq, mu) {
  -exp(x) + x * n - 0.5 / tauSq * (x - mu)^2
}
################################################################################



################################################################################
log_proposalGauss = function(x, tauSq, xi) {
  - 0.5 / tauSq * (x - xi)^2
}
################################################################################



################################################################################
quad_form = function(M, x) drop(crossprod(M %*% x, x))
################################################################################



################################################################################
logt_kernel_b = function(w, b, eta2, phi1, phi2, T) {
  r = 1 + b
  sumS_wsq = sum(w[-c(1, T)]^2)
  sumH_Lww = sum(w[-T] * w[-1])
  sumT_wsq = sumS_wsq + sum(w[c(1, T)]^2)
  sumSsq_w = sum(w[-c(1, T)])^2
  sumST_w = sum(w[-c(1, T)]) * sum(w)
  sumTsq_w = sum(w)^2
  rsq = r^2
  quad_eq_no_w = rsq * (eta2 * (T - 2) - 1) - 2 * r * eta2 * (T - 1) +
    eta2 * T + 1
  (1 - 0.5 * T) * log1p(-rsq) - 0.5 * log1p(quad_eq_no_w - 1) -
    (phi1 + 0.5 * T) * log1p(
      (rsq * sumS_wsq - 2 * r * sumH_Lww + sumT_wsq) / (2 * phi2 * (1 - rsq)) -
        (eta2 * (1 - rsq) * (rsq * sumSsq_w - 2 * r * sumST_w + sumTsq_w)) / (
          (1 + r)^2 * 2 * phi2 * quad_eq_no_w
        )
    )
}
################################################################################



################################################################################
logGauss_kernel_theta = function(z, b, theta1, theta2, T = length(z)) {
  colSums(rbind(
    c(dnorm(z[1, ], mean = theta1, sd = sqrt(theta2), log = TRUE)),
    dnorm(
      x = z[2:T, ], mean = -b * theta1 + (1 + b) * z[1:(T - 1), ],
      sd = sqrt(-b * (2 + b) * theta2), log = TRUE
    )
  ))
}
################################################################################



################################################################################
logt_kernel_b2 = function(w, b, eta2, phi1, phi2, T) {
  invB = matrix(0, nrow = T, ncol = T)
  invB[1, 1] = 1
  invB[1, 2] = -(1 + b)
  invB[2, 1] = -(1 + b)
  invB[T, T] = 1
  for (t in 2:(T - 1)) {
    invB[t, t] = 1 + (1 + b)^2
    invB[t, t + 1] = -(1 + b)
    invB[t + 1, t] = -(1 + b)
  }
  invB = invB / (1 - (1 + b)^2)
  logdetB = (T - 1) * log1p(-(1 + b)^2)
  suminvB = sum(invB)
  rowSums_invB = c(rowSums(invB)) # it is B^-1 %*% 1_T

  drop(
    -0.5 * (logdetB + log1p(eta2 * suminvB)) - (phi1 + 0.5 * T) * log1p(
      0.5 * quad_form(invB, w) / phi2 - eta2 / (1 + eta2 * suminvB) * (
        (t(w) %*% rowSums_invB)^2
      ) / (2 * phi2)
    )
  )
}
################################################################################



################################################################################
logt_kernel_b3 = function(w, b, eta2, phi1, phi2, T) {
  res = -0.5 * determinant(
    eta2 * matrix(1, nrow = T, ncol = T) +
      BayesRGMM::AR1.cor(n = T, rho = 1 + b),
    logarithm = TRUE
  )$modulus - 0.5 * (2 * phi1 + T) * log1p(
    0.5 / phi2 * t(w) %*% chol2inv(chol(
      eta2 * matrix(1, nrow = T, ncol = T) +
        BayesRGMM::AR1.cor(n = T, rho = 1 + b)
    )) %*% (w)
  )
  attributes(res) = NULL
  res
}
################################################################################



################################################################################
ellipSlice_b = function(b, sub_nsim, Z, eta1, eta2, phi1, phi2, T) {
  beta = log(- b / (2 + b))
  for (isub_nsim in 1:sub_nsim) {
    V = c(rpostlogiskolmo(1, x = beta))
    strike = logt_kernel_b(
      w = Z - eta1, b = b, eta2 = eta2, phi1 = phi1, phi2 = phi2, T = T
    ) - rexp(n = 1, rate = 1)
    delta = runif(1, min = 0, max = 2 * pi)
    delta_min = delta - 2 * pi
    delta_max = delta
    betaTilde = sqrt(V) * rnorm(1)
    while (TRUE) {
      betaProp = beta * cos(delta) + betaTilde * sin(delta)
      bProp = -2 / (1 + exp(-betaProp))
      loglike_bProp = logt_kernel_b(
        w = Z - eta1, b = bProp, eta2 = eta2, phi1 = phi1, phi2 = phi2, T = T
      )
      if (loglike_bProp >= strike) {
        beta = betaProp
        b = bProp
        break
      } else {
        if (delta < 0) delta_min = delta else delta_max = delta
        delta = runif(1, min = delta_min, max = delta_max)
      }
    }
  }
  return(b)
}
################################################################################


################################################################################
GompPois_mcmc_Z = function(nsim, T, Nstar, theta1, theta2, b) {

  mcmc_sample = matrix(nrow = T, ncol = nsim)
  a = -b * theta1
  sigmaSq = -b * (2 + b) * theta2

  # sample Z from importance density
  IS_obj = c(GompPois_likelihood(
    Nstar, theta1 = theta1, theta2 = theta2, b = b,
    nsim = 1e+4, log = TRUE, details = TRUE
  ))
  Z = IS_obj$Z[, sample(
    c(1:1e+4), size = 1, replace = FALSE, prob = IS_obj$normWeights
  )]

  for (insim in 1:nsim) {
    for (t in 1:T) {

      if (t == 1) {
        mu = (theta1 * sigmaSq + (1 + b) * (Z[2] - a) * theta2) / (
          sigmaSq + (1 + b)^2 * theta2
        )
        tauSq = sigmaSq * theta2 / (sigmaSq + (1 + b)^2 * theta2)
      } else if (t == T) {
        mu = a + (1 + b) * Z[T - 1]
        tauSq = sigmaSq
      } else {
        mu = (a + (1 + b) * (Z[t + 1] + Z[t - 1] - a)) / (1 + (1 + b)^2)
        tauSq = sigmaSq / (1 + (1 + b)^2)
      }

      xi = Nstar[t] * tauSq + mu - lambertW_expArg(
        log(tauSq) + Nstar[t] * tauSq + mu
      )

      logC = -exp(xi) + xi * Nstar[t] - 0.5 / tauSq * (xi - mu)^2

      while (TRUE) {
        x = xi + sqrt(tauSq) * rnorm(100)
        check = -rexp(100) <= -exp(x) + x * Nstar[t] -
          0.5 / tauSq * (x - mu)^2 +
          0.5 / tauSq * (x - xi)^2 -
          logC
        if (sum(check) > 0) {
          Z[t] = x[check][1]
          break
        }
      }

    }
    mcmc_sample[, insim] = Z
  }

  return(mcmc_sample)
}
################################################################################