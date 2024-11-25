#' @title Gompertz Model with Poisson Sampling Error with Stochastic
#'  Expectation Maximization: Empirical Bayes

#' @description
#' Calibration of hyper-parameters of the NIG prior for Gompertz model with
#'  Poisson sampling error distribution using stochastic EM algorithm.

#' @usage
#' GompPois_stEM_eb = function(nsim, Nstar, starter = NULL, burn = 1e+3,
#'                             thin = 1, verbose = +Inf)

#' @details
#' This function calibrates the hyper-parameters of the NIG prior for the
#'  following model:
#'  \deqn{
#'  N^\star_t | N_t \sim \mathrm{Poisson}(N_t) \, , \\
#'  N_t = \exp(Z_t) \, , \\
#'  Z_t = a + (1 + b) Z_{t - 1} + \varepsilon_t \, , \\
#'  \varepsilon_t \overset{i.i.d.}{\sim} N(0, \sigma^2) \, .
#'  }
#'
#' A different parametrization is used, say:
#'  \deqn{
#'  \theta_1 = -\frac{a}{b} \, , \, \theta_2 = -\frac{\sigma^2}{b(2 + b)}
#'  \, ,
#'  }
#' which implies
#' \deqn{
#'  Z_1, Z_2, ..., Z_T \sim N_T ( \theta_1 1_T, \theta_2 \bar{\Sigma})  \, , \\
#'  \bar{\Sigma}_{ij} = (1 + b ) ^ {\vert i - j \vert} \, .
#' }
#' The priors are:
#' \deqn{
#'  b \sim U(-1, 0) \, , \\
#'  \theta_2 \sim Inv.Gamma(\phi_1, \phi_2) \, , \\
#'  \theta_1 \vert \theta_2 \sim N(\eta_1, \eta_2 \theta_2) \, .
#' }
#' Setting \eqn{\beta = \log(-b / (1 + b))}, it is obtained \eqn{\beta \sim
#' Logis(0, 1)}, that is \eqn{\beta \vert V \sim N(0, V)} and \eqn{V} is a
#' logistic Kolmogorov distribution (four times the square of a Kolmogorov
#' distribution).
#' Thus, a stochastic expectation maximization algorithm is used to calibrate
#' \eqn{\phi_1, \phi_2, \eta_1, \eta_2}.
#' See [INSERT REFERENCE] for more details.

#' @references
#' [INSERT REFERENCE]

#' @param nsim The number of simulations for the Markov chain.
#' @param Nstar The vector of the data.
#' @param starter List for starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Return a matrix containing a sequence of optimizers of the
#'         parameters.

#' @export
GompPois_stEM_eb = function(
  nsim, Nstar, starter = NULL, burn = 1e+3, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("GompPois_stEM_eb: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn < 0) burn = 0




  ########## initialization

  ### sample size
  T = length(Nstar)

  ### get method of moments estimate of parameters
  res_MoM = GompPois_MoM(Nstar)
  # check if estimated values belong to the support of the parameters
  if (res_MoM$theta2 < 0.01) starter$theta2 = 0.01
  if (res_MoM$b > -0.01) starter$b = -0.01
  if (res_MoM$b < -0.99) starter$b = -0.99



  ##### get starting values

  ### if starting values are not provided then use method of moments
  if (is.null(starter)) {

    # get starting points of parameters from MoM estimates
    starter$theta1 = res_MoM$theta1
    starter$theta2 = res_MoM$theta2
    starter$b = res_MoM$b

    # get starting points of hyper-parameters from MoM estimates ### TO BE CHANGED WITH APPROPRIATE FORMULA FOR CORRELATED DATA
    starter$eta1 = starter$theta1
    starter$eta2 = 1 / T
    starter$phi1 = 0.5 * (T - 1)
    starter$phi2 = 0.5 * (T - 1) / starter$theta2
  }

  ### initialize parameter values
  theta1 = starter$theta1
  theta2 = starter$theta2
  b = starter$b

  phi1 = starter$phi1
  phi2 = starter$phi2
  eta1 = starter$eta1
  eta2 = starter$eta2

  a = -b * theta1
  sigmaSq = -b * (2 + b) * theta2

  ### initialize inverse of correlation matrix
  invB = matrix(0, nrow = T, ncol = T)

  ### initialize sub_mcmc results object
  sub_mcmc = list(
    sample = matrix(
      c(theta1, theta2, b), nrow = 3,
      dimnames = list(c("theta1", "theta2", "b"), NULL)
    ),
    lastZ = Nstar # just for initialization, it is not used
  )



  ### storing vectors
  keep_eta1 = double(nsim)
  keep_eta2 = double(nsim)
  keep_phi1 = double(nsim)
  keep_phi2 = double(nsim)



  ########## draw the chain
  for (insim in (1 - burn):nsim) {

    ithin = 0

    ##### iterate up to thinning value
    while (ithin < thin) {

      ##########################################################################
      #                                                                        #
      #                             sample Z and b                             #
      #                                                                        #
      ##########################################################################

      ### get new starting point for theta1, theta2, and b
      theta1 = 0.25 * mean(sub_mcmc$sample["theta1", ]) + res_MoM$theta1 * 0.75
      theta2 = 0.25 * mean(sub_mcmc$sample["theta2", ]) + res_MoM$theta2 * 0.75
      b = 0.25 * mean(sub_mcmc$sample["b", ]) + res_MoM$b * 0.75

      ### sample starting point of Z using importance sampling resampling
      # get weighted sample
      Zis = GompPois_likelihood(
        Nstar, theta1 = theta1, theta2 = theta2, b = b,
        nsim = 1e+4, log = TRUE, details = TRUE
      )
      # get one draw from the weighted sample
      Z = Zis$Z[, sample(
        c(1:1e+4), size = 1, replace = FALSE, prob = Zis$normWeights
      )]

      ### get one sample from Z
      sub_mcmc = GompPois_flatNIG_mcmc(
        100, Nstar = Nstar, phi1 = phi1, phi2 = phi2, eta1 = eta1, eta2 = eta2,
        starter = list(theta1 = theta1, theta2 = theta2, b = b, Z = Z),
        lastZ = TRUE
      )

      ### get last value of Z
      Z = sub_mcmc$lastZ

      ### get last value of b
      b = tail(sub_mcmc$sample["b", ], 1)

      ### get inverse of B
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

      ### get constant values with respect to inverse of B
      sum_invB = c(sum(invB)) # 1'_T %*% B^-1 %*% 1_T
      invB_ones_sq = tcrossprod(rowSums(invB)) # B^-1 %*% 1_T %*% 1'_T %*% B^-1



      ##########################################################################
      #                                                                        #
      #                  optimize phi1, phi2, eta1, and eta2                   #
      #                                                                        #
      ##########################################################################

      ### optimize
      opt_obj = tryCatch(
        expr = {
          optim(
            par = c(phi1, phi2, eta1, eta2),
            fn = function(xfun) {
              - (
                lgamma(xfun[1] + 0.5 * T) - lgamma(xfun[1]) -
                  0.5 * T * log(xfun[2]) - 0.5 * log(1 + xfun[4] * sum_invB) -
                  (xfun[1] + 0.5 * T) * log(1 + 0.5 * quad_form(
                    invB - xfun[4] * invB_ones_sq / (1 + xfun[4] * sum_invB),
                    Z - xfun[3]
                  ) / xfun[2])
              )
            },
            method = "L-BFGS-B",
            lower = c(0.01, 0.01, -Inf, 0.01), upper = rep(+Inf, 4)
          )
        },
        error = function(err) return(bCMLE)
      )

      ### get values
      phi1 = opt_obj$par[1]
      phi2 = opt_obj$par[2]
      eta1 = opt_obj$par[3]
      eta2 = opt_obj$par[4]



      ##### end of single iteration
      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }

    ### check if outside of burn-in period
    if (insim > 0) {

      ### keep values
      keep_eta1[insim] = eta1
      keep_eta2[insim] = eta2
      keep_phi1[insim] = phi1
      keep_phi2[insim] = phi2

      ### print status of the chain
      if (insim %% verbose == 0) {
        print(paste(
          "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
          sep = ""
        ))
      }

    }

  }



  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("GompPois_stEM_eb: end time at ", Sys.time(), sep = ""))
  }

  ### get results
  res = matrix(
    c(keep_eta1, keep_eta2, keep_phi1, keep_phi2), byrow = TRUE,
    nrow = 4, ncol = nsim, dimnames = list(
      c("eta1", "eta2", "phi1", "phi2"), NULL
    )
  )

  ### return results
  return(res)

}