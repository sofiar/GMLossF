#' @title Gompertz Model with Poisson Sampling Error with Stochastic
#'  Expectation Maximization: Maximum Likelihood Estimator

#' @description
#' Computation of the maximum likelihood estimator for Gompertz model with
#'  Poisson sampling error distribution using stochastic EM algorithm.

#' @usage
#' GompPois_stEM_mle = function(nsim, Nstar, starter = NULL, burn = 1e+3,
#'                              thin = 1, verbose = +Inf)

#' @details
#' This function computes the maximum likelihood estimator for the following
#'  model:
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
#' Thus, a stochastic expectation maximization algorithm is used.
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
GompPois_stEM_mle = function(
  nsim, Nstar, starter = NULL, burn = 1e+3, thin = 1, verbose = +Inf
) {

  #--------------
  sub_Gomp_loglike_b = function(b, Z, theta1, theta2) {
    res = 0
    for (t in 2:T) {
      res = res + dnorm(
        Z[t], mean = -b * theta1 + (1 + b) * Z[t - 1],
        sd = sqrt(-b * (2 + b) * theta2), log = TRUE
      )
    }
    return(res)
  }
  #--------------

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("GompPois_stEM_mle: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn < 0) burn = 0




  ########## initialization

  ##### get starting values

  ### if starting values are not provided then use method of moments
  if (is.null(starter)) {
    starter = GompPois_MoM(Nstar)
    # check if estimated values belong to the support of the parameters
    if (starter$theta2 < 0.01) starter$theta2 = 0.01
    if (starter$b > -0.01) starter$b = -0.01
    if (starter$b < -1.99) starter$b = -1.99
  }

  ### initialize parameter values
  theta1 = starter$theta1
  theta2 = starter$theta2
  b = starter$b

  a = -b * theta1
  sigmaSq = -b * (2 + b) * theta2

  ### sample size
  T = length(Nstar)

  ### initialize inverse of correlation matrix
  invB = matrix(0, nrow = T, ncol = T)



  ### storing vectors
  keep_theta1 = double(nsim)
  keep_theta2 = double(nsim)
  keep_b = double(nsim)



  ########## draw the chain
  for (insim in (1 - burn):nsim) {

    ithin = 0

    ##### iterate up to thinning value
    while (ithin < thin) {

      ##########################################################################
      #                                                                        #
      #                                sample Z                                #
      #                                                                        #
      ##########################################################################

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

      ### get one sample from Z using Gibbs sampler
      for (i_sub_mcmc in 1:10) {

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

      }



      ##########################################################################
      #                                                                        #
      #                     optimize theta1, theta2, and b                     #
      #                                                                        #
      ##########################################################################

      ### get initial value for b using CMLE
      bCMLE = (
        (T - 1) * sum(Z[1:(T - 1)] * Z[2:T]) - sum(Z[1:(T - 1)]) * sum(Z[2:T])
      ) / (
        (T - 1) * sum(Z[1:(T - 1)]^2) - sum(Z[1:(T - 1)])^2
      ) - 1

      ### check if bounds are respected
      if (bCMLE < -1.99) bCMLE = -1.99
      if (bCMLE > -0.01) bCMLE = -0.01

      b = tryCatch(
        expr = {
          optim(
            par = bCMLE,
            fn = function(b) {
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

              rowSums_invB = c(rowSums(invB)) #  B^-1 %*% 1_T
              sum_invB = sum(invB) # 1'_T %*% B^-1 %*% 1_T

              theta1MLE = c(t(rowSums_invB) %*% Z) / sum_invB
              theta2MLE = c(quad_form(invB, Z - theta1MLE)) / T

              return(-sub_Gomp_loglike_b(
                b, Z = Z, theta1 = theta1MLE, theta2 = theta2MLE
              ))
            },
            method = "L-BFGS-B", lower = -1.99, upper = 0.01
          )$par
        },
        error = function(err) return(bCMLE)
      )

      ### get theta1 and theta2
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

      rowSums_invB = c(rowSums(invB)) #  B^-1 %*% 1_T
      sum_invB = sum(invB) # 1'_T %*% B^-1 %*% 1_T

      theta1 = c(t(rowSums_invB) %*% Z) / sum_invB
      theta2 = c(quad_form(invB, Z - theta1)) / T



      ##### end of single iteration
      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }

    ### check if outside of burn-in period
    if (insim > 0) {

      ### keep values
      keep_b[insim] = b
      keep_theta1[insim] = theta1
      keep_theta2[insim] = theta2

      ### print status of the chain
      if (insim %% verbose == 0) {
        print(paste(
          "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
          sep = ""
        ))
      }

    } else {

      ### print status of the chain during burn-in
      if ((insim + burn) %% verbose == 0) {
        print(paste(
          "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
          sep = ""
        ))
      }

    }

  }



  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("GompPois_stEM_mle: end time at ", Sys.time(), sep = ""))
  }

  ### get results
  res = matrix(
    c(keep_theta1, keep_theta2, keep_b), byrow = TRUE,
    nrow = 3, ncol = nsim, dimnames = list(c("theta1", "theta2", "b"), NULL)
  )

  ### return results
  return(res)

}