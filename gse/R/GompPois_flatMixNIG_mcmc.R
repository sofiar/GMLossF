#' @title Gompertz Model with Poisson Sampling Error under Flat and Mixture of
#' Normal-Inverse Gamma Prior: Markov Chain Monte Carlo Algorithm

#' @description
#' Posterior simulation using elliptical slice sampling within Gibbs for
#'  Gompertz model with Poisson sampling error distribution under flat
#'  and mixture of normal-inverse gamma priors.

#' @usage
#' GompPois_flatMixNIG_mcmc = function(nsim, Nstar, kappa, psi1, psi2, nu,
#'                                     zeta1, zeta2, starter = NULL, burn = 1,
#'                                     thin = 1, verbose = +Inf)

#' @details
#' This function sample from the posterior distribution of the following model:
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
#'  U \sim Gamma(\psi_1, \psi_2) \, , \\
#'  W \sim Inv.Gamma(\nu/2, \nu/2) \, , \\
#'  \theta_2 \vert U \sim Inv.Gamma(\kappa, U) \, , \\
#'  \theta_1 \vert W, \theta_2 \sim N(\zeta_1, \zeta_2 W \theta_2) \, .
#' }
#' This implies
#' \deqn{
#'  \theta_2 \sim Comp.Gamma(\psi_1, \kappa, \psi_2) \, , \,
#'  \theta_1 \vert \theta_2 \sim Student(\nu, \zeta_1, \zeta_2 \theta_2) \, .
#' }
#' Setting \eqn{\beta = \log(-b / (1 + b))}, it is obtained \eqn{\beta \sim
#' Logis(0, 1)}, that is \eqn{\beta \vert V \sim N(0, V)} and \eqn{V} is a
#' logistic Kolmogorov distribution (four times the square of a Kolmogorov
#' distribution). The following updating scheme is used:
#' \enumerate{
#'  \item sample \eqn{V \vert \beta} (posterior of a logistic Kolmogorov),
#'  \item sample \eqn{U \vert \theta_2} (gamma),
#'  \item sample \eqn{W \vert \theta_1, \theta2} (inverse gamma),
#'  \item for \eqn{t = 1, 2, ..., T}:
#'        sample \eqn{Z_t \vert N^\star_t, Z_{-t}, b, \theta_1, \theta_2}
#'        (using acceptance-rejection),
#'  \item sample \eqn{\beta \vert Z, V, \theta_1, \theta_2} (using elliptical
#'        slice sampling),
#'  \item sample \eqn{\theta_2, \theta_1 \vert Z, \beta, U, W}
#'       (inverse gamma and normal distribution respectively).
#' }
#' See [INSERT REFERENCE] for more details regarding the updating scheme.
#'
#' \code{starter}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{starter == NULL} then:
#'  \enumerate{
#'    \item set \code{b, theta2, theta1} equal to method of moments
#'     estimates,
#'    \item check if \code{b} is in the support and then get \code{beta},
#'    \item sample \code{Z} using importance sampling resampling.
#'  }
#'  \item \code{else} it must be a named list with names equal to the variables
#'    to be initialized, say \code{starter$theta1, starter$theta2, starter$b,
#'    starter$Z}.
#' }
#'
#' Only one value every \code{thin} values is kept in the chain, so the true
#'  number of complete scans will be \code{nsim * thin + burn}. By default
#'  \code{thin = 1}, that is no thinning.
#'
#' The current time and the current number of iteration are printed one every
#' \code{verbose} iterations. Furthermore:
#' \itemize{
#'  \item if \code{verbose == +-Inf} then there is no printing,
#'  \item if \code{verbose != +-Inf} then at least start and end of simulation
#'   are reported.
#' }

#' @param nsim The number of complete scans.
#' @param Nstar The vector of the data.
#' @param kappa Shape parameter at denominator for compounded gamma in mixture
#'              of NIGs prior.
#' @param psi1 Shape parameter at numerator for compounded gamma in mixture of
#'              NIGs prior.
#' @param psi2 Rate parameter for inverse compounded gamma in mixture of NIGs
#'             prior.
#' @param nu Degrees of freedom for Student-t in mixture of NIGs prior.
#' @param zeta1 Location parameter for Student-t in mixture of NIGs prior.
#' @param zeta2 Scale parameter for Student-t in mixture of NIGs prior.
#' @param starter List for starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Return a matrix containing a posterior sample of the parameters.

#' @export
GompPois_flatMixNIG_mcmc = function(
  nsim, Nstar, kappa, psi1, psi2, nu, zeta1, zeta2,
  starter = NULL, burn = 1, thin = 1, verbose = +Inf
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
    print(
      paste("GompPois_flatMixNIG_mcmc: start time at ", Sys.time(), sep = "")
    )
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
    if (starter$b < -0.99) starter$b = -0.99
    # sample Z from importance density
    starter$IS = c(GompPois_likelihood(
      Nstar, theta1 = starter$theta1, theta2 = starter$theta2, b = starter$b,
      nsim = 1e+4, log = TRUE, details = TRUE
    ))
    starter$Z = starter$IS$Z[, sample(
      c(1:1e+4), size = 1, replace = FALSE, prob = starter$IS$normWeights
    )]
  }

  ### initialize parameter values
  theta1 = starter$theta1
  theta2 = starter$theta2
  b = starter$b
  Z = starter$Z

  a = -b * theta1
  sigmaSq = -b * (2 + b) * theta2
  beta = log(-b / (1 + b))
  # set V, U, and W just for initialization, they are not used
  V = pi^2 / 3
  U = 1
  W = 1



  ### sample size
  T = length(Nstar)

  ### initialize inverse of B
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
      #                                update V                                #
      #                                                                        #
      ##########################################################################

      V = c(rpostlogiskolmo(1, x = beta)) #/ sqrt(c)))



      ##########################################################################
      #                                                                        #
      #                                update U                                #
      #                                                                        #
      ##########################################################################

      U = rgamma(1, shape = psi1 + kappa, rate = psi2 + 1 / theta2)



      ##########################################################################
      #                                                                        #
      #                                update W                                #
      #                                                                        #
      ##########################################################################

      W = 1 / rgamma(
        1, shape = 0.5 * (nu + 1),
        rate = 0.5 * (nu + (theta1 - zeta1)^2 / (theta2 * zeta2))
      )



      ##########################################################################
      #                                                                        #
      #                                update Z                                #
      #                                                                        #
      ##########################################################################

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



      ##########################################################################
      #                                                                        #
      #                                update b                                #
      #                                                                        #
      ##########################################################################

      strike = sub_Gomp_loglike_b(b, Z = Z, theta1 = theta1, theta2 = theta2) -
        rexp(1)
      delta = runif(1, min = 0, max = 2 * pi)
      delta_min = delta - 2 * pi
      delta_max = delta
      betaTilde = sqrt(V) * rnorm(1) #sqrt(c * V) * rnorm(1)

      while (TRUE) {
        betaProp = beta * cos(delta) + betaTilde * sin(delta)
        bProp = -1 / (1 + exp(-betaProp))
        loglike_bProp = sub_Gomp_loglike_b(
          bProp, Z = Z, theta1 = theta1, theta2 = theta2
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



      ##########################################################################
      #                                                                        #
      #                        update theta1 and theta2                        #
      #                                                                        #
      ##########################################################################

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

      commonFactor = 1 + zeta2 * W * sum(invB)
      rowSums_invB = c(rowSums(invB)) # it is B^-1 %*% 1_T

      theta2 = 1 / rgamma(
        1, shape = kappa + 0.5 * T, rate = U + 0.5 * quad_form(
          invB - zeta2 * W * tcrossprod(rowSums_invB) / commonFactor, Z - zeta1
        )
      )

      theta1 = rnorm(
        1,
        mean = (zeta1 + zeta2 * W * crossprod(rowSums_invB, Z)) / commonFactor,
        sd = sqrt(zeta2 * W * theta2 / commonFactor)
      )

      a = -b * theta1
      sigmaSq = -b * (2 + b) * theta2



      ##### end of single complete scan
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

    }

  }



  ### print end time if required
  if (!is.infinite(verbose)) {
    print(
      paste("GompPois_flatMixNIG_mcmc: end time at ", Sys.time(), sep = "")
    )
  }

  ### get results
  res = matrix(
    c(keep_theta1, keep_theta2, keep_b), byrow = TRUE,
    nrow = 3, ncol = nsim, dimnames = list(c("theta1", "theta2", "b"), NULL)
  )

  ### return results
  return(res)

}