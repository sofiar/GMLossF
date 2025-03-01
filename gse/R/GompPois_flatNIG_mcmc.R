#' @title Gompertz Model with Poisson Sampling Error under Flat and
#' Normal-Inverse Gamma Prior: Markov Chain Monte Carlo Algorithm

#' @description
#' Posterior simulation using elliptical slice sampling within Gibbs for
#'  Gompertz model with Poisson sampling error distribution under flat
#'  and normal-inverse gamma priors.

#' @usage
#' GompPois_flatNIG_mcmc = function(nsim, Nstar, phi1, phi2, eta1, eta2,
#'                                  starter = NULL, burn = 1, thin = 1,
#'                                  verbose = +Inf, lastZ = FALSE)

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
#'  b \sim U(-2, 0) \, , \\
#'  \theta_2 \sim Inv.Gamma(\phi_1, \phi_2) \, , \\
#'  \theta_1 \vert \theta_2 \sim N(\eta_1, \eta_2 \theta_2) \, .
#' }
#' Setting \eqn{\beta = \log(-b / (1 + b))}, it is obtained \eqn{\beta \sim
#' Logis(0, 1)}, that is \eqn{\beta \vert V \sim N(0, V)} and \eqn{V} is a
#' logistic Kolmogorov distribution (four times the square of a Kolmogorov
#' distribution). The following updating scheme is used:
#' \enumerate{
#'  \item sample \eqn{V \vert \beta} (posterior of a logistic Kolmogorov),
#'  \item for \eqn{t = 1, 2, ..., T}:
#'        sample \eqn{Z_t \vert N^\star_t, Z_{-t}, b, \theta_1, \theta_2}
#'        (using acceptance-rejection),
#'  \item sample \eqn{\beta \vert Z, V, \theta_1, \theta_2} (using elliptical
#'        slice sampling),
#'  \item sample \eqn{\theta_2, \theta_1 \vert Z, \beta} (inverse gamma and
#'        normal distribution respectively).
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
#' @param phi1 Shape parameter of inverse gamma in NIG prior.
#' @param phi2 Rate parameter of inverse gamma in NIG prior.
#' @param eta1 Mean parameter of normal in NIG prior.
#' @param eta2 Scale parameter of normal in NIG prior.
#' @param starter List for starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @param lastZ Logical, if true the last value of Z is reported.
#' @return Return a matrix containing a posterior sample of the parameters or a
#'         list containing the matrix of posterior sample and vector of last
#'         value of Z.

#' @export
GompPois_flatNIG_mcmc = function(
  nsim, Nstar, phi1, phi2, eta1, eta2, starter = NULL,
  burn = 1, thin = 1, verbose = +Inf, lastZ = FALSE
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
    print(paste("GompPois_flatNIG_mcmc: start time at ", Sys.time(), sep = ""))
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
  beta = log(-b / (2 + b))
  V = pi^2 / 3 # just for initialization, it is not used



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
        bProp = -2 / (1 + exp(-betaProp))
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

      commonFactor = 1 + eta2 * sum(invB) # sum(invB) is 1'_T %*% B^-1 %*% 1_T
      rowSums_invB = c(rowSums(invB)) # it is B^-1 %*% 1_T
      #matrixQuadForm_theta2 = invB - eta2 * tcrossprod(rowSums_invB) /
      #  commonFactor

      theta2 = 1 / rgamma(
        1, shape = phi1 + 0.5 * T, rate = phi2 + 0.5 * quad_form(
          invB - eta2 * tcrossprod(rowSums_invB) / commonFactor, Z - eta1
        )
      )

      theta1 = rnorm(
        1, mean = (eta1 + eta2 * crossprod(rowSums_invB, Z)) / commonFactor,
        sd = sqrt(eta2 * theta2 / commonFactor)
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
    print(paste("GompPois_flatNIG_mcmc: end time at ", Sys.time(), sep = ""))
  }

  ### get results
  res = matrix(
    c(keep_theta1, keep_theta2, keep_b), byrow = TRUE,
    nrow = 3, ncol = nsim, dimnames = list(c("theta1", "theta2", "b"), NULL)
  )

  ### return results
  if (lastZ == FALSE) {
    return(res)
  } else {
    return(list(
      sample = res,
      lastZ = Z
    ))
  }

}