#' @title Gompertz Model with Poisson Sampling Error under Flat and
#' Normal-Inverse Gamma Prior: Markov Chain Monte Carlo Algorithm

#' @description
#' Posterior simulation using elliptical slice sampling within Gibbs for
#'  Gompertz model with Poisson sampling error distribution under flat
#'  and normal-inverse gamma priors.

#' @usage
#' GompPois_discrFlatNIG_mcmc = function(nsim, Nstar, phi1, phi2, eta1, eta2,
#'                                  starter = NULL, burn = 0, thin = 1,
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
GompPois_discrFlatNIG_mcmc = function(
  nsim, Nstar, phi1, phi2, eta1, eta2, starter = NULL,
  burn = 0, thin = 1, verbose = +Inf, lastZ = FALSE
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "GompPois_discrFlatNIG_mcmc: start time at ", Sys.time(), sep = ""
    ))
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



  ### discrete prior for b
  b_space = c(1:199) / 100 - 2
  b_weights = double(length = length(b_space))

  ### sample size
  T = length(Nstar)

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
      #                                update Z                                #
      #                                                                        #
      ##########################################################################

      ### get values of alternative parametrization
      a = -b * theta1
      sigmaSq = -b * (2 + b) * theta2

      ### sample sequentially new values of Z
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
      #                      update b, theta1, and theta2                      #
      #                                                                        #
      ##########################################################################

      ##### update b

      ### compute un-normalized log-probabilities
      b_weights = logt_kernel_b(
        w = Z - eta1, b = b_space, eta2 = eta2, phi1 = phi1, phi2 = phi2, T = T
      )

      ### compute un-normalized probabilities rescaled by maximum
      b_weights = exp(b_weights - max(b_weights))

      ### sample new value of b
      b = sample(x = b_space, size = 1, prob = b_weights)



      ##### update theta2

      ### compute residuals of Z from the prior mean eta1
      z = Z - eta1

      ### compute 1 + eta2 * 1'_T %*% B^-1 %*% 1_T
      one_p_eta2_suminvB = (
        (1 + b)^2 * (eta2 * (T - 2) - 1) - 2 * (1 + b) * eta2 * (T - 1) +
          eta2 * T + 1
      ) / (1 - (1 + b)^2)

      ### compute (Z - eta1 * 1_T)' %*% invB %*% (Z - eta1 * 1_T)
      quad_form_z = (
        (1 + b)^2 * sum(z[-c(1, T)]^2) - 2 * (1 + b) * sum(z[-T] * z[-1]) +
          sum(z^2)
      ) / (1 - (1 + b)^2) - eta2 / one_p_eta2_suminvB *
        (
          (1 + b)^2 * sum(z[-c(1, T)])^2 - 2 * (1 + b) * sum(z[-c(1, T)]) *
            sum(z) + sum(z)^2
        ) / (2 + b)^2

      ### sample new value of theta2
      theta2 = 1 / rgamma(
        1, shape = phi1 + 0.5 * T, rate = phi2 + 0.5 * quad_form_z
      )



      ##### update theta1

      ### compute 1'_T %*% invB %*% Z
      rowsumsInvB_Z = (sum(Z) - (1 + b) * sum(Z[-c(1, T)])) /
        (2 + b)

      ### sample new value of theta1
      theta1 = rnorm(
        1, mean = (eta2 * rowsumsInvB_Z + eta1) / one_p_eta2_suminvB,
        sd = sqrt(eta2 * theta2 / one_p_eta2_suminvB)
      )



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
    print(paste(
      "GompPois_discrFlatNIG_mcmc: end time at ", Sys.time(), sep = ""
    ))
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