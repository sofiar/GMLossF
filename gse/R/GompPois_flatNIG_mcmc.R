#' @title Gompertz Model with Poisson Sampling Error under Flat and
#' Normal-Inverse Gamma Prior: Markov Chain Monte Carlo Algorithm

#' @description
#' Posterior simulation using elliptical slice sampling within Gibbs for
#'  Gompertz model with Poisson sampling error distribution under flat
#'  and normal-inverse gamma priors.

#' @usage
#' GompPois_flatNIG_mcmc = function(nsim, Nstar, phi1, phi2, eta1, eta2,
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
GompPois_flatNIG_mcmc = function(
  nsim, Nstar, phi1, phi2, eta1, eta2, starter = NULL,
  burn = 0, thin = 1, verbose = +Inf, lastZ = FALSE
) {

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
  
  

  ### discrete prior for b (to get starting point for optimizer)
  b_space = c(1:199) / 100 - 2
  #b_weights = double(length = length(b_space))

  ### sample size
  T = length(Nstar)

  ### storing vectors
  keep_theta1 = double(nsim)
  keep_theta2 = double(nsim)
  keep_b = double(nsim)




  ########## draw the chain
  result <- mcmc_loop_cpp(
    burn, nsim, thin, verbose,
    Nstar, T, b_space,
    eta1, eta2, phi1, phi2,
    theta1, theta2, b, Z,
    update_Z_cpp, update_pars_cpp
  )

  nZ          <- length(Z)
  keep_b      <- result[1:nsim]
  keep_theta1 <- result[(nsim + 1):(2 * nsim)]
  keep_theta2 <- result[(2 * nsim + 1):(3 * nsim)]
  # b           <- result[3 * nsim + 1]
  # theta1      <- result[3 * nsim + 2]
  # theta2      <- result[3 * nsim + 3]
  Z           <- result[(3 * nsim + 4):(3 * nsim + 3 + nZ)]



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