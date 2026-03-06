#' @title Gompertz Model with Poisson Sampling Error with
#'  Composite Likelihood: Empirical Bayes

#' @description
#' DO NOT USE THIS FUNCTION: CALIBRATION PRODUCES DEGENERATE PRIOR!!!!!!!!!!
#'  Calibration of hyper-parameters of the NIG prior for Gompertz model with
#'  Poisson sampling error distribution using composite likelihood approach.

#' @usage
#' GompPois_cl_eb = function(nsim, Nstar, starter = NULL, burn = 1e+3,
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
#'  b \sim U(-2, 0) \, , \\
#'  \theta_2 \sim Inv.Gamma(\phi_1, \phi_2) \, , \\
#'  \theta_1 \vert \theta_2 \sim N(\eta_1, \eta_2 \theta_2) \, .
#' }
#' The true likelihood is replaced with a composite likelihood, that is, TO BE FINISHED.
#' Thus, a stochastic expectation maximization algorithm is used to calibrate
#' \eqn{\phi_1, \phi_2, \eta_1, \eta_2}.
#' HOWEVER DO NOT USE THIS FUNCTION. See [INSERT REFERENCE] for more details.

#' @references
#' [INSERT REFERENCE]

#' @param nsim The number of simulations for the Markov chain.
#' @param Nstar The vector of the data.
#' @param starter List for starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @return Return a matrix containing a sequence of optimizers of the
#'         parameters.

#' @export
GompPois_cl_eb = function(
  nsim, Nstar, starter = NULL, burn = 1e+3, thin = 1, verbose = +Inf
) {

  # #--------------
  # sub_Gomp_logCL_1 = function(b, Z, theta1, theta2) {
  #   res = 0
  #   for (t in 2:T) {
  #     res = res + dnorm(
  #       Z[t], mean = -b * theta1 + (1 + b) * Z[t - 1],
  #       sd = sqrt(-b * (2 + b) * theta2), log = TRUE
  #     )
  #   }
  #   return(res)
  # }
  # #--------------

  #stop("Do not use this function! It produces a degenerate prior.")

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn < 0) burn = 0




  ########## initialization

  ### sample size
  T = length(Nstar)

  ### sample size for Monte Carlo
  N = 1e+4

  ### get

  ### get method of moments estimate of parameters
  res_MoM = GompPois_MoM(Nstar)
  # check if estimated values belong to the support of the parameters
  if (res_MoM$theta2 < 0.01) starter$theta2 = 0.01
  if (res_MoM$b > -0.01) starter$b = -0.01
  if (res_MoM$b < -1.99) starter$b = -1.99



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
    starter$phi2 = starter$theta2 * (0.5 * (T - 1) - 1)
  }

  # ### initialize parameter values
  # theta1 = starter$theta1
  # theta2 = starter$theta2
  # b = starter$b

  phi1 = starter$phi1
  phi2 = starter$phi2
  eta1 = starter$eta1
  eta2 = starter$eta2

  # a = -b * theta1
  # sigmaSq = -b * (2 + b) * theta2

  # ### initialize inverse of correlation matrix
  # invB = matrix(0, nrow = T, ncol = T)

  # ### initialize sub_mcmc results object
  # sub_mcmc = list(
  #   sample = matrix(
  #     c(theta1, theta2, b), nrow = 3,
  #     dimnames = list(c("theta1", "theta2", "b"), NULL)
  #   ),
  #   lastZ = Nstar # just for initialization, it is not used
  # )



  ### storing vectors
  # keep_phi1 = double(nsim)
  # keep_phi2 = double(nsim)
  # keep_eta1 = double(nsim)
  # keep_eta2 = double(nsim)



  ##### run optimization using univariate composite likelihood
  res = optim(
    par = c(phi1, phi2, eta1, eta2),
    fn = function(xfun) {
      theta2 = 1 / rgamma(N, shape = xfun[1], rate = xfun[2])
      theta1 = rnorm(N, mean = xfun[3], sd = sqrt(theta2 * xfun[4]))
      Z = matrix(0, nrow = T, ncol = N)
      Z[1, ] = sqrt_theta2 * sqrt(xfun[4] + 1) * rnorm(N)
      for (t in 1:T) {
        Z[t, ] = sqrt_theta2 * 0
      }
      Z = xfun[3] + sqrt((xfun[4] + 1) * xfun[2] / xfun[1]) *
        qt(u, df = 2 * xfun[1])
      Z_mat = matrix(Z, nrow = T, ncol = N, byrow = TRUE)
      return(
        -sum(
          log(
            rowSums(
              exp(Z_mat * Nstar - exp(Z_mat))
            )
          )
        )
      )
    },
    method = "L-BFGS-B", control = list(maxit = 1000),
    lower = c(0.01, 0.01, -Inf, 0.01), upper = rep(+Inf, 4)
  )

  ###
  names(res$par) = c("phi1", "phi2", "eta1", "eta2")

  ### return results
  return(res)

}
