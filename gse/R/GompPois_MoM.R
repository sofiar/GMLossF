#' @title Gompertz Model with Poisson Sampling Error: Method of Moments

#' @description
#' Point estimator using method of moments for Gompertz model with Poisson
#'  sampling error distribution

#' @usage
#' GompPois_MoM = function(Nstar)

#' @details
#' This function estimates the parameters of the following model:
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
#' The method of moments is used. In particular, the parameters are estimated
#' such that the mean, variance, and autocovariance at lag 1 are matched to the
#' empirical one.

#' @references
#' [INSERT REFERENCE]

#' @param Nstar The vector of the data.
#' @return Named list of estimates.

#' @export
GompPois_MoM = function(Nstar) {

  ### get empirical mean
  t1 = mean(Nstar)

  ### get empirical variance minus empirical mean
  t2prime = var(Nstar) - t1

  ### get empirical autocovariance at lag 1
  t3 = acf(Nstar, lag.max = 1, plot = FALSE, type = "covariance")$acf[2]

  ### get MoM estimates of the parameters
  theta2hat = log(1 + t2prime / t1^2)
  theta1hat = log(t1) - 0.5 * theta2hat
  theta3hat = log(1 + t3 * exp(-2 * theta1hat - theta2hat))

  ### get MoM estimates of the original parametrization
  b =  theta3hat / theta2hat - 1
  a = -theta1hat * b
  sigmaSq = -theta2hat * b * (2 + b)

  ### return the results
  return(list(
    theta1 = theta1hat,
    theta2 = theta2hat,
    a = a,
    b = b,
    sigmaSq = sigmaSq
  ))

}