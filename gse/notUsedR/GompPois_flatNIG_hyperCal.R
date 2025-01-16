#' @title Gompertz Model with Poisson Sampling Error under Flat and
#' Normal-Inverse Gamma Prior: Hyper-Parameters Calibration

#' @description
#' Hyper-parameters calibration of normal-inverse gamma prior using method of
#'  moments for Gompertz model with Poisson sampling error distribution.

#' @usage
#' GompPois_flatNIG_hyperCal = function(Nstar)

#' @details
#' This function calibrates the hyper-parameters of prior distribution.
#'  The model is:
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
#' It is set:
#' \deqn{
#'  \phi_1 = 2 \, , \, \phi_2 = \widehat{\theta}_2 \, , \, \eta_1 =
#'  \widehat{\theta}_1 \, , \, \eta_2 =
#'  \frac{T (1 - (1 + \widehat{b})^2)}{T + (T-2)(1+\widehat{b})^2 -2(T-1)
#'   (1+\widehat{b})}
#' }
#' The method of moments is used for \eqn{\widehat{\theta}_1,
#'  \widehat{\theta}_2, \widehat{b}}. In this way, the prior means are centered
#'  to an unbiased estimator, the prior variance of \eqn{\theta_2} is infinite,
#'  and the variance of \eqn{\theta_1} arises from Zellner-g prior with \eqn{
#'  g = T}.

#' @references
#' [INSERT REFERENCE]

#' @param Nstar The vector of the data.
#' @return Named list of estimates.

#' @export
GompPois_flatNIG_hyperCal = function(Nstar) {

  ### initialize res
  res = double(4)
  # give names
  names(res) = c("phi1", "phi2", "eta1", "eta2")

  ### get the sample size
  T = length(Nstar)

  ### get method of moments estimator
  res_mom = GompPois_MoM(Nstar)

  ### fill results
  res["phi1"] = 2
  res["phi2"] = res_mom$theta2
  res["eta1"] = res_mom$theta1
  res["eta2"] = (
    T + (T - 2) * (1 + res_mom$b)^2 - 2 * (T - 1) * (1 + res_mom$b)
  ) / (
    1 - (1 + res_mom$b)^2
  )

  ### return results
  return(res)

}