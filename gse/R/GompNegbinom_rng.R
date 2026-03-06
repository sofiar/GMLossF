#' @title Gompertz Model with Negative Binomial Sampling Error: Random Number
#'  Generator

#' @description
#' Generation of independent draws from Gompertz model with negative binomial
#'  sampling error distribution.

#' @details This function sample from the following model: \deqn{ N^\star_t |
#' N_t \sim \mathrm{negbinom}(p = prob, mu = N_t) \, , \\
#'  N_t = \exp(Z_t) \, , \\
#'  Z_t = a + (1 + b) Z_{t - 1} + \varepsilon_t \, , \\
#'  \varepsilon_t \overset{i.i.d.}{\sim} N(0, \sigma^2) \, . }
#'
#' A different parametrization is used, say: \deqn{ \theta_1 = -\frac{a}{b} \, ,
#'  \, \theta_2 = -\frac{\sigma^2}{b(2 + b)} \, , } which implies \deqn{ Z_1,
#'  Z_2, ..., Z_T \sim N_T ( \theta_1 1_T, \theta_2 \bar{\Sigma})  \, , \\
#'  \bar{\Sigma}_{ij} = (1 + b ) ^ {\vert i - j \vert} \, . }
#' 
#' As \code{prob} goes to one, the model converges to the Gompertz Poisson
#' sampling error distribution

#' @usage
#' GompNegbinom_rng = function(T, theta1, theta2, b, prob)

#' @references
#' [INSERT REFERENCE]

#' @param T The length of the chain, say the sample size.
#' @param theta1 The mean of the stationary distribution.
#' @param theta2 The variance of the stationary distribution.
#' @param b The autoregressive coefficient.
#' @param prob Probability parameter for negative binomial.
#' @return A vector containing the noised chain.

#' @export
GompNegbinom_rng = function(T, theta1, theta2, b, prob) {

  ### initialize result object
  res = double(T)

  ### get the intercept
  a = -b * theta1

  ### get the conditional variance
  sigmaSq = -b * (2 + b) * theta2
  
  ### get size for negative binomial
  size = prob / (1 - prob) * exp(theta1 + 0.5 * theta2)

  ### sample the first draw from the stationary distribution
  res[1] = theta1 + sqrt(theta2) * rnorm(1)

  ### sample the remaining draws from the conditional distribution
  for (t in 2:T) {
    res[t] = a + (1 + b) * res[t - 1] + sqrt(sigmaSq) * rnorm(1)
  }

  ### sample observed variables from Poisson
  res = rnbinom(T, size = size, mu = exp(res))

  return(res)

}