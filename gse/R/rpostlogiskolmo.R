#' @title Posterior Logistic Kolmogorov Distribution: Random Number Generator

#' @description
#' Exact random number generator for posterior logistic Kolmogorov distribution.

#' @usage rpostlogiskolmo = function(n, x)

#' @details
#' It is said \eqn{W \sim \mathrm{LK}(\cdot)}, say a logistic Kolmogorov
#'  distribution, if \eqn{W = 4 V^2} and \eqn{V} follows a Kolmogorov
#'  distribution.
#'
#' The posterior logistic Kolmogorov distribution is the random variable
#'  \eqn{W \vert X} where \eqn{W \sim \mathrm{LK}(\cdot)} and
#'  \eqn{X \vert W \sim N(0, W)}.
#'
#' It is used the acceptance-rejection scheme described in Onorati and Liseo
#'  (2022)

#' @references
#' Onorati, P., & Liseo, B. (2022). An extension of the unified skew-normal
#'  family of distributions and application to Bayesian binary regression.
#'  arXiv preprint arXiv:2209.03474.

#' @param n The number of simulations.
#' @param x Vector of conditioning values.
#' @return Vector of draws from posterior logistic Kolmogorov distribution.

#' @export
rpostlogiskolmo = function(n, x) {

  ##### initialization

  ### vector of results
  res = double(length = n)

  ### status of each simulation
  isToBeSampled = rep(TRUE, n)

  ### number of each simulation to be sampled
  sum_isToBeSampled = n

  ### vector of rate parameters of prior proposal
  ratePriorProp = 0.5 * pi^2 # NOT change from 0.5 * pi^2, otherwise see theory

  ### vector of logarithm of pre_iM
  logPre_iM = double(length = n)



  ##### get coherence between dimension

  ### work with a vector of absolute value of x
  x0 = as.vector(abs(x))

  ### make compatible length of x0 with n
  # length of x0 i.e. lenX0
  lenX0 = length(x0)
  # add alements to x0 such that lenX0 = n
  if (n >= lenX0) {
    # integer division
    nINTlenX0 = n %/% lenX0
    # modulus
    nMODlenX0 = n %% lenX0
    # add x0 nINTlenX0 times
    x0 = rep(x0, nINTlenX0)
    # if modulus greater than 0, then add the first nMODlenX0
    if (nMODlenX0 > 0) x0 = c(x0, x0[1:nMODlenX0])
    # if n > lenX0 then take only the first n elements of x0
  } else {
    x0 = x0[1:n]
  }




  ########## setup parameters of acceptance-rejection algorithm

  ### shape parameters of prior proposal
  # pre calculus
  shapePriorProp = 1 + exp(-x0)
  # compute posterior means
  shapePriorProp = x0 * shapePriorProp + shapePriorProp * (
    log(shapePriorProp) + exp(
      x0 + log(log(shapePriorProp))
    )
  )
  # compute shape parameters
  shapePriorProp = 0.5 * (1 + (pi^2 + x0^2) / shapePriorProp)



  ##### compute logPre_iM

  ### get vEqual
  vEqual = 1.9834 # NOT change from 1.9834, otherwise see theory

  ### get vStar i.e. arg max secondMaj assuming that shapePriorProp >= pi -1
  vStar = 1 + shapePriorProp + sqrt((1 + shapePriorProp)^2 - pi^2)

  ### compute first majorant i.e. firstMaj
  firstMaj = 2.1811794892486414 + # log(sqrt(2 * pi) * pi^2 * vEqual^(-1.5))
    lgamma(shapePriorProp) - shapePriorProp * log(ratePriorProp) +
    shapePriorProp*log(vEqual)

  ### compute second majorant i.e. secondMaj
  secondMaj = lgamma(shapePriorProp) - shapePriorProp * log(ratePriorProp) +
    (1 + shapePriorProp) * log(vStar) + ratePriorProp / vStar - 0.5 * vStar

  ### get the maximum
  logPre_iM = -pmax(firstMaj, secondMaj)




  ########## sampling until desired number is accepted
  while (sum_isToBeSampled > 0) {

    ##### select setup parameters of draws to be sampled

    ### x0
    selX0 = x0[isToBeSampled]

    ### shape prior proposal
    selShapePriorProp = shapePriorProp[isToBeSampled]

    ### logPre_iM
    selLogPre_iM = logPre_iM[isToBeSampled]



    ##### proposal distribution

    ### draws from posterior proposal i.e. propDraws
    propDraws = 1 / rgamma(
      sum_isToBeSampled,
      shape = selShapePriorProp + 0.5,
      rate = ratePriorProp + 0.5 * selX0^2
    )

    ### compute prior proposal log-density
    logPropDIGamma = dgamma(
      1 / propDraws,
      shape = selShapePriorProp,
      rate = ratePriorProp,
      log = TRUE
    ) - 2 * log(propDraws)

    ### compute prior target log-density
    logTarDLogisKolmo = dkolmo(
      0.5 * sqrt(propDraws), log = TRUE
    ) - 1.386294361119890572454 - # it's log(4)
      0.5 * log(propDraws)



    ##### accept-reject

    ### compute logarithmic accept-reject ratio
    logARratio = selLogPre_iM + logTarDLogisKolmo - logPropDIGamma

    ### sample log-uniform i.e. sample from negative exponential with rate 1
    logU = -rexp(sum_isToBeSampled, rate = 1)

    ### identify success
    succ = logU <= logARratio

    ### store success
    res[isToBeSampled][succ] = propDraws[succ]

    ### update status for each simulation
    isToBeSampled[isToBeSampled] = !succ

    ### update residual number of simulations
    sum_isToBeSampled = sum_isToBeSampled - sum(succ)

  }




  ########## organize and take res

  ##### if number of simulations less than number of x 
  if (n < lenX0) {

    ### case x is a vector
    if (is.null(dim(x))) {
      # fill missing simulations with NA
      res = c(res, rep(NA, lenX0 - n))
      # report a warning
      warning(paste(
        "length of x is greater than the number of simulations, ",
        "there are some NA's in the results vector", sep = ""
      ))
      # return res
      return(res)

      ### case x is an array
    } else {
      # fill missing simulations with NA
      res = array(c(res, rep(NA, lenX0 - n)), dim = dim(x))
      # report a warning
      warning(paste(
        "length of x is greater than the number of simulations, ",
        "there are some NA's in the results array", sep = ""
      ))
      # check if res is a 2D array
      if (length(dim(x)) == 2) {
        res = matrix(res, nrow = dim(x)[1], ncol = dim(x)[2])
      }
      # rerutn res
      return(res)
    }

    ##### if number of simulations equal or greater than number of x 
  } else {

    ### case number of x is not a multiple of number of simulations
    if (nMODlenX0 > 0) {
      # case x is a vector
      if (is.null(dim(x))) {
        # fill missing simulations with NA
        res = matrix(
          c(res, rep(NA, lenX0 - nMODlenX0)), nrow = lenX0, ncol = nINTlenX0 + 1
        )
        # report a warning
        warning(paste(
          "length of x is not a multiple of the number of simulations, ",
          "there are some NA's in the results matrix", sep = ""
        ))
        # return res
        return(res)

        # case x is an array
      } else {
        # fill missing simulation with NA
        res = array(
          c(res, rep(NA, lenX0 - nMODlenX0)), dim = c(dim(x), nINTlenX0 + 1)
        )
        # report a warning
        warning(paste(
          "length of x is not a multiple of the number of simulations, ",
          "there are some NA's in the results array", sep = ""
        ))
        # check if res is a 2D array
        if (length(dim(x)) == 2) {
          res = matrix(res, nrow = lenX0, ncol = nINTlenX0 + 1)
        }
        # return res
        return(res)
      }

      ### case number of x is a multiple of number of simulations
    } else {
      # case x is a scalar
      if (length(x) == 1) {
        return(res)
      } else {
        # case x is a vector
        if (is.null(dim(x))) {
          # fill matrix
          res = matrix(res, nrow = lenX0, ncol = nINTlenX0)
          # return res
          return(res)
          # case x is an array
        } else {
          # fill array
          res = array(res, dim = c(dim(x), nINTlenX0))
          # check if res is a 2D array
          if (length(dim(x)) == 2) {
            res = matrix(res, nrow = lenX0, ncol = nINTlenX0)
          }
          # return res
          return(res)
        }
      }
    }

  }

}