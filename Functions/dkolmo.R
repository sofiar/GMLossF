################################################################################
#
# Kolmogorov Distribution: Density Function
#
################################################################################


############### description
# This function approximates the Kolmogorov distribution density. 2 different
# representations are used: one for 0 < x < 1.207, the other for x >= 1.207.
# In both cases the series is truncated up to first 15 terms.


############### function
dkolmo = function(x, log = FALSE) {

  # x is the argument whose we want to compute the density;
  # log, if it is true then the log-density is returned.

  # no extra notes


  ###############--------------- function script

  ### get argument of Kolmogorov distribution as vector
  x0 = as.vector(x)

  ### compute number of element of x0
  n0 = length(x0)

  ### initialize vector of results
  res = rep(0, times = n0)

  ### find values less than 1.207
  i0 = x0 < 1.207

  ### compute for x < 1.207
  if (sum(i0) > 0) {

    ### get series's terms up to 15 and except the first
    ki0 = c(
      9, 25, 49, 81, 121, 169, 225, 289, 361, 441, 529, 625, 729, 841
    ) #(2 * c(2:15) - 1)^2

    ### compute 4 * x^2
    fourXsq = 4 * x0[i0]^2

    ### compute log-density
    res[i0] = -0.4673558279152179029126 - #0.5 * log(2 * pi) - log(4)
      4 * log(x0[i0]) + log(pi^2 - fourXsq) - pi^2 / (2 * fourXsq) +
      log1p(rowSums(
        outer(-fourXsq, ki0 * pi^2, FUN = "+") / (pi^2 - fourXsq) *
          exp(outer(-1 / (2 * fourXsq), (ki0 - 1) * pi^2, FUN = "*"))
      ))
  }

  ### compute for x >= 1.207
  if (sum(!i0) > 0) {

    ### get series's terms up to 15 and except the first
    kni0 = c(
      4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225
    ) #c(2:15)^2

    ### compute 2 * x^2
    twoXsq = 2 * x0[!i0]^2

    ### compute log-density
    res[!i0] = 2.079441541679835747658 + #log(8)
      log(x0[!i0]) - twoXsq + log1p(colSums(
      c(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1) * kni0 *
        exp( outer(kni0-1, -twoXsq, FUN = "*") )
    ))
  }

  ### if log is false then get the exponential of log-density
  if (log == FALSE) {
    res = exp(res)
  }

  ### get results in the same object form of x
  if (is.null(dim(x))) {
    return(res)
  } else {
    return(array(data = res, dim = dim(x)))
  }

}
