{
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  ### uncomment this for case T = 30
  # load("./SimStudy_correctModel_T30/SimStudy_correctModel_T30.RData")
  
  ### uncomment this for case T = 100
  # load("./SimStudy_correctModel_T100/SimStudy_correctModel_T100.RData")
  
  get_coverage = function(case) {
    if (case == "case1") {
      thetaTRUE = theta_case1
      CI95 = CI95_case1
    }
    if (case == "case2") {
      thetaTRUE = theta_case2
      CI95 = CI95_case2
    }
    ordered_thetaTRUE = rep(x = thetaTRUE, times = 3)
    ordered_thetaTRUE = ordered_thetaTRUE[c(1, 4, 7, 2, 5, 8, 3, 6, 9)]
    coverage = CI95[, "2.5%", ] <= ordered_thetaTRUE &
      CI95[, "97.5%", ] >= ordered_thetaTRUE
    return(matrix(
      rowMeans(coverage), nrow = 3, ncol = 3,
      dimnames = list(c("mle", "Gibbs", "stan"), c("theta1", "theta2", "b"))
    ))
  }
  
  get_times = function(case) {
    if (case == "case1") compTimes = compTimes_case1
    if (case == "case2") compTimes = compTimes_case2
    return(rowMeans(compTimes))
  }
  
  get_mse = function(case) {
    if (case == "case1") {
      thetaTRUE = theta_case1
      pointEst = pointEst_case1
    }
    if (case == "case2") {
      thetaTRUE = theta_case2
      pointEst = pointEst_case2
    }
    ordered_thetaTRUE = rep(x = thetaTRUE, times = 3)
    ordered_thetaTRUE = ordered_thetaTRUE[c(1, 4, 7, 2, 5, 8, 3, 6, 9)]
    return(rowMeans((pointEst - ordered_thetaTRUE)^2, dims = 2))
  }
  
  get_bias = function(case) {
    if (case == "case1") {
      thetaTRUE = theta_case1
      pointEst = pointEst_case1
    }
    if (case == "case2") {
      thetaTRUE = theta_case2
      pointEst = pointEst_case2
    }
    ordered_thetaTRUE = rep(x = thetaTRUE, times = 3)
    ordered_thetaTRUE = ordered_thetaTRUE[c(1, 4, 7, 2, 5, 8, 3, 6, 9)]
    return((rowMeans(pointEst, dims = 2) - ordered_thetaTRUE))
  }
  
  get_eff = function(case) {
    if (case == "case1") effSize = effSize_case1
    if (case == "case2") effSize = effSize_case2
    return(rowMeans(effSize, dims = 2))
  }
  
  get_effps = function(case) {
    if (case == "case1") {
      effSize = effSize_case1
      compTimes = compTimes_case1
    }
    if (case == "case2") {
      effSize = effSize_case2
      compTimes = compTimes_case2
    }
    tmp = matrix(
      0, nrow = 2, ncol = 3,
      dimnames = list(c("Gibbs", "stan"), c("theta1", "theta2", "b"))
    )
    for (i_alg in c("Gibbs", "stan")) {
      tmp[i_alg, "theta1"] = mean(effSize[i_alg, "theta1", ] /
                                    compTimes[i_alg, ])
      tmp[i_alg, "theta2"] = mean(effSize[i_alg, "theta2", ] /
                                    compTimes[i_alg, ])
      tmp[i_alg, "b"] = mean(effSize[i_alg, "b", ] / compTimes[i_alg, ])
    }
    return(tmp)
  }
}

theta_case1
theta_case2

get_mse("case1")
get_mse("case2")

get_bias("case1")
get_bias("case2")

get_coverage("case1")
get_coverage("case2")

get_times("case1")
get_times("case2")

get_eff("case1")
get_eff("case2")

get_effps("case1")
get_effps("case2")
