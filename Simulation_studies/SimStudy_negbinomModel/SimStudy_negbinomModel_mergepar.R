##### merge results from different workers for theta_case1

### initialize merged objects
compTimes_case1 = matrix(
  nrow = 3, ncol = niter,
  dimnames = list(c("mle", "Gibbs", "stan"))
)

pointEst_case1 = array(
  dim = c(3, 3, niter), dimnames = list(
    c("mle", "Gibbs", "stan"),
    c("theta1", "theta2", "b")
  )
)

effSize_case1 = array(
  dim = c(2, 3, niter), dimnames = list(
    c("Gibbs", "stan"),
    c("theta1", "theta2", "b")
  )
)

CI95_case1 = array(
  dim = c(9, 2, niter), dimnames = list(
    c(
      "mle_theta1", "Gibbs_theta1", "stan_theta1",
      "mle_theta2", "Gibbs_theta2", "stan_theta2",
      "mle_b", "Gibbs_b", "stan_b"
    ),
    c("2.5%", "97.5%")
  )
)

### fill merged objects
for (i_niter in 1:niter) {
  compTimes_case1[, i_niter] = res_case1[[i_niter]]$compTimes
  pointEst_case1[, , i_niter] = res_case1[[i_niter]]$pointEst
  effSize_case1[, , i_niter] = res_case1[[i_niter]]$effSize
  CI95_case1[, , i_niter] = res_case1[[i_niter]]$CI95
}



##### merge results from different workers for theta_case2

### initialize merged objects
compTimes_case2 = matrix(
  nrow = 3, ncol = niter,
  dimnames = list(c("mle", "Gibbs", "stan"))
)

pointEst_case2 = array(
  dim = c(3, 3, niter), dimnames = list(
    c("mle", "Gibbs", "stan"),
    c("theta1", "theta2", "b")
  )
)

effSize_case2 = array(
  dim = c(2, 3, niter), dimnames = list(
    c("Gibbs", "stan"),
    c("theta1", "theta2", "b")
  )
)

CI95_case2 = array(
  dim = c(9, 2, niter), dimnames = list(
    c(
      "mle_theta1", "Gibbs_theta1", "stan_theta1",
      "mle_theta2", "Gibbs_theta2", "stan_theta2",
      "mle_b", "Gibbs_b", "stan_b"
    ),
    c("2.5%", "97.5%")
  )
)

### fill merged objects
for (i_niter in 1:niter) {
  compTimes_case2[, i_niter] = res_case2[[i_niter]]$compTimes
  pointEst_case2[, , i_niter] = res_case2[[i_niter]]$pointEst
  effSize_case2[, , i_niter] = res_case2[[i_niter]]$effSize
  CI95_case2[, , i_niter] = res_case2[[i_niter]]$CI95
}