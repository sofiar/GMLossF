########## do iterations in parallel for theta_case1
res_case1 = parallel::clusterApply(
  cl = cluster, x = matrix(1:niter), fun = function(xapp) {

    ### print start of iteration
    print(paste(
      "theta_case1: iteration ", xapp, " start at time ", Sys.time()
    ))

    ### stopping variable for while loop
    toDo = TRUE

    ### run until get all successes
    while (toDo) {

      ### draw values
      Nstar = gse::GompNegbinom_rng(
        T, theta1 = theta_case1[1], theta2 = theta_case1[2], b = theta_case1[3],
        prob = prob
      )

      ### maximum likelihood estimator
      time = Sys.time()
      res_mle = tryCatch(
        gse::GompPois_mcmcEM(
          Nstar = Nstar, maxit = maxit_em, max_nsim = max_nsim_em,
          nsim = nsim_em
        ),
        error = function(err) return(NA)
      )
      time_mle = difftime(Sys.time(), time, units = "secs")
      if (is.na(res_mle)[1]) toDo = FALSE

      ### Bayesian approach (Gibbs)
      if (toDo) {
        time = Sys.time()
        res_Gibbs = tryCatch(
          gse::GompPois_flatNIG_mcmc(
            nsim, Nstar = Nstar, phi1 = 0.1, phi2 = 0.1, eta1 = 0, eta2 = 100
          ),
          error = function(err) return(NA)
        )
        time_Gibbs = difftime(Sys.time(), time, units = "secs")
        if (is.na(res_Gibbs[1])) toDo = FALSE
      }

      ### Bayesian approach (stan)
      if (toDo) {
        object_stan = tryCatch(
          suppressWarnings(suppressMessages(rstan::sampling(
            object = stan_model, data = list(M = length(Nstar), Nstar = Nstar),
            chain = 1, cores = 1, iter = nsim + 1e+3, warmup = 1e+3,
            show_messages = FALSE, refresh = 0
          ))),
          error = function(err) return(NA)
        )
        if (class(object_stan) == "logical") {
          toDo = FALSE
        } else {
          time_stan = rstan::get_elapsed_time(object_stan)
          sample_stan = rstan::extract(object_stan, permuted = FALSE)
          res_stan = matrix(
            nrow = 3, ncol = nsim, dimnames = list(c("theta1", "theta2", "b"))
          )
          res_stan[1, ] = sample_stan[, 1, "theta1"]
          res_stan[2, ] = sample_stan[, 1, "theta2"]
          res_stan[3, ] = sample_stan[, 1, "b"]
        }
      }

      toDo = !toDo

    }



    ##### get results from the different methods

    ### computational times
    compTimes = double(3)
    names(compTimes) = c("mle", "Gibbs", "Stan")
    compTimes["mle"] = time_mle
    compTimes["Gibbs"] = time_Gibbs
    compTimes["Stan"] = sum(time_stan)

    ### point estimator
    pointEst = matrix(
      nrow = 3, ncol = 3,
      dimnames = list(c("mle", "Gibbs", "stan"), c("theta1", "theta2", "b"))
    )
    pointEst[1, ] = res_mle$par
    pointEst[2, ] = rowMeans(res_Gibbs)
    pointEst[3, ] = rowMeans(res_stan)

    ### effective sample size
    effSize = matrix(
      nrow = 2, ncol = 3,
      dimnames = list(c("Gibbs", "stan"), c("theta1", "theta2", "b"))
    )
    effSize[1, ] = sns::ess(t(res_Gibbs))
    effSize[2, ] = sns::ess(t(res_stan))



    ### CI at 95%
    CI95 = matrix(
      nrow = 9, ncol = 2,
      dimnames = list(c(
        "mle_theta1", "Gibbs_theta1", "stan_theta1",
        "mle_theta2", "Gibbs_theta2", "stan_theta2",
        "mle_b", "Gibbs_b", "stan_b"
      ))
    )
    CI95["mle_theta1", ] = pointEst[1, 1] + qnorm(c(0.025, 0.975)) * sqrt(
      res_mle$invFIM[1, 1]
    )
    CI95["mle_theta2", ] = pointEst[1, 2] + qnorm(c(0.025, 0.975)) * sqrt(
      res_mle$invFIM[2, 2]
    )
    CI95["mle_b", ] = pointEst[1, 3] + + qnorm(c(0.025, 0.975)) * sqrt(
      res_mle$invFIM[3, 3]
    )
    CI95["Gibbs_theta1", ] = quantile(
      res_Gibbs["theta1", ], probs = c(0.025, 0.975)
    )
    CI95["Gibbs_theta2", ] = quantile(
      res_Gibbs["theta2", ], probs = c(0.025, 0.975)
    )
    CI95["Gibbs_b", ] = quantile(
      res_Gibbs["b", ], probs = c(0.025, 0.975)
    )
    CI95["stan_theta1", ] = quantile(
      res_stan["theta1", ], probs = c(0.025, 0.975)
    )
    CI95["stan_theta2", ] = quantile(
      res_stan["theta2", ], probs = c(0.025, 0.975)
    )
    CI95["stan_b", ] = quantile(
      res_stan["b", ], probs = c(0.025, 0.975)
    )



    ### return results
    return(list(
      compTimes = compTimes,
      pointEst = pointEst,
      effSize = effSize,
      CI95 = CI95
    ))

  }
)




### flag for completion of theta_case1es and start theta_case2
parallel::clusterEvalQ(
  cl = cluster, print(
    "Simulation for theta_case1 is completed, start simulation for theta_case2"
  )
)




########## do iterations in parallel for theta_case2
res_case2 = parallel::clusterApply(
  cl = cluster, x = matrix(1:niter), fun = function(xapp) {

    ### print start of iteration
    print(paste(
      "theta_case2: iteration ", xapp, " start at time ", Sys.time()
    ))

    ### stopping variable for while loop
    toDo = TRUE

    ### run until get all successes
    while (toDo) {

      ### draw values
      Nstar = gse::GompPois_rng(
        T, theta1 = theta_case2[1], theta2 = theta_case2[2], b = theta_case2[3]
      )

      ### maximum likelihood estimator
      time = Sys.time()
      res_mle = tryCatch(
        gse::GompPois_mcmcEM(
          Nstar = Nstar, maxit = maxit_em, max_nsim = max_nsim_em,
          nsim = nsim_em
        ),
        error = function(err) return(NA)
      )
      time_mle = difftime(Sys.time(), time, units = "secs")
      if (is.na(res_mle)[1]) toDo = FALSE

      ### Bayesian approach (Gibbs)
      if (toDo) {
        time = Sys.time()
        res_Gibbs = tryCatch(
          gse::GompPois_flatNIG_mcmc(
            nsim, Nstar = Nstar, phi1 = 0.1, phi2 = 0.1, eta1 = 0, eta2 = 100
          ),
          error = function(err) return(NA)
        )
        time_Gibbs = difftime(Sys.time(), time, units = "secs")
        if (is.na(res_Gibbs[1])) toDo = FALSE
      }

      ### Bayesian approach (stan)
      if (toDo) {
        object_stan = tryCatch(
          rstan::sampling(
            object = stan_model, data = list(M = length(Nstar), Nstar = Nstar),
            chain = 1, cores = 1, iter = nsim + 1e+3, warmup = 1e+3,
            show_messages = FALSE, refresh = 0
          ),
          error = function(err) return(NA)
        )
        if (class(object_stan) == "logical") {
          toDo = FALSE
        } else {
          time_stan = rstan::get_elapsed_time(object_stan)
          sample_stan = rstan::extract(object_stan, permuted = FALSE)
          res_stan = matrix(
            nrow = 3, ncol = nsim, dimnames = list(c("theta1", "theta2", "b"))
          )
          res_stan[1, ] = sample_stan[, 1, "theta1"]
          res_stan[2, ] = sample_stan[, 1, "theta2"]
          res_stan[3, ] = sample_stan[, 1, "b"]
        }
      }

      toDo = !toDo

    }



    ##### get results from the different methods

    ### computational times
    compTimes = double(3)
    names(compTimes) = c("mle", "Gibbs", "Stan")
    compTimes["mle"] = time_mle
    compTimes["Gibbs"] = time_Gibbs
    compTimes["Stan"] = sum(time_stan)

    ### point estimator
    pointEst = matrix(
      nrow = 3, ncol = 3,
      dimnames = list(c("mle", "Gibbs", "stan"), c("theta1", "theta2", "b"))
    )
    pointEst[1, ] = res_mle$par
    pointEst[2, ] = rowMeans(res_Gibbs)
    pointEst[3, ] = rowMeans(res_stan)

    ### effective sample size
    effSize = matrix(
      nrow = 2, ncol = 3,
      dimnames = list(c("Gibbs", "stan"), c("theta1", "theta2", "b"))
    )
    effSize[1, ] = sns::ess(t(res_Gibbs))
    effSize[2, ] = sns::ess(t(res_stan))



    ### CI at 95%
    CI95 = matrix(
      nrow = 9, ncol = 2,
      dimnames = list(c(
        "mle_theta1", "Gibbs_theta1", "stan_theta1",
        "mle_theta2", "Gibbs_theta2", "stan_theta2",
        "mle_b", "Gibbs_b", "stan_b"
      ))
    )
    CI95["mle_theta1", ] = pointEst[1, 1] + qnorm(c(0.025, 0.975)) * sqrt(
      res_mle$invFIM[1, 1]
    )
    CI95["mle_theta2", ] = pointEst[1, 2] + qnorm(c(0.025, 0.975)) * sqrt(
      res_mle$invFIM[2, 2]
    )
    CI95["mle_b", ] = pointEst[1, 3] + + qnorm(c(0.025, 0.975)) * sqrt(
      res_mle$invFIM[3, 3]
    )
    CI95["Gibbs_theta1", ] = quantile(
      res_Gibbs["theta1", ], probs = c(0.025, 0.975)
    )
    CI95["Gibbs_theta2", ] = quantile(
      res_Gibbs["theta2", ], probs = c(0.025, 0.975)
    )
    CI95["Gibbs_b", ] = quantile(
      res_Gibbs["b", ], probs = c(0.025, 0.975)
    )
    CI95["stan_theta1", ] = quantile(
      res_stan["theta1", ], probs = c(0.025, 0.975)
    )
    CI95["stan_theta2", ] = quantile(
      res_stan["theta2", ], probs = c(0.025, 0.975)
    )
    CI95["stan_b", ] = quantile(
      res_stan["b", ], probs = c(0.025, 0.975)
    )



    ### return results
    return(list(
      compTimes = compTimes,
      pointEst = pointEst,
      effSize = effSize,
      CI95 = CI95
    ))

  }
)