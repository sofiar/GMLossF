###############################################################################
#
#
# Gompertz Model with Poisson Sampling Error: Bayesian vs Frequentist
#
#
###############################################################################





############################### Initialization ################################

### set the seed
set.seed(1994)



##### set environment for simulation study

### number of iterations in simulation study
niter = 100

### number of scans in MCMC
nsim = 1e+4

### sample size
T = 30



##### get true parameter values
theta_Bay = c(1.90, 0.55, -0.20)
theta_freq = c(1.90, 0.20, -0.35)





################## Simulation Study #####################

##### parallel environment setup

### get cluster
cluster = parallel::makeCluster(24)

### get true number of iterations according to cluster size
niter = ceiling(niter / length(cluster)) * length(cluster)

### export variables
parallel::clusterExport(cl = cluster, varlist = ls(all.names = TRUE))

### identify cores
for (id_cluster in 1:length(cluster)) {
  parallel::clusterExport(cl = cluster[id_cluster], "id_cluster")
}

### create folder for reports file
suppressWarnings(dir.create("./reports"))

### open report file
parallel::clusterEvalQ(
  cl = cluster, sink(paste("./reports/report", id_cluster, ".txt", sep = ""))
)

### set seed
parallel::clusterSetRNGStream(cl = cluster)




########## do iterations in parallel for theta_Bay
res_Bay = parallel::clusterApply(
  cl = cluster, x = matrix(1:niter), fun = function(xapp) {

    ### print start of iteration
    print(paste(
      "iteration ", xapp, " start at time ", Sys.time()
    ))

    ### stopping variable for while loop
    toDo = TRUE

    ### run until get all successes
    while (toDo) {

      ### draw values
      Nstar = gse::GompPois_rng(
        T, theta1 = theta_Bay[1], theta2 = theta_Bay[2], b = theta_Bay[3]
      )

      ### method of moments
      res_mom = gse::GompPois_MoM(Nstar)
      res_mom = unlist(res_mom)[c(1, 2, 4)]

      ### maximum likelihood estimator
      res_mle = tryCatch(
        gse::GompPois_stEM_mle(nsim, Nstar = Nstar),
        error = function(err) return(NA)
      )
      if (is.na(res_mle)[1]) toDo = FALSE

      ### sparse normal inverse gamma prior
      if (toDo) {
        res_sp = tryCatch(
          gse::GompPois_flatNIG_mcmc(
            nsim, Nstar = Nstar, phi1 = 0.1, phi2 = 0.1, eta1 = 0, eta2 = 100
          ),
          error = function(err) return(NA)
        )
        if (is.na(res_sp[1])) toDo = FALSE
      }

      ### mixture of normal inverse gamma priors
      if (toDo) {
        res_mix = tryCatch(
          gse::GompPois_flatMixNIG_mcmc(
            nsim, Nstar = Nstar, kappa = 2, psi1 = 0.5, psi2 = 0.5,
            nu = 2, zeta1 = 0, zeta2 = 2
          ),
          error = function(err) return(NA)
        )
        if (is.na(res_mix[1])) toDo = FALSE
      }

      ### empirical Bayes for scale parameters of normal inverse gamma prior
      if (toDo) {
        hyper_eb = tryCatch(
          gse::GompPois_stEM_eb(nsim, Nstar = Nstar, phi1 = 2, eta1 = 0),
          error = function(err) return(NA)
        )
        if (is.na(hyper_eb[1])) {
          toDo = FALSE
        } else {
          res_eb = tryCatch(
            gse::GompPois_flatNIG_mcmc(
              nsim, Nstar = Nstar,
              phi1 = hyper_eb$phi1, phi2 = hyper_eb$phi2,
              eta1 = hyper_eb$eta1, eta2 = hyper_eb$eta2
            ),
            error = function(err) return(NA)
          )
          if (is.na(res_eb[1])) toDo = FALSE
        }
      }

      toDo = !toDo

    }



    ##### get values from the different methods

    ### results for theta1
    theta1 = c(
      mom = res_mom["theta1"], mle = rowMeans(res_mle)["theta1"],
      sp = rowMeans(res_sp)["theta1"], mix = rowMeans(res_mix)["theta1"],
      eb = rowMeans(res_eb)["theta1"]
    )
    names(theta1) = c("mom", "mle", "sp", "mix", "eb")

    ### results for theta2
    theta2 = c(
      mom = res_mom["theta2"], mle = rowMeans(res_mle)["theta2"],
      sp = rowMeans(res_sp)["theta2"], mix = rowMeans(res_mix)["theta2"],
      eb = rowMeans(res_eb)["theta2"]
    )
    names(theta2) = c("mom", "mle", "sp", "mix", "eb")

    ### results for b
    b = c(
      mom = res_mom["b"], mle = rowMeans(res_mle)["b"],
      sp = rowMeans(res_sp)["b"], mix = rowMeans(res_mix)["b"],
      eb = rowMeans(res_eb)["b"]
    )
    names(b) = c("mom", "mle", "sp", "mix", "eb")



    ### return results
    return(list(theta1 = theta1, theta2 = theta2, b = b))

  }
)




### flag for completion of theta_Bayes and start theta_freq
parallel::clusterEvalQ(
  cl = cluster, print(
    "Simulation for theta_Bayes completed, start simulation for theta_freq"
  )
)




########## do iterations in parallel for theta_freq
res_freq = parallel::clusterApply(
  cl = cluster, x = matrix(1:niter), fun = function(xapp) {

    ### print start of iteration
    print(paste(
      "iteration ", xapp, " start at time ", Sys.time()
    ))

    ### stopping variable for while loop
    toDo = TRUE

    ### run until get all successes
    while (toDo) {

      ### draw values
      Nstar = gse::GompPois_rng(
        T, theta1 = theta_freq[1], theta2 = theta_freq[2], b = theta_freq[3]
      )

      ### method of moments
      res_mom = gse::GompPois_MoM(Nstar)
      res_mom = unlist(res_mom)[c(1, 2, 4)]

      ### maximum likelihood estimator
      res_mle = tryCatch(
        gse::GompPois_stEM_mle(nsim, Nstar = Nstar),
        error = function(err) return(NA)
      )
      if (is.na(res_mle)[1]) toDo = FALSE

      ### sparse normal inverse gamma prior
      if (toDo) {
        res_sp = tryCatch(
          gse::GompPois_flatNIG_mcmc(
            nsim, Nstar = Nstar, phi1 = 0.1, phi2 = 0.1, eta1 = 0, eta2 = 100
          ),
          error = function(err) return(NA)
        )
        if (is.na(res_sp[1])) toDo = FALSE
      }

      ### mixture of normal inverse gamma priors
      if (toDo) {
        res_mix = tryCatch(
          gse::GompPois_flatMixNIG_mcmc(
            nsim, Nstar = Nstar, kappa = 2, psi1 = 0.5, psi2 = 0.5,
            nu = 2, zeta1 = 0, zeta2 = 2
          ),
          error = function(err) return(NA)
        )
        if (is.na(res_mix[1])) toDo = FALSE
      }

      ### empirical Bayes for scale parameters of normal inverse gamma prior
      if (toDo) {
        hyper_eb = tryCatch(
          gse::GompPois_stEM_eb(nsim, Nstar = Nstar, phi1 = 2, eta1 = 0),
          error = function(err) return(NA)
        )
        if (is.na(hyper_eb[1])) {
          toDo = FALSE
        } else {
          res_eb = tryCatch(
            gse::GompPois_flatNIG_mcmc(
              nsim, Nstar = Nstar,
              phi1 = hyper_eb$phi1, phi2 = hyper_eb$phi2,
              eta1 = hyper_eb$eta1, eta2 = hyper_eb$eta2
            ),
            error = function(err) return(NA)
          )
          if (is.na(res_eb[1])) toDo = FALSE
        }
      }

      toDo = !toDo

    }



    ##### get values from the different methods

    ### results for theta1
    theta1 = c(
      mom = res_mom["theta1"], mle = rowMeans(res_mle)["theta1"],
      sp = rowMeans(res_sp)["theta1"], mix = rowMeans(res_mix)["theta1"],
      eb = rowMeans(res_eb)["theta1"]
    )
    names(theta1) = c("mom", "mle", "sp", "mix", "eb")

    ### results for theta2
    theta2 = c(
      mom = res_mom["theta2"], mle = rowMeans(res_mle)["theta2"],
      sp = rowMeans(res_sp)["theta2"], mix = rowMeans(res_mix)["theta2"],
      eb = rowMeans(res_eb)["theta2"]
    )
    names(theta2) = c("mom", "mle", "sp", "mix", "eb")

    ### results for b
    b = c(
      mom = res_mom["b"], mle = rowMeans(res_mle)["b"],
      sp = rowMeans(res_sp)["b"], mix = rowMeans(res_mix)["b"],
      eb = rowMeans(res_eb)["b"]
    )
    names(b) = c("mom", "mle", "sp", "mix", "eb")



    ### return results
    return(list(theta1 = theta1, theta2 = theta2, b = b))

  }
)



##### parallel environment shut down

### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cl = cluster)



##### merge results from different workers

### initialize results
# Bayesian case
theta1_Bay_T30 = matrix(
  nrow = length(res_Bay[[1]]$theta1), ncol = niter,
  dimnames = list(c("mom", "mle", "sp", "mix", "eb"), NULL)
)
theta2_Bay_T30 = theta1_Bay_T30
b_Bay_T30 = theta1_Bay_T30
# frequentist case
theta1_freq_T30 = matrix(
  nrow = length(res_freq[[1]]$theta1), ncol = niter,
  dimnames = list(c("mom", "mle", "sp", "mix", "eb"), NULL)
)
theta2_freq_T30 = theta1_freq_T30
b_freq_T30 = theta1_freq_T30

### fill objects
for (iiter in 1:niter) {
  theta1_Bay_T30[, iiter] = res_Bay[[iiter]]$theta1
  theta2_Bay_T30[, iiter] = res_Bay[[iiter]]$theta2
  b_Bay_T30[, iiter] = res_Bay[[iiter]]$b

  theta1_freq_T30[, iiter] = res_freq[[iiter]]$theta1
  theta2_freq_T30[, iiter] = res_freq[[iiter]]$theta2
  b_freq_T30[, iiter] = res_freq[[iiter]]$b
}




################################### Closure ###################################
save(list = ls(), file = "SimStudy_Bay_vs_freq_T30.RData")