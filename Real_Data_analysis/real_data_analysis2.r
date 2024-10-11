#######################################################
################## Real data analysis #################
#######################################################
#library(readr)
#library(ggplot2)
#library(PVAClone)
#source('Extrafunctions.R')
#source('mcmc_alg.R')

set.seed(1984)

#load data
American_redStart = readr::read_csv("./Real_Data_analysis/American_redStart2.csv")


#str(American_redStart)
Nstar = American_redStart$redstart

# set.seed(20)
# set.seed(40)
# Nstar = GompPois_simulate(
#           30, theta1 = 1.9244, theta2 = 0.4726 ^ 2  , b = -0.24
#         )
        

source_dir = '/u/ruizsuar/GMLossF/Functions'

files = list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
for (ifun in files) source(ifun)

# Try hyperTuner approach. 

set.seed(5656)

outputEB = GompPois_hyperTuner(Nstar, maxit = 1e+3, abstol = sqrt(.Machine$double.eps),
  nsim = 1e+4, starter = NULL, verbose = 1)

nsim = 1e+4
phi1 = outputEB[3]
phi2 = outputEB[4]
zeta1 = outputEB[1]
zeta2 = outputEB[2]
c = 1
verbose = 100

postSample = GompPois_quasiGibbs(nsim = nsim, Nstar = Nstar, psi1 = phi1, 
psi2 = phi2, eta1 = zeta1, eta2 = zeta2, c = c, verbose = verbose)

mean(postSample$theta1)
mean(postSample$theta2)
mean(postSample$b)

# Lele's value
# theta1 = 1.9244
# theta2 = 0.2234
# b = -0.24

###############################################################################
nsim = 1e+4
phi1 = 0.5
phi2 = 0.5
nu = 2
zeta1 = 0
zeta2 = 1
c = 2
verbose = 100

# # new priors
# true_likelihood= GompPois_newPriors_quasiGibbs(
#   nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
#   zeta1 = zeta1, zeta2 = zeta2, c, verbose = verbose
# )

# Em algorithm
em_model= GompPois_mcmc_em(Nstar = Nstar, maxit = 1e+3, nsim = 1e+4, verbose =1) 


#
time_start_likelihood = Sys.time()
true_likelihood = GompPois_UScaledPriors_quasiGibbs(
  nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
  zeta1 = zeta1, zeta2 = zeta2, c, verbose = verbose
)
time_end_likelihood = Sys.time()

post_sample_likelihood = matrix(nrow = 3, ncol = nsim)
post_sample_likelihood[1, ] = true_likelihood$theta1
post_sample_likelihood[2, ] = true_likelihood$theta2
post_sample_likelihood[3, ] = true_likelihood$b

median(sqrt(post_sample_likelihood[2, ]))
mean(sqrt(post_sample_likelihood[2, ]))

median((post_sample_likelihood[2, ]))

mean((post_sample_likelihood[1, ]))


plot(post_sample_likelihood[2, ],type='p')
plot(density(post_sample_likelihood[2, ]))

#rowMeans(post_sample_likelihood)
#is.na(post_sample_likelihood[2, ])

### to avoid infinite values
post_sample_likelihood[2, post_sample_likelihood[2, ] == 0] = 0 + .Machine$double.eps
post_sample_likelihood[3, post_sample_likelihood[3, ] == -1] = -1 + .Machine$double.eps
post_sample_likelihood[3, post_sample_likelihood[3, ] == 0] = 0 - .Machine$double.eps


transf_post_sample = post_sample_likelihood
transf_post_sample[2 ,] = log(post_sample_likelihood[2, ])
transf_post_sample[3, ] = log(
  -post_sample_likelihood[3, ] / (1 + post_sample_likelihood[3, ])
)
delta = diff(t(transf_post_sample))
Omega = t(delta) %*% delta / (nsim - 1)
chol(Omega)

unlist(GompPois_MoM(Nstar)[c(1, 2, 3)])

time_start_complike = Sys.time()
post_sample_complike = GompPois_UScaledPriors_compositeLike(
  nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
  zeta1 = zeta1, zeta2 = zeta2, c, Omega = Omega, verbose = verbose
)
time_end_complike = Sys.time()

difftime(time_end_complike, time_start_complike, "min")
mean(sqrt(post_sample_complike$theta2))
mean(post_sample_complike$theta1)
mean(post_sample_complike$b)

rowMeans(post_sample_likelihood)
c(
  mean(post_sample_complike$theta1),
  mean(post_sample_complike$theta2),
  mean(post_sample_complike$b)
)
unlist(GompPois_MoM(Nstar)[c(1,2,3)])


# save.image("/u/ruizsuar/GMLossF/Rdata/realDataResults2.RData")

################### Check results ##########################

load("/u/ruizsuar/GMLossF/Rdata/realDataResults2.RData")

# plot(post_sample_likelihood[1, ], type = "l")
# plot(post_sample_likelihood[2, ], type = "l")
# plot(post_sample_likelihood[3, ], type = "l")

# plot(post_sample_complike$theta1, type = "l")
# points(post_sample_complike$theta1, type = "l", col = "red")
# plot(post_sample_complike$theta2, type = "l")
# plot(post_sample_complike$b, type = "l")



# difftime(time_end_complike, time_start_complike, "min")
# mean(sqrt(post_sample_complike$theta2))
# mean((post_sample_complike$theta2))
# mean(post_sample_complike$theta1)
# mean(post_sample_complike$b)


# difftime(time_end_likelihood, time_start_likelihood, "min")
# mean(sqrt(post_sample_likelihood[2,]))
# mean((post_sample_likelihood[2,]))
# mean(post_sample_likelihood[1,])
# mean(post_sample_likelihood[3,])

# par(mfrow = c(1,2))
# dev.off()

# sns::ess(post_sample_likelihood[1 ,])
# sns::ess(post_sample_complike$theta1)

# sns::ess(post_sample_likelihood[2 ,])
# sns::ess(post_sample_complike$theta2)

# sns::ess(post_sample_likelihood[3 ,])
# sns::ess(post_sample_complike$b)

# acceptance rate for MH in CL approach

d = diff(post_sample_complike$theta1)
sum(d!=0)/length(d)

