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
American_redStart <- readr::read_csv("./Real_Data_analysis/American_redStart2.csv")
#ggplot(American_redStart)+geom_point(aes(years,redstart))+theme_bw()

#str(American_redStart)
Nstar = American_redStart$redstart

source_dir = '/u/ruizsuar/GMLossF/Functions'

files = list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
for (ifun in files) source(ifun)

nsim = 1e+5
phi1 = 0.5
phi2 = 0.5
nu = 2
zeta1 = 0
zeta2 = 1
c = 2
verbose = 100

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

#mean(sqrt(post_sample_likelihood[2, ]))

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

# load("/u/ruizsuar/GMLossF/Rdata/realDataResults2.RData")

# plot(post_sample_likelihood[1, ], type = "l")
# plot(post_sample_likelihood[2, ], type = "l")
# plot(post_sample_likelihood[3, ], type = "l")

# plot(post_sample_complike$theta1, type = "l")
# points(post_sample_complike$theta1, type = "l", col = "red")
# plot(post_sample_complike$theta2, type = "l")
# plot(post_sample_complike$b, type = "l")



# difftime(time_end_complike, time_start_complike, "min")
# mean(sqrt(post_sample_complike$theta2))
# mean(post_sample_complike$theta1)
# mean(post_sample_complike$b)


# mean(sqrt(post_sample_likelihood[2,]))
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

