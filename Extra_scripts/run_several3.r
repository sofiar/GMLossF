nsims=c(30,50,100)
nreps=100
c=1
b=-0.24
b.init=b
source('like_lele_sim.r')
save.image('./like_lele_mcmc_v2_024FB.RData') 
