##############################################################################
################# simulation analysis 2: replicate Lele 2006 #################
##############################################################################
source('Extrafunctions.R')
source('mcmc_alg.R')

# b=-.24 # |1+b|<1 to be stationary
# b=-.8 

theta1 = 1.9244
theta2 = (0.4726)^2
theta2sqrt =sqrt(theta2)
a = -b*theta1
sigma2 = -b*(b+2)*theta2

set.seed(1606)
#nsims=c(30,50,100,500)
#nsims=50
#nreps=20

est_theta1=c()
est_theta2sqrt=c()
est_b=c()
sample_size=c()


start=proc.time()[3]

for(n in nsims)
{
 for(i in 1:nreps)
 {
 Nt.obs=Simulate_Gompertz(TT=n,a=a,b=b,sigma2=sigma2)
#  results=mcmc.quasiGibbs(iters=10000, burn=1000, n.chains=1, #theta1.init, theta2.init,
#                  b.init=c(-0.5,-0.1),
#                  Nt.obs,
#                  phi.1=0.5, phi.2=0.5, eta.1=0, eta.2=10,c=c)
                 
results=mcmc.quasiGibbs_v2(iters=10000, burn=1000, n.chains=1, #theta1.init, theta2.init,
                 b.init=b.init,theta1.init=c(0,0),theta2.init=c(1,1),
                 Nt.obs,
                 phi.1=0.5, phi.2=0.5, eta.1=0, eta.2=10,c=c)
est_theta1=c(est_theta1,mean(results$Keep.theta1[,1]))
est_theta2sqrt=c(est_theta2sqrt,mean(sqrt(results$Keep.theta2[,1])))
est_b=c(est_b,mean(results$Keep.b[,1]))
sample_size=c(sample_size,n)
}   
}

end=proc.time()[3]
time = end-start
time
#save.image('./like_lele_sims3.RData') # mcmc.quasiGibbs: priors phi.1=0.5, phi.2=0.5, eta.1=0, eta.2=10 
#save.image('./like_lele_sims4.RData') # mcmc.quasiGibbs_gprior
#save.image('./like_lele_sims5.RData') # mcmc.quasiGibbs: priors phi.1=0.5, phi.2=0.5, eta.1=0, eta.2=10 ,n=500
#save.image('./like_lele_sims6.RData') # mcmc.quasiGibbs: priors phi.1=0.5, phi.2=0.5, eta.1=0, eta.2=10 ,n=1000
