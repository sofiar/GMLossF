###############################################################
##################### Simulation analysis #####################
###############################################################

library(sns)
source('Extrafunctions.R')
source('mcmc_alg.R')


set.seed(99)

#set paramters 
#b_s=c(-0.2,-0.5,-0.8)
#a_s=c(2.1,5.24,8.4)
#sigma2_s=c(0.01,0.06,0.12)

b_s=c(-0.2,-0.5,-0.8)
a_s=c(1.2,2.1,3.2)
sigma2_s=c(0.01,0.06,0.12)


keep.time=c()
keep.mse.b=c()
keep.mse.t1=c()
keep.mse.t2=c()
keep.ess.b=c()
keep.ess.t1=c()
keep.ess.t2=c()
keep.cov.b=c()
keep.cov.t1=c()
keep.cov.t2=c()
keep.b=c()
keep.t1=c()
keep.t2=c()

iters=1000
burn=500

for (bi in b_s)
{
b=bi    
    for (ai in a_s)
    {
    a=ai
        for (si in sigma2_s)
        {
        sigma2=si
        theta1=-a/b
        theta2=-sigma2/(b*(b+2))
            for (n in 1:1)
            {
            # simulate time series 
            Nt.obs=Simulate_Gompertz(TT=30,a=a,b=b,sigma2=sigma2)
            # run mcmc 
            results=mcmc.quasiGibbs(iters=iters, burn=burn, n.chains=1,
            b.init=-runif(1), Nt.obs, phi.1=0.1, phi.2=0.1, eta.1=0, eta.2=100)
            ## save values 
            # time
            keep.time=c(keep.time,results$time)            
            # mse
            keep.mse.b=c(keep.mse.b,mean((results$Keep.b[burn:iters]-bi)^2))
            keep.mse.t1=c(keep.mse.t1,mean((results$Keep.theta1[burn:iters]-theta1)^2))
            keep.mse.t2=c(keep.mse.t2,mean((results$Keep.theta2[burn:iters]-theta2)^2))
            # ess 
            keep.ess.b=c(keep.ess.b,sns::ess(results$Keep.b))
            keep.ess.t1=c(keep.ess.t1,sns::ess(results$Keep.theta1))
            keep.ess.t2=c(keep.ess.t2,sns::ess(results$Keep.theta2))
            # cov 
            # b
            low=quantile(results$Keep.b[burn:iters],probs = c(0.025))
            high=quantile(results$Keep.b[burn:iters],probs = c(0.975))
            keep.cov.b=c(keep.cov.b, as.numeric((b > low & b < high)))
            # theta1
            low=quantile(results$Keep.theta1[burn:iters],probs = c(0.025))
            high=quantile(results$Keep.theta1[burn:iters],probs = c(0.975))
            keep.cov.t1=c(keep.cov.t1, as.numeric((theta1 > low & theta1 < high)))
            # theta2
            low=quantile(results$Keep.theta2[burn:iters],probs = c(0.025))
            high=quantile(results$Keep.theta2[burn:iters],probs = c(0.975))
            keep.cov.t2=c(keep.cov.t2, as.numeric((theta2 > low & theta2 < high)))

            keep.b=c(keep.b,b)
            keep.t1=c(keep.t1,theta1)
            keep.t2=c(keep.t2,theta2)

            }
        }
    }
}

save.image("./simulation_results.RData")

# 
#