source('Extrafunctions.R')
source('mcmc_alg.R')

set.seed(1984)
### Try MCMC 
#a = 0.461856
#b =-.8 # |1+b|<1 to be stationary
#sigma2=0.1200734

a= 2.1
b=-.8 # |1+b|<1 to be stationary
sigma2=0.1200734

theta1=-a/b
theta2=-sigma2/(b*(b+2))

source('simulate_Gompertz.R')
Nt.obs=Nt.obs[,1]

#plot(1:TT,Nt.obs[,1])
#lines(1:TT,Nt.obs[,1])
#acf(Nt.obs[,1])

source('mcmc_alg.R')
results=mcmc.quasiGibbs(iters=10000, burn=1000, n.chains=2, #theta1.init, theta2.init,
                  b.init=c(-0.2,-0.1),
                  # Observations
                  Nt.obs,
                  # prior parameters
                  phi.1=theta2^2/1+2, phi.2=(theta2^2+1)*theta2, eta.1=theta1, eta.2=1)

save.image("./data_results.RData")
#results$time
#results$Keep.b
#results$Keep.theta1
#results$Keep.theta2
#load("./data_results.RData")

#plot(results$Keep.b[,1],type='l',col='blue',xlab='b',ylab='')
#lines(results$Keep.b[,2],type='l',col='red')
#abline(h = b, col="black", lwd=3, lty=2)

#plot(results$Keep.theta1[,1],type='l',col='blue',xlab='theta1',ylab='')
#lines(results$Keep.theta1[,2],type='l',col='red')
#abline(h = theta1, col="black", lwd=3, lty=2)

#plot(results$Keep.theta2[,1],type='l',col='blue',xlab='theta2',ylab='')
#lines(results$Keep.theta2[,2],type='l',col='red')
#abline(h = theta2, col="black", lwd=3, lty=2)
