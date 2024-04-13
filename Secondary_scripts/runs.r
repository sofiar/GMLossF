source('Extrafunctions.R')
source('mcmc_alg.R')

set.seed(1984)
### Try MCMC 
#a = 0.461856
#b =-.8 # |1+b|<1 to be stationary
#sigma2=0.1200734

a= 0.461856
b=-.24 # |1+b|<1 to be stationary
sigma2=0.09434338


theta1=-a/b
theta2=-sigma2/(b*(b+2))

Nt.obs=Simulate_Gompertz(TT=50,a=a,b=b,sigma2=sigma2)
#Nt.obs=Nt.obs[1:30]
plot(Nt.obs,type='l')
lines(Nt.obs,type='p')


results=mcmc.quasiGibbs_gprior(iters=10000, burn=1000, n.chains=1, #theta1.init, theta2.init,
                  b.init=c(-0.2,-0.1),
                  Nt.obs,
                  phi.1=0.01, phi.2=0.01, eta.1=0)

#results=mcmc.quasiGibbs(iters=10000, burn=1000, n.chains=1, #theta1.init, theta2.init,
#                  b.init=c(-0.2,-0.1),
#                  # Observations
#                  Nt.obs,
#                  # prior parameters
#                  #phi.1=theta2^2/1+2, phi.2=(theta2^2+1)*theta2, eta.1=theta1, eta.2=1
#                  phi.1=0.1, phi.2=0.1, eta.1=0, eta.2=100)

save.image("./data_results.RData")
#load("./data_results.RData")
#results$time
#results$Keep.b
#results$Keep.theta1
#results$Keep.theta2
#load("./data_results.RData")

plot(results$Keep.b[,1],type='l',col='blue',xlab='b',ylab='')
#lines(results$Keep.b[,2],type='l',col='red')
mean(results$Keep.b[,1])
#quantile(results$Keep.b[,1], probs = c(0.025, 0.975))

#abline(h = b, col="black", lwd=3, lty=2)
#acf(results$Keep.theta1[,1])
#acf(results$Keep.theta1[,2])
#sns::ess(results$Keep.theta1[,1])
#sns::ess(results$Keep.theta1[,2])


plot(results$Keep.theta1[,1],type='l',col='blue',xlab='theta1',ylab='')
mean(results$Keep.theta1[,1])
#lines(results$Keep.theta1[,2],type='l',col='red')
#abline(h = theta1, col="black", lwd=3, lty=2)
#acf(results$Keep.theta2[,1])
#acf(results$Keep.theta2[,2])
#sns::ess(results$Keep.theta2[,1])
#sns::ess(results$Keep.theta2[,2])

plot(results$Keep.theta2[,1],type='l',col='blue',xlab='theta2',ylab='')
mean(results$Keep.theta2[,1])
#lines(results$Keep.theta2[,2],type='l',col='red')
#abline(h = theta2, col="black", lwd=3, lty=2)
#acf(results$Keep.b[,1])
#acf(results$Keep.b[,2])
#sns::ess(results$Keep.b[,1])
#sns::ess(results$Keep.b[,2])


