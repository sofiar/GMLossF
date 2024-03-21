source('Extrafunctions.R')
source('mcmc_alg.R')


# prior parameters

 #theta1t0 = log(-a/b)
 #theta2t0 = log(sigma2/(-b(2+b)))
 #theta3t0 = log((1+b/2)/(1-b/2))

#theta1t0=log(1.9244)
#theta2t0=log(0.4726^2)
#theta3t0=log((0.12)/(1-0.12))

nrep=20
a= 2.1
b=-.8 # |1+b|<1 to be stationary
sigma2=0.09434


theta1t0 = log(-a/b)
theta2t0 = log(sigma2/(-b*(2+b)))
theta3t0 = log((-b/2)/(1+b/2))


mu.theta1=theta1t0
sd.theta1=0.7
mu.theta2=theta2t0
sd.theta2=0.5
#mu.theta3=0
#sd.theta3=1

#b=-.8
set.seed(100)
#thetas1t=rnorm(nrep,mean=mu.theta1,0.6)
thetas1t=rep(theta1t0, nrep)
#thetas2t=rep(theta2t0, nrep)
thetas2t=rnorm(nrep,mean=mu.theta2,sd.theta2)
thetas3t=rep(theta3t0, nrep)
#thetas3t=rnorm(nrep,mean=mu.theta3,sd.theta3)

covers1=numeric(nrep)
covers2=numeric(nrep)
times=numeric(nrep)
ar1s=numeric(nrep)
ar2s=numeric(nrep)

Mtheta1t=numeric(nrep)
Mtheta2t=numeric(nrep)
Stheta1t=numeric(nrep)
Stheta2t=numeric(nrep)

chains1=list()
chains2=list()


thetas1t[1]=theta1t0
thetas2t[1]=theta2t0
thetas3t[1]=theta3t0

nchains=2

for (ii in 1:nrep)
{
  
   theta1t=thetas1t[ii]
   theta2t=thetas2t[ii]
   theta3t=thetas3t[ii]
  #theta1t=log(1.9244)
  #theta2t=log(0.4726^2)
  #theta3t=log((0.12)/(1-0.12))
  source('simulate_Gompertz.R')
  
  # set inits
  #theta1.inits=rnorm(nchains,mean=mu.theta1,0.3)
  theta2.inits=rnorm(nchains,mean=mu.theta2,0.5)
  #theta3.inits=rnorm(nchains,mean=mu.theta3,0.5)
  # set real value
  theta1.inits=rep(theta1t,nchains)
  theta2.inits=rep(theta2t0,nchains)
  theta3.inits=rep(theta3t,nchains)

  aa=mcmc.al2(iters=3500, burn=1000, n.chains=nchains, theta1t.init=theta1.inits, 
  theta2t.init=theta2.inits,theta3t.init=theta3.inits,M=3500,
              # Observations
              Nt.obs=Nt.obs,
              # real parameters
              theta1t=theta1t,theta2t=theta2t,theta3t = theta3t,
              # prior parameters
              mu.theta1t=mu.theta1, sd.theta1t=sd.theta1, mu.theta2t=mu.theta2, sd.theta2t=sd.theta2,
              # proposal parameters
              p.theta1_1=0.05, p.theta2_1=0.03, #0.005; 0.01
	      p.theta1_2=0.01, p.theta2_2=0.01)
  
  #bb=mcmc.al1(iters=3000, burn=1000, n.chains=nchains, theta1t.init=theta1.inits, 
  #theta2t.init=theta2.inits,theta3t.init=theta3.inits,M=5000,
  #            # Observations
  #            Nt.obs=Nt.obs,
  #            # real parameters
  #            theta1t=theta1t,theta2t=theta2t,theta3t = theta3t,
  #            # prior parameters
  #            mu.theta1t=mu.theta1, sd.theta1t=sd.theta1, mu.theta2t=mu.theta2, sd.theta2t=sd.theta2,
  #            # proposal parameters
  #            p.theta1_1=0.02, p.theta2_1=0.05, 
	#      p.theta1_2=0.005, p.theta2_2=0.05)
  
  covers1[ii]=aa$cover1
  covers2[ii]=aa$cover2
  times[ii]=aa$time
  ar1s[ii]=aa$ar1
  ar2s[ii]=aa$ar2
  Mtheta1t[ii]=mean(aa$pseudo1)
  Mtheta2t[ii]=mean(aa$pseudo2)
  Stheta1t[ii]=sd(aa$pseudo1)
  Stheta2t[ii]=sd(aa$pseudo2)
  
 chains1[[ii]]=aa$pseudo1
 chains2[[ii]]=aa$pseudo2

 save.image("./LF1.RData")

  }  


save.image("./LF1.RData")

### compare times and resutls 
# start_time <- Sys.time()
# cu=approx.CL1(Nt.obs,Theta1,Theta2,M)
# end_time <- Sys.time()
# end_time - start_time
# 
# 
# start_time2 <- Sys.time()
# cu2=bis.approx.CL1(Nt.obs,Theta1,Theta2,M)
# end_time2 <- Sys.time()
# end_time2 - start_time2
# 
# start_time3 <- Sys.time()
# cu3=apprxCL1_cpp(Nt.obs, Theta1, Theta2,M)
# end_time3 <- Sys.time()
# end_time3 - start_time3



