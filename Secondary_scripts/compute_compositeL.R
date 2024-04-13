####### compute composite likelihood using rcpp 

source('./Main_scripts/Extrafunctions.R')


theta1 = 1.9244
theta2 = (0.4726)^2
b=-.24 # |1+b|<1 to be stationary
TT=100
a = -b*theta1
sigma2 = -b*(b+2)*theta2
# Simulate process
Nsims=Simulate_Gompertz(TT=TT,a=a,b=b,sigma2=sigma2)



# try with rcpp
Rcpp::sourceCpp('./Functions/GompPois_composite_likelihood.cpp')

system.time(composite_cpp<-GompPois_composite_likelihood_cpp(Nstar=Nsims, theta1=theta1, 
            theta2=theta2, b=b, nsim=10000,d=3))



system.time(composite_r <- GompPois_composite_likelihood(Nstar=Nsims,theta1,
                         theta2,b,nsim=10000,d=3)) 
composite_r
exp()
  


