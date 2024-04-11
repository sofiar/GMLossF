################################################################################
############################ Extra functions ###################################
################################################################################

library(Rfast)
library(foreach)
library(future.apply)
library(doParallel)
library(future.batchtools)
library(MASS)
#library(Rcpp)
#library(RcppArmadillo)
#library(RcppEigen)
library(lamW)

#sourceCpp("optimfunctions.cpp")
#subfun2_bis_wrap(Ns = c(10, 20, 30), M = 5000, theta1 = 1.5, theta2 = 1, b = -0.8)
#subfun2(Ns = c(10, 20, 30), M = 5000, theta1 = 1.5, theta2 = 1, b = -0.8)


approx.CL2=function(Nobs,theta1,theta2,b,M,d=NULL)# d<length(Nobs)
{
curr=0
Times=length(Nobs)
mu=rep(theta1, 2)
if(is.null(d))
{d=Times}
if(d>Times)
{cat('d must be less than the length of the time series')}
for (t1 in 1:(Times-1)) {
limit=min(t1+d,Times)    
for (t2 in (t1+1):limit) {
B = (1+b)^(abs(outer(c(1,abs(t2-t1+1)),c(1,abs(t2-t1+1)) , "-"))) # get B matrix
cov_matrix=theta2*B
N.sample <- exp(mvrnorm(M, mu = mu, Sigma = cov_matrix))
logval = colSums(dpois(Nobs[c(t1,t2)],t(N.sample),log=TRUE))
curr = curr-log(M) + log( sum(exp(logval)))
}    
}
return(curr)
}


get.l=function(theta1t,theta2t,theta3t,Nobs,Nt) 
{
abs=to_abs(theta1t,theta2t,theta3t)
nreps=dim(Nobs)[2]
TT=dim(Nobs)[1]
out=0
for(j in 1:nreps)
{
N1obs=Nobs[,j]
N1t=Nt[,j]

z_hat=numeric(TT)
vars=numeric(TT)
z_hat[1]=theta2*N1obs[1]+theta1-lambertW0(theta2*exp(theta2*N1obs[1]+theta1))
vars[1]=theta2/(1+theta2*exp(-z_hat[1]))

for(n in 2:TT)
{
 z_hat[n]=abs[['sigma2']]*N1obs[n]+abs[['a']]+(1+abs[['b']])*z_hat[n-1]-lambertW0(theta2*exp(theta2*N1obs[n]+(abs[['a']]+(a+abs[['b']]*z_hat[n-1]))))
 vars[n]=abs[['sigma2']]/(1+abs[['sigma2']]*exp(-z_hat[n]))
}
out=out+sum(dnorm(N1t,mean=z_hat,sd=sqrt(vars),log=TRUE))

}
return(out)
}


lambertW_expArg = function(x)
{
if(x<=709)
{
out=lamW::lambertW0(exp(x))
}else{
  out=x - log(x) + log(x)/x + (log(x)*(log(x)-2))/(2*x^2) +
    log(x)*(2*log(x)^2-9*log(x)+6)/(6*x^3) +
    log(x)*(-12+36*log(x)-22*log(x)^2+3*log(x)^3)/(12*x^4)
}
return(out)
}

log_targetPoissonGauss = function(x,n,tauSq,mu){
  -exp(x) + x*n - 0.5/tauSq * (x-mu)^2
}

log_proposalGauss = function(x,tauSq,xi){
  - 0.5/tauSq * (x-xi)^2
}


Simulate_Gompertz=function(TT,a,b,sigma2)
{
theta1=-a/b
theta2=-sigma2/(b*(b+2))

# B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
# Z=theta1+sqrt(theta2)*t(chol(B))%*%rnorm(TT)
# Nt.obs = rpois(TT, exp(Z)) 

Nt=numeric(TT)
Nt[1]=exp(rnorm(1,theta1,sqrt(theta2)))
Es=rnorm(TT,0,sqrt(sigma2))
for(t in 2:TT)
{
Nt[t]=Nt[t-1]*exp(a+b*log(Nt[t-1])+Es[(t-1)])  
}

Nt.obs=numeric(TT)

for(t in 1:TT)
{Nt.obs[t]=rpois(1,lambda=Nt[t])}

return(Nt.obs)

}


##### Simulation Gompertz with log normal error

Simulate_Gompertz_ln = function(TT,a,b,sigma2,tau2)
{
  Y = numeric(TT)
  mu = -a/b
  var = sigma2/(1-((1+b)^2))+tau2
  Y[1] = rnorm(1,mu,sqrt(var))
  for(i in 2:length(Y))
  {
  mu = a+(1+b)*(mu+(var-tau2)/var*(Y[i-1]-mu))
  var = (1+b)^2*(var-tau2)/var*tau2+sigma2+tau2
  Y[i] = rnorm(1,mu,sqrt(var))
  } 
return(exp(Y))
}

## comment to try pull request