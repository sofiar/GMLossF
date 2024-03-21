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


approx.CL1=function(Nobs,theta1,theta2,M)
{
nsims=dim(Nobs)[2]
Times=dim(Nobs)[1]
curll=0
#N.sample=exp(rnorm(M,mean=theta1,sd=sqrt(theta2)))
for (n in 1:nsims) {
 for (t in 1:Times) {
    N.sample=exp(Rnorm(M,m=theta1,s=sqrt(theta2)))
    dp=dpois(Nobs[t,n],lambda = N.sample)
    # to avoid -Inf values 
    dp[dp==0]=0.0001
    curll=curll+log(mean(dp))
  #print(curll)

    }
}
return(curll)
}

approx.CL2=function(Nobs,theta1,theta2,theta3,M)
{
nsims=dim(Nobs)[2]
Times=dim(Nobs)[1]
curll=0
b=-2*(theta3)/(1+theta3)
cov_matrix=theta2*matrix(c(1,1+b,1+b,1),nrow=2)
mu=rep(theta1, 2)
for (n in 1:nsims) {
 for (t in 1:(Times-1)) {
    N.sample <- exp(mvrnorm(M, mu = mu, Sigma = cov_matrix))
    positive_rows <- apply(N.sample, 1, function(row) all(row > 0))
    N.sample <- N.sample[positive_rows, , drop = FALSE]
    dp1=dpois(Nobs[t,n], lambda = N.sample[,1])
    dp2=dpois(Nobs[t+1,n], lambda = N.sample[,2])
    #dp=cbind(dp1,dp2)
    dp=dp1*dp2
    dp[dp == 0] <- 0.0001
    
    curll <- curll + log(mean(dp))
   
  #print(curll)

    }
}
return(curll)
}

subfun=function(Ns,M,theta1,theta2)
{
  aa=0
  Times=length(Ns)
  for (t in 1:Times) {
  N.sample=exp(Rnorm(M,m=theta1,s=sqrt(theta2)))
  dp=dpois(Ns[t],lambda = N.sample)
  # to avoid -Inf values 
  dp[dp==0]=0.0001
  aa=aa+log(mean(dp))
  }
return(aa)
}

subfun2=function(Ns,M,theta1,theta2,b)
{
  aa=0
  Times=length(Ns)
  cov_matrix=theta2*matrix(c(1,1+b,1+b,1),nrow=2)
  mu=rep(theta1, 2)
  #foreach(t = 1:(Times-1)) %do% {
  for (t in 1:(Times-1)) {
    N.sample = exp(mvrnorm(M, mu = mu, Sigma = cov_matrix))
    positive_rows <- apply(N.sample, 1, function(row) all(row > 0))
    N.sample <- N.sample[positive_rows, , drop = FALSE]
    dp1=dpois(Ns[t], lambda = N.sample[,1])
    dp2=dpois(Ns[t+1], lambda = N.sample[,2])
    #dp=cbind(dp1,dp2)
    dp=dp1*dp2
    dp[dp == 0] <- 0.0001
    aa=aa+log(mean(dp))
  }
return(aa)
}



bis.approx.CL1=function(Nobs,theta1,theta2,M)
{
  nsims=dim(Nobs)[2]
  Times=dim(Nobs)[1]
  plan(multicore,workers=10)
  #plan(multisession,workers = 100)
  #plan(multisession,workers = 8)
  a=future_apply(Nobs,2,subfun,M,theta1,theta2,future.seed = NULL)
  curll=sum(a)
  
  return(curll)
}



bis.approx.CL2=function(Nobs,theta1,theta2,theta3,M)
{
  nsims=dim(Nobs)[2]
  Times=dim(Nobs)[1]
  plan(multicore,workers=128)
  #plan(multisession,workers = 100)
  #plan(multisession,workers = 10)
  b=-2*(theta3)/(1+theta3)
  a=future_apply(Nobs,2,subfun2,M,theta1,theta2,b,future.seed = NULL)
  curll=sum(a)
  
  return(curll)
}


bis.approx.cpp.CL2=function(Nobs,theta1,theta2,b,M)
{
  nsims=dim(Nobs)[2]
  Times=dim(Nobs)[1]
  #plan(multicore,workers=10)
  #plan(multisession,workers = 100)
  #plan(multisession,workers = 10)
  #b=-2*(theta3)/(1+theta3)
  a=apply(Nobs,2,subfun2_bis_wrap,M,theta1,theta2,b)
  curll=sum(a)
  
  return(curll)
}
to_abs=function(theta1t,theta2t,theta3t)
{
  expt1=exp(theta1t)
  expt2=exp(theta2t)
  expt3=exp(theta3t)

  b=-2*(expt3)/(1+expt3)
  a=-b*expt1
  sigma2=-expt2*b*(2+b)
  out=c(a=a,b=b,sigma2=sigma2)
  return(out)
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

Nt=numeric(TT)
theta1=-a/b
theta2=-sigma2/(b*(b+2))

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

