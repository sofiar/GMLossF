
library(Rfast)
library(foreach)
library(future.apply)
library(doParallel)
library(future.batchtools)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(lamW)


approx.CL2=function(Nobs,theta1,theta2,b,M)
{
#nsims=dim(Nobs)[2]
curr=0
Times=length(Nobs)
mu=rep(theta1, 2)
for (t1 in 1:(Times-1)) {
for (t2 in (t1+1):(Times)) {
B = (1+b)^(abs(outer(c(1,abs(t2-t1+1)),c(1,abs(t2-t1+1)) , "-"))) # get B matrix
cov_matrix=theta2*B
#cov_matrix=theta2*matrix(c(1,(1+b)^abs(t2-t1),(1+b)^(t2-t1),1),nrow=2)
N.sample <- exp(mvrnorm(M, mu = mu, Sigma = cov_matrix))
#logval =dpois(Nobs[t1],N.sample[,1],log=TRUE)+dpois(Nobs[t2],N.sample[,2],log=TRUE)
logval = colSums(dpois(Nobs[c(t1,t2)],t(N.sample),log=TRUE))
curr = curr-log(M) + log( sum(exp(logval)))
}    
}
return(curr)
}


theta1 = 1.9244
theta2 = (0.4726)^2
b=-.24 # |1+b|<1 to be stationary
TT=1000
B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt(theta2)*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )
eZstar = exp(Zstar)

b_s=-c(1:1000)/(1001)
ys=numeric(length(b_s))

for(i in 1:length(b_s))
{
ys[i]=approx.CL2(Nobs=eZstar,theta1,theta2,b_s[i],M=5000)
print(i)
}
plot(b_s,ys,type='l')


save.image('./compute_Pslikelihood.RData')
