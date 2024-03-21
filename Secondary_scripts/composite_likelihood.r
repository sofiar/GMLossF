
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
bs=c(-.24,-0.8) # |1+b|<1 to be stationary
TT=1000

b_s=-c(1:200)/(201)
y_s=matrix(NA,ncol=2,nrow=length(b_s))

for(j in 1:length(bs))
{
set.seed(1606)
b=bs[j]    
B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt(theta2)*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )
eZstar = exp(Zstar)

# parallelize code 
### get cluster
Sys.setenv(OPENBLAS_NUM_THREADS=1)
Sys.setenv(MKL_NUM_THREADS=1)
Sys.setenv(OMP_NUM_THREADS=1)

cluster = parallel::makeCluster(50)
### export environment to workers
parallel::clusterExport( cl = cluster, varlist = ls(all.names = TRUE) )

### set cluster seed
parallel::clusterSetRNGStream(cl = cluster, iseed = 1984)
parallel::clusterEvalQ(cl = cluster, expr = {library(MASS)})

y_s[,j]=parallel::parApply(cl=cluster,X=matrix(1:length(b_s)),MARGIN=1,FUN=function(xapp){
approx.CL2(Nobs=eZstar,theta1,theta2,b_s[xapp],M=5000)   
})    
### leave cluster
parallel::stopCluster(cluster)

# for(i in 1:length(b_s))
# {
# ys[i,j]=approx.CL2(Nobs=eZstar,theta1,theta2,b_s[i],M=5000)
# }
}

#save.image('./compute_Pslikelihood_parallel.RData')

# to plot 
load('./compute_Pslikelihood_parallel.RData')
plot(b_s[b_s< -0.05],y_s[,1][b_s< -0.05],type='l',col='red')
abline(v=bs[1],col='red')

plot(b_s[b_s< -0.05],y_s[,2][b_s< -0.05],type='l',col='blue')
abline(v=bs[2],col='blue')


# plot(b_s[b_s< -0.05],ys[b_s< -0.05],type='l')
# abline(v=b,col='red')

#save.image('./compute_Pslikelihood2.RData')
