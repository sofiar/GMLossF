#######
rm(list=ls())
library(MASS)

setwd('/u/ruizsuar/GMLossF')
source('GompUnkGP_likelihood2.R')
source('approx_cl2.r')

load('./all_likelihoods_parallel4.RData') 
#################
nreps=500
b_s=-c(1:1000)/(1001)
y_stl=matrix(NA,ncol=nreps,nrow=length(b_s)) # true likelihood
y_scl=matrix(NA,ncol=nreps,nrow=length(b_s)) # composite likelihood
theta1 = 1.9244
theta2 = (0.4726)^2
#bs=c(-.8) # |1+b|<1 to be stationary
b=-.8
TT=500
B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix

print(Sys.time())
### get cluster
set.seed(1606)
Sys.setenv(OPENBLAS_NUM_THREADS=1)
Sys.setenv(MKL_NUM_THREADS=1)
Sys.setenv(OMP_NUM_THREADS=1)

cluster = parallel::makeCluster(50)
### export environment to workers
parallel::clusterExport( cl = cluster, varlist = ls(all.names = TRUE) )

### set cluster seed
parallel::clusterSetRNGStream(cl = cluster, iseed = 1984)
parallel::clusterEvalQ(
  cl = cluster, expr = {
    source('GompUnkGP_likelihood2.R') 
    source('approx_cl2.r')
}
)
for(j in 1:nreps)
{
#Generate data
#B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt(theta2)*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )
eZstar = exp(Zstar)
parallel::clusterExport( cl = cluster, varlist = ls() )

# compute true likelihood
y_stl[,j]=parallel::parApply(cl=cluster,X=matrix(1:length(b_s)),MARGIN=1,FUN=function(xapp){
GompUnkGP_likelihood2(theta1, theta2, b_s[xapp], exp(Zstar),10000, log = TRUE, details = FALSE)    
})

# compute composite likelihood of dimension 2
y_scl[,j]=parallel::parApply(cl=cluster,X=matrix(1:length(b_s)),MARGIN=1,FUN=function(xapp){
approx.CL2(Nobs=eZstar,theta1,theta2,b_s[xapp],M=5000)   
})    
print(paste(j,Sys.time(),sep='-'))

if (j%%10==0)
{
#save.image('./all_likelihoods_parallel2.RData')  # nreps=1500, TT=30   
#save.image('./all_likelihoods_parallel3.RData')  # nreps=1500, TT=100   
save.image('./all_likelihoods_parallel4.RData')  # nreps=500, TT=500   
}
}
### leave cluster
parallel::stopCluster(cluster)

save.image('./all_likelihoods_parallel4.RData') 
###################

# load('./all_likelihoods_parallel4.RData') 

# l_mean=apply(y_stl,1,mean,na.rm=TRUE)
# cl_mean=apply(y_scl,1,mean,na.rm=TRUE)

# png(file="compare_plots4.png")
# par(mfrow=c(1,2))
#  plot(b_s,l_mean,type='l',col='blue',main='Likelihood')
#  abline(v=b,col='blue')

#  plot(b_s,cl_mean,type='l',col='red',main='Composite Likelihood')
#  abline(v=b,col='red')

# dev.off()

