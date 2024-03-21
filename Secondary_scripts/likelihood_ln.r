rm(list=ls())
setwd('/u/ruizsuar/GMLossF')
source('GompUnkGL_likelihood.r')
source('Extrafunctions.R')

nreps=1000
b_s=-c(1:200)/(201)
b_real=c(-.8,-0.24)
TTs=c(30,50,100,500,1000)
save.likelihoods= array(NA,dim = c(length(b_s),2,5))

for (i in 1:length(b_real))
{
b=b_real[i]
for (k in 1:length(TTs))
{
TT=TTs[k]
#b=-.8
#TT=50
y_stl=matrix(NA,ncol=nreps,nrow=length(b_s)) # likelihood

theta1 = 1.9244
a = -b * theta1
theta2 = (0.4726)^2
sigma2 = -b*(2+b)*theta2
#Pseudo=Simulate_Gompertz(TT=10000,a=a,b=b,sigma2=sigma2)    
tau2 = 0.2315

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
  cl = cluster, expr = {source('GompUnkGL_likelihood.r')}
)
for(j in 1:nreps)
{
#Generate data
#B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
eZstar = Simulate_Gompertz_ln(TT,a,b,sigma2,tau2)
parallel::clusterExport( cl = cluster, varlist = c('eZstar'))


# compute true likelihood
y_stl[,j]=parallel::parApply(cl=cluster,X=matrix(1:length(b_s)),MARGIN=1,FUN=function(xapp){
GompUnkGL_likelihood(a, b_s[xapp], sigma2, tau2, eZstar, log = TRUE)    
})
}
### leave cluster
parallel::stopCluster(cluster)

l_mean=apply(y_stl,1,mean,na.rm=TRUE)
save.likelihoods[,i,k]=l_mean
save.image('./likelihood_allLN.RData')
print(paste('b=',b,'and  TT=',TT,sep=''))
}
}

save.image('./likelihood_allLN.RData')
###################################################################################
# plots from runs
load('./likelihood_allLN.RData') 

dim(save.likelihoods)
# bs, b_true_ , TT

for(i in 1:2)
{
png(file=paste("lnError_likelihood",as.character(b_real[i]),".png",sep=''))
par(mfrow=c(2,3),cex.axis=1,cex.lab=1,cex=1)
for(k in 1:5)
{
plot(b_s,save.likelihoods[,i,k],type='l',
col='blue',main=paste('Likelihood T=',TTs[k]),cex=10)
abline(v=b_real[i],col='blue')  
}
dev.off()
}



# l_mean=apply(y_stl,1,mean,na.rm=TRUE)

#  png(file="lnError_likelihood30.png")
#  plot(b_s,l_mean,type='l',col='blue',main='Likelihood')
#  abline(v=b,col='blue')
#  dev.off()







