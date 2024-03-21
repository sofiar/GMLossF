#######
rm(list=ls())
source('GompUnkGP_likelihood2.R')
setwd('/u/ruizsuar/GMLossF')
#load("all_likelihoodb05.RData")
#load('./all_likelihoods2.RData') 

#theta1_s=c(1:1000)/(1001)*2 + 1
b_s=-c(1:200)/(201)
#y_values=c()
y_s=matrix(NA,ncol=2,nrow=length(b_s))
theta1 = 1.9244
theta2 = (0.4726)^2
bs=c(-.8,-0.24) # |1+b|<1 to be stationary
TT=1000

print(Sys.time())

for(j in 1:length(bs))
{
#Generate data
b=bs[j]
set.seed(1606)
B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt(theta2)*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )
eZstar = exp(Zstar)

### get cluster
Sys.setenv(OPENBLAS_NUM_THREADS=1)
Sys.setenv(MKL_NUM_THREADS=1)
Sys.setenv(OMP_NUM_THREADS=1)

cluster = parallel::makeCluster(50)
### export environment to workers
parallel::clusterExport( cl = cluster, varlist = ls(all.names = TRUE) )

### set cluster seed
parallel::clusterSetRNGStream(cl = cluster, iseed = 1984)
parallel::clusterEvalQ(
  cl = cluster, expr = {source('GompUnkGP_likelihood2.R')}
)

y_s[,j]=parallel::parApply(cl=cluster,X=matrix(1:length(b_s)),MARGIN=1,FUN=function(xapp){
GompUnkGP_likelihood2(theta1, theta2, b_s[xapp], exp(Zstar),10000, log = TRUE, details = FALSE)    
})    
### leave cluster
parallel::stopCluster(cluster)
}

save.image('./all_likelihoods_parallel.RData') 

# to check and plot results 
# load('./all_likelihoods_parallel.RData') 
#plot(b_s[b_s < -0.05],y_s[,1][b_s < -0.05],type='l',col='red')
#abline(v=bs[1],col='red', lty = 2)

# plot(b_s[b_s < -0.05],y_s[,2][b_s < -0.05],type='l',col='blue')
# abline(v=bs[2],col='blue', lty = 2)

# code without parallelization 
#y_s=numeric(length(b_s))    
# for(i in 1:length(b_s))
# {
# y_s[i]=GompUnkGP_likelihood2(theta1, theta2, b_s[i], exp(Zstar),10000, log = TRUE, details = FALSE)
# }


#y_values=c(y_values,y_s) 
#b_values=c(b_values,b_s) 
#print(Sys.time())

#}   
#print(Sys.time())
#save.image('./all_likelihoodb05.RData') 

# plots for likelihood theta1
# plot(theta1_s,y_s,type='l')
# abline(v=theta1,col='red', lty = 2)

# plots for likelihood bs 

# load('./all_likelihoods_parallel.RData') 
#  plot(b_s[b_s < -0.05],y_s[b_s < -0.05],type='l',col='red')
#  abline(v=b,col='red', lty = 2)
# lines(b_s[b_s < -0.05],y_values[1001:2000][b_s < -0.05],type='l',col='blue')
#  abline(v=bss[2],col='blue', lty = 2)



# b_s2 = b_s^2
# b_s3 = b_s^3
# lm_fit=lm(y_s[b_s < -0.05] ~ b_s[b_s < -0.05] + b_s2[b_s < -0.05]+b_s3[b_s < -0.05])
# cs=coef(lm_fit)
# plot(x = b_s,y = y_s,type='l')
# plot(x = b_s[b_s < -0.05], y = exp(y_s[b_s < -0.05]),type='l')
# plot(x = b_s[b_s < -0.05], y = y_s[b_s < -0.05],type='l')
# lines(x = b_s[b_s < -0.05],cs[1]+cs[2]*b_s[b_s < -0.05]+cs[3]* b_s2[b_s < -0.05]+cs[4]* b_s3[b_s < -0.05],col='blue')


# plot(x = b_s[b_s < -0.05],exp(cs[1]+cs[2]*b_s[b_s < -0.05]+cs[3]* b_s2[b_s < -0.05]+cs[4]* b_s3[b_s < -0.05])/1.583974e-320,col='blue',type='l')

# abline(v=b,col='red', lty = 2)
# abline(h = 1.7e-62, col = "green3", lty = 3)

# integrate(f=function(x){exp(x^3*cs[4]+exp(x^2)*cs[3]+exp(x)*cs[2]+cs[1])},lower=-1,upper=0, subdivisions = 1e+4)
# x = b_s
# sum(exp(x^3*cs[4]+exp(x^2)*cs[3]+exp(x)*cs[2]+cs[1]))

