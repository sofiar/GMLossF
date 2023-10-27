#########################################
######## Simulate Gompertz model ########
#########################################

library(tidyverse)
library(ggplot2)

# Nt+1=Nt exp(a+log(Nt)+et+1) #

# set parameters
# a= 0.46
# b=-.24 # |1+b|<1 to be stationary
# sigma2=0.19
#a= .46
#b=-.6 # |1+b|<1 to be stationary
# a=0.46
# b=-.8
# sigma2=0.19
# a= .3
# b=-.7 # |1+b|<1 to be stationary
# sigma2=0.1

# a=-theta1*b
# sigma2=theta2*(-b*(2+b))

#expt1=exp(theta1t)
#expt2=exp(theta2t)
#expt3=exp(theta3t)

#b=-2*(expt3)/(1+expt3)
#a=-b*expt1
#sigma2=-expt2*b*(2+b)

TT=500
nsims=1

Nt=matrix(NA,nrow=TT,ncol=nsims)
Nt[1,]=rep(10,nsims)
Es=rnorm(TT*nsims,0,sqrt(sigma2))

for (i in 1:nsims)
{
for(t in 2:TT)
  {
    Nt[t,i]=Nt[t-1,i]*exp(a+b*log(Nt[t-1,i])+Es[i*(t-1)])  
  }
  
}

#### Including sampling error 
Nt.obs=matrix(NA,nrow=TT,ncol=nsims)
for (i in 1:nsims)
{
for(t in 1:TT)
{Nt.obs[t,i]=rpois(1,lambda=Nt[t,i])}
} 

# to avoid -Inf values 
Nt.obs[Nt.obs==0]=0.01
### Transform data 
lN=log(Nt)
lN.obs=log(Nt.obs)

# theta1=-a/b
# theta2=sigma2/(-b*(2+b))

### Create data frame
sims=rep(1:nsims,each=TT)
times=rep(1:TT,nsims)
data=data.frame(c(Nt),c(Nt.obs),c(lN),c(lN.obs),times,sims)
names(data)=c('Nt','Nt.obs','lN','lN.obs','times','sims')
#data$lN.obs[!is.finite(data$lN.obs)] = NA

### Create 2-D pairs 

lN.obs[!is.finite(lN.obs)]=NA

# |i-j|=1

d1.obs=list()
d1=list()
for (i in 1:nsims)
{
  #d1[[i]]=log(cbind(Nt.obs[1:(TT-1),i],Nt.obs[2:(TT),i]))
  d1.obs[[i]]=cbind(lN.obs[1:(TT-1),i],lN.obs[2:(TT),i])
  d1[[i]]=cbind(lN[1:(TT-1),i],lN[2:(TT),i])
  
  # |i-j|=2
  #d2=cbind(Nt.obs[1:(TT-2)],Nt.obs[3:(TT)])
  
}



# ggplot(data %>% filter(sims==1))+geom_point(aes(x=times,y=Nt),col='red')+
#   geom_point(aes(x=times,y=Nt.obs))+theme_bw()+facet_wrap(~sims,scale='free')

Nt.obs[Nt.obs==0.01]=0
#### Check autocorrelations 

# acf(data %>% filter(sims==5) %>% select(Nt))
# acf(data %>% filter(sims==4) %>% select(Nt))
# acf(data %>% filter(sims==10) %>% select(Nt))


