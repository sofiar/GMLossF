#######################################################
############ Prior predictive check for b #############
#######################################################
library(readr)
library(ggplot2)
library(truncnorm)
library(sn)

source('./Functions/Extrafunctions.R')

# load real data
American_redStart = read_csv("./Real_Data_analysis/American_redStart2.csv")
ggplot(American_redStart)+geom_point(aes(years,redstart))+
geom_line(aes(years,redstart))+theme_bw()
n=length(American_redStart$redstart)
nsims=1000

# 1. Fix theta 1 and theta 2
theta1_s = rep(1.9244,nsims) 
theta2_s = rep((0.4726)^2,nsims)

# 2. Sampling theta 1 and theta 2
phi.1= 3 #1 #0.1
phi.2= 2 #0.1

phi.2/(phi.1-1)

eta.1= 3.5 #0
eta.2= 2.2635 #100

set.seed(100)
theta2_s = 1/rgamma(nsims,shape=phi.1,rate=phi.2)
theta1_s = rnorm(nsims,mean=eta.1,sd=sqrt(theta2_s*eta.2))
#theta1_s = rsn(nsims,xi=eta.1,omega=sqrt(eta.2*theta2_s),alpha=100)
#theta1_s= rtruncnorm(nsims, a=0, b=Inf, mean = eta.1, sd = sqrt(theta2_s*eta.2))
# sample b
bs = -runif(nsims)

sum(theta1_s<0)

# histogram of the priors
hist(bs)
hist(theta2_s)
hist(theta1_s)

Nsims=matrix(NA,ncol=n,nrow=nsims)
for (i in 1:nsims)
{
b=bs[i]
theta1=theta1_s[i]
theta2=theta2_s[i]
# compute parameter transformation
a = -b*theta1
sigma2 = -b*(b+2)*theta2
# Simulate process
Nsims[i,]=Simulate_Gompertz(TT=n,a=a,b=b,sigma2=sigma2)
}

# plot some of them 
q=sample(1:nsims,8)

par(mfrow=c(3,3))
plot(American_redStart$redstart,type='l',col='red',main='real data')
for(i in 1:8)
{
main=paste('b=',as.character(round(bs[q[i]],digits=2)),'\n',
'theta2=',as.character(round(theta2_s[q[i]],digits=2)),'\n',
'theta1=',as.character(round(theta1_s[q[i]],digits=2)))
plot(Nsims[q[i],],type='l',main=main)}
dev.off()

### Compute acf
par(mfrow=c(3,3))
acf(American_redStart$redstart,main='real data',col='red')
for(i in 1:8)
{acf(Nsims[q[i],],main=paste('b=',as.character(bs[q[i]])))}
dev.off()

# acf 1
tot_acf1=apply(Nsims, 1, acf,lag=1,plot=FALSE)
acf1=numeric(nsims)
for(i in 1:nsims)
{acf1[i]=tot_acf1[[i]]$acf[2]}
hist(acf1)
abline(v=acf(American_redStart$redstart,lag=1,plot=FALSE)$acf[2],lty=2,col='red',main=mean)

### Compute max, min, sd and mean 
tot_mean=apply(Nsims, 1, mean)
tot_max=apply(Nsims, 1, max)
tot_min=apply(Nsims, 1, min)
tot_sd=apply(Nsims, 1, sd)

par(mfrow=c(2,2))
boxplot(tot_mean)
abline(h=mean(American_redStart$redstart),lty=2,col='red',main=min)
boxplot(tot_min)
abline(h=min(American_redStart$redstart),lty=2,col='red',main=min)
boxplot(tot_sd)
abline(h=sd(American_redStart$redstart),lty=2,col='red',main=sd)
boxplot(tot_max)
abline(h=max(American_redStart$redstart),lty=2,col='red',main=max)
dev.off()

summary(tot_mean)
mean(American_redStart$redstart)
