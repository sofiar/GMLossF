##################################################################################
####################### get simulation data results ##############################
##################################################################################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyverse)


#### Lele's analysis 

load('./like_lele_mcmc_v2_024FB.RData')

data.results=as.data.frame(cbind(est_theta1,est_theta2sqrt,est_b,sample_size))

res_b=data.results%>%group_by(sample_size)%>%
summarize(q1=quantile(est_b, probs = c(0.025)),mean=mean(est_b),q3=quantile(est_b, probs = c(0.975)),
MAD=mean(abs(est_b-b)))

res_theta1=data.results%>%group_by(sample_size)%>%
summarize(q1=quantile(est_theta1, probs = c(0.025)),mean=mean(est_theta1),q3=quantile(est_theta1, probs = c(0.975)),
MAD=mean(abs(est_theta1-theta1)))

res_theta2=data.results%>%group_by(sample_size)%>%
summarize(q1=quantile(est_theta2sqrt, probs = c(0.025)),mean=mean(est_theta2sqrt),q3=quantile(est_theta2sqrt, probs = c(0.975)),
MAD=mean(abs(est_theta2sqrt-theta2sqrt)))

plot(results$Keep.b,type='l',col='blue',ylim=c(-0.5,0))
abline(h=b,col='blue')
hist(results$Keep.b)

plot(sqrt(results$Keep.theta2),type='l')
abline(h=sqrt(theta2),col='red')

plot((results$Keep.theta1),type='l')
abline(h=theta1,col='red')



ggplot(data.results)+geom_histogram(aes(est_b))+
facet_grid(~sample_size)+theme_bw()+geom_vline(xintercept=b,linetype=2)+
theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))

ggplot(data.results)+geom_histogram(aes(est_theta1))+
facet_grid(~sample_size)+theme_bw()+geom_vline(xintercept=theta1,linetype=2)+
theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))

ggplot(data.results)+geom_histogram(aes(est_theta2sqrt))+
facet_grid(~sample_size)+theme_bw()+geom_vline(xintercept=theta2sqrt,linetype=2)+
theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))



#### Analysis 2 
load("./simulation_results.RData")


res_data=as.data.frame(cbind(keep.time,keep.mse.b,keep.mse.t1,keep.mse.t2,
keep.ess.b,keep.ess.t1,keep.ess.t2,keep.cov.b,keep.cov.t1,
keep.cov.t2,keep.b,keep.t1,keep.t2))

names(res_data)=c('time','mse.b','mse.t1','mse.t2','ess.b',
'ess.t1','ess.t2','cov.b','cov.t1','cov.t2','b','theta1','theta2')
View(res_data)
max(res_data$time)
mean(res_data$time)


names(res_data)
ggplot(res_data)
