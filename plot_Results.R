### plot results
library(ggplot2)
library(tidyverse)
library(patchwork)
source('Extrafunctions.R')
load("LF1.RData")

nrep
ii
nrep=ii
burn=2000
n=1
### create df
iters=length(chains1[[n]][,1])
true_values=rep(thetas1t,each=nrep*nchains)
true_values2=rep(thetas2t,each=nrep*nchains)
all_values=unlist(chains1,recursive = FALSE)
all_values2=unlist(chains2,recursive = FALSE)
n_chain=rep(rep(1:2,each=iters),nrep)
n_rep=rep(rep(1:nrep),each=iters*2)
n_iter=rep(rep(1:iters,2),nrep)
all_theta1=data.frame(true_values,all_values,n_chain,n_rep,n_iter)
all_theta2=data.frame(true_values2,all_values2,n_chain,n_rep,n_iter)
all_theta1$n_chain=as.factor(all_theta1$n_chain)
all_theta2$n_chain=as.factor(all_theta2$n_chain)

## Plots - theta 1 - theta 2
## MCMC
p_mcmc=list()
p_mcmc2=list()
for (i in 1:nrep)
{
  p_mcmc[[i]]=ggplot(all_theta1 %>% filter(n_rep==i))+geom_line(aes(x=n_iter,y=all_values,col=n_chain),alpha=.8, show.legend = F)+
    theme_bw()+ ylim(thetas1t[i]-0.1,thetas1t[i]+0.1)+ 
  geom_hline(yintercept=thetas1t[i], linetype="dashed", color = "black")+xlab('iteration')+ylab('')

  p_mcmc2[[i]]=ggplot(all_theta2 %>% filter(n_rep==i))+geom_line(aes(x=n_iter,y=all_values2,col=n_chain),alpha=.8, show.legend = F)+
    theme_bw()+ ylim(thetas2t[i]-0.25,thetas2t[i]+0.25)+ 
  geom_hline(yintercept=thetas2t[i], linetype="dashed", color = "black")+xlab('iteration')+ylab('')

}

# theta1
p_mcmc[[1]] + p_mcmc[[2]]+p_mcmc[[3]]+p_mcmc[[4]]+p_mcmc[[5]]+p_mcmc[[6]]+p_mcmc[[7]]+
  p_mcmc[[8]]+p_mcmc[[9]]+p_mcmc[[10]]+p_mcmc[[11]]+p_mcmc[[12]]

# theta2
  p_mcmc2[[1]] + p_mcmc2[[2]]+p_mcmc2[[3]]+p_mcmc2[[4]]+p_mcmc2[[5]]+p_mcmc2[[6]]+p_mcmc2[[7]]+
    p_mcmc2[[8]]+p_mcmc2[[9]]+p_mcmc2[[10]]+p_mcmc2[[11]]+p_mcmc2[[12]]+p_mcmc2[[13]]+p_mcmc2[[14]]+
  p_mcmc2[[15]]+p_mcmc2[[16]]+p_mcmc2[[17]]+p_mcmc2[[18]]+p_mcmc2[[19]]+p_mcmc2[[20]]



## density 
p_density=list()
for (i in 1:nrep)
{
  p_density[[i]]=ggplot()+geom_density(data=all_theta1 %>% filter(n_rep==i,n_iter>burn),aes(x=all_values),alpha=.8, show.legend = F)+
  theme_bw() + geom_vline(xintercept=thetas1t[i], linetype="dashed", color = "black")+xlab('theta1')+ylab('')
}

p_density[[1]] + p_density[[2]]+p_density[[3]]+p_density[[4]]+p_density[[5]]+p_density[[6]]+
  p_density[[7]]+  p_density[[8]]+p_density[[9]]+p_density[[10]]+p_density[[11]]+p_density[[12]]


## histogram 
p_hist=list()
for (i in 1:nrep)
{
  p_hist[[i]]=ggplot()+geom_histogram(data=all_theta1 %>% filter(n_rep==i,n_iter>burn),aes(x=all_values),alpha=.8, show.legend = F)+
    theme_bw() + geom_vline(xintercept=thetas1t[i], linetype="dashed", color = "black")+xlab('theta1')+ylab('')
}

p_hist[[1]] +p_hist[[2]] +p_hist[[3]] +p_hist[[4]] +p_hist[[5]] +p_hist[[6]] +
  p_hist[[7]] +p_hist[[8]] +p_hist[[9]] +p_hist[[10]] +p_hist[[11]] +p_hist[[12]] 


# coverage
mean(covers1[1:ii])
mean(covers2[1:ii])


a_s=c()
b_s=c()
sigma2_s=c()

for (i in 1:ii)
{
  a_s[i]=to_abs(thetas1t[i],thetas2t[i],thetas3t[i])[1]
  b_s[i]=to_abs(thetas1t[i],thetas2t[i],thetas3t[i])[2]
  sigma2_s[i]=to_abs(thetas1t[i],thetas2t[i],thetas3t[i])[3]
}



covers1[1:ii]
thetas1t[1:ii]
a_s



plot(covers1[1:ii], thetas1t[1:ii])
plot(covers1[1:ii],a_s)
plot(covers1[1:ii],b_s)
plot(covers1[1:ii],sigma2_s)



# MSE
mse1=numeric(nrep)
mse2=numeric(nrep)

for(n in 1:ii)
{
  mse1[n]=mean((chains1[[n]][1000:3000]-thetas1t[n])^2)#/var(chains1[[n]][1000:3500])
  mse2[n]=mean((chains2[[n]][1000:3000]-thetas2t[n])^2)#/var(chains1[[n]][1000:3500])
  
}

boxplot(mse1)
boxplot(mse1,mse2)


## Aceptance ratio
ac1_a=numeric(nrep)
ac1_b=numeric(nrep)
ac1_c=numeric(nrep)

ac2_a=numeric(nrep)
ac2_b=numeric(nrep)
ac2_c=numeric(nrep)

for(n in 1:ii)
{
ar1=as.numeric(!chains1[[n]][1:(3500-1)]-chains1[[n]][2:3500]==0)
ar2=as.numeric(!chains2[[n]][1:(3500-1)]-chains2[[n]][2:3500]==0)
# first part
ac1_a[n]=mean(ar1[1:500])
ac2_a[n]=mean(ar2[1:1500])
# second part
ac1_b[n]=mean(ar1[501:1500])
ac2_b[n]=mean(ar2[501:1500])
# last part
ac1_c[n]=mean(ar1[1501:(3000-1)])
ac2_c[n]=mean(ar2[1501:(3500-1)])
}


boxplot(ac2_a)
boxplot(ac2_b)
boxplot(ac2_c)


