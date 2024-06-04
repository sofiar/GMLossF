#######################################################
################## Real data analysis #################
#######################################################
library(readr)
library(ggplot2)
library(PVAClone)

source('Extrafunctions.R')
source('mcmc_alg.R')

#load data 
American_redStart <- read_csv("./Real_Data_analysis/American_redStart2.csv")
ggplot(American_redStart)+geom_line(aes(years,redstart))+theme_bw()

# fit model 

phi.1=0.5
phi.2=0.5 
eta.1=0
eta.2=10

#phi.1=theta2^2/1+2, phi.2=(theta2^2+1)*theta2, eta.1=theta1, eta.2=1

results=mcmc.quasiGibbs(iters=20000, burn=1000, n.chains=2, #theta1.init, theta2.init,
                  b.init=c(-0.5,-0.2),
                  # Observations
                  American_redStart$redstart,
                  # prior parameters
                  phi.1=phi.1, phi.2=phi.2, eta.1=eta.1, eta.2=eta.2)

plot(results$Keep.theta1[5000:6000,1],type='l',col='blue')
lines(results$Keep.theta1[,2],type='l',col='red')

plot(results$Keep.theta2[,1],type='l',col='blue')
lines(results$Keep.theta2[,2],type='l',col='red')

plot(results$Keep.b[5000:6000,1],type='l',col='blue')
lines(results$Keep.b[5000:6000,2],type='l',col='red')

mean(results$Keep.theta1[5000:6000,])
sqrt(mean(results$Keep.theta2[3000:10000,]))
mean(results$Keep.b[3000:10000,])




save.image("./realDataResults.RData")



