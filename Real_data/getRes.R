setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("./RealData_AmericanRedstart2.RData")



res_mle$par
rowMeans(res_Gibbs)
rowMeans(res_stan)

time_mle
time_Gibbs
time_stan

sns::ess(t(res_Gibbs))
sns::ess(t(res_stan))

sns::ess(t(res_Gibbs)) / as.numeric(time_Gibbs)
sns::ess(t(res_stan)) / as.numeric(time_stan)

par(mfrow = c(2, 3))

acf(res_Gibbs["theta1", ])
acf(res_Gibbs["theta2", ])
acf(res_Gibbs["b", ])

acf(res_stan["theta1", ])
acf(res_stan["theta2", ])
acf(res_stan["b", ])
