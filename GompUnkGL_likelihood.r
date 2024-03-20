### Compute likelihood for the log normal model ####

##########################################################################################
#
# Gompertz Model with Unknown Response II Version: Gaussian White
#  Noise and log Normal Sampling Distribution
# 
##########################################################################################

############### function
GompUnkGL_likelihood = function(
    a, b,sigma2,tau2, Nstar,
    log = FALSE){
  
  # a whose likelihood must be computed;
  # sigma2 whose likelihood must be computed;
  # tau2 whose likelihood must be computed; 
  # b whose likelihood must be computed;
  # Nstar is the vector of the observed variables, say Nstar_t, t = 1,2,...,T;
  # log is a logical variable. Should the log-likelihood be computed?

  ###############--------------- function script 
  
  ########## initialization
  
  ### sample size
  T = length(Nstar)

  Y=log(Nstar)

  mu = -a/b
  var = sigma2/(1-((1+b)^2))+tau2

  logres = dnorm(Y[1],mu,sqrt(var),log= TRUE)

for (i in 2:T)
{ 
  mu = a+(1+b)*(mu+(var-tau2)/var*(Y[i-1]-mu))
  var = (1+b)^2*(var-tau2)/var*tau2+sigma2+tau2
  logres = logres + dnorm(Y[i],mu,sqrt(var),log= TRUE)
}
if(log)
{return(logres)
}else{
return(exp(logres))
}
}


############### appendix

if(FALSE){ # appendix script is not run
  
  # debugging
  rm(list = ls(all.names = TRUE)); gc(); cat("\14")
  set.seed(1994)
  a = 0
  sigma2 =0.24
  tau2=0.5
  b = -0.8
  Y = numeric(50)
  mu = -a/b
  var = sigma2/(1-((1+b)^2))+tau2
  Y[1] = rnorm(1,mu,sqrt(var))
  for(i in 2:length(Y))
  {
  mu = a+(1+b)*(mu+(var-tau2)/var*(Y[i-1]-mu))
  var = (1+b)^2*(var-tau2)/var*tau2+sigma2+tau2
  Y[i] = rnorm(1,mu,sqrt(var))
  } 
 
  Nstar = exp(Y)
  log = TRUE
  trial = GompUnkGL_likelihood(a, b, sigma2,tau2, Nstar, log)
  
   
  ###############
  
  
  # trials  
  
}
