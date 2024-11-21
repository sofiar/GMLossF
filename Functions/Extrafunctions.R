################################################################################
############################ Extra functions ###################################
################################################################################

################################################################################
lambertW_expArg = function(x) {

  ctrl = x <= 709
  res = double(length = length(x))

  if (sum(ctrl) > 0) res[ctrl] = lamW::lambertW0(exp(x[ctrl]))
  if (sum(!ctrl) > 0) {
    res[!ctrl] = x[!ctrl] - log(x[!ctrl]) + log(x[!ctrl]) / x[!ctrl] +
      (log(x[!ctrl]) * (log(x[!ctrl]) - 2)) / (2 * x[!ctrl]^2) +
      log(x[!ctrl]) * (2 * log(x[!ctrl])^2 - 9 * log(x[!ctrl]) + 6) /
        (6 * x[!ctrl]^3) + log(x[!ctrl]) *
      (
       -12 + 36 * log(x[!ctrl]) - 22 * log(x[!ctrl])^2 +
         3 * log(x[!ctrl])^3) / (12 * x[!ctrl]^4
    )
  }

  return(res)

}
################################################################################



################################################################################
GompPois_simulate = function(T, theta1, theta2, b) {

  # This function simulates Gompertz model by transformed parameters

  # T is the length of the chain

  res = double(T)
  a = -b * theta1
  sigmaSq = -b * (2 + b) * theta2

  res[1] = theta1 + sqrt(theta2) * rnorm(1)
  for (t in 2:T) {
    res[t] = a + (1 + b) * res[t - 1] + sqrt(sigmaSq) * rnorm(1)
  }

  res = rpois(T, lambda = exp(res))

  return(res)

}
################################################################################



################################################################################
Simulate_Gompertz=function(TT,a,b,sigma2){
  
  # This function simulates Gompertz model by the original paramters 

  theta1=-a/b
  theta2=-sigma2/(b*(b+2))

  Nt=numeric(TT)
  Nt[1]=exp(rnorm(1,theta1,sqrt(theta2)))
  Es=rnorm(TT,0,sqrt(sigma2))
  
  for(t in 2:TT){
    Nt[t]=Nt[t-1]*exp(a+b*log(Nt[t-1])+Es[(t-1)])  
  }

  Nt.obs=numeric(TT)

  for(t in 1:TT){
    Nt.obs[t]=rpois(1,lambda=Nt[t])
  }

return(Nt.obs)

}
################################################################################



################################################################################
GompPois_simulate_ln = function(TT,a,b,sigma2,tau2) {
  
  # Simulation Gompertz with log normal error 

  Y = numeric(TT)
  mu = -a/b
  var = sigma2/(1-((1+b)^2))+tau2
  Y[1] = rnorm(1,mu,sqrt(var))
  
  for(i in 2:length(Y)){
    mu = a+(1+b)*(mu+(var-tau2)/var*(Y[i-1]-mu))
    var = (1+b)^2*(var-tau2)/var*tau2+sigma2+tau2
    Y[i] = rnorm(1,mu,sqrt(var))
  } 

  return(exp(Y))

}
################################################################################



################################################################################
GompPois_MoM = function(Nstar) {

  # Gompertz Model with Poisson Sampling Error: Method of Moments
  #
  # the parameters of the model are chosen in order to match the empirical mean,
  #  variance and autocovariance at lag 1.

  t1 = mean(Nstar)
  t2prime = var(Nstar) - t1
  t3 = acf(Nstar, lag.max = 1, plot = FALSE, type = "covariance")$acf[2]

  theta2hat = log(1 + t2prime / t1^2)
  theta1hat = log(t1) - 0.5 * theta2hat
  theta3hat = log(1 + t3 * exp(-2 * theta1hat - theta2hat))

  b =  theta3hat / theta2hat - 1
  a = -theta1hat * b
  sigmaSq = -theta2hat * b * (2 + b)

  return(list(
    theta1 = theta1hat,
    theta2 = theta2hat,
    b = b,
    a = a,
    sigmaSq = sigmaSq
  ))

}
################################################################################



################################################################################
# old versio of approx.CL2
GompPois_composite_likelihood=function(Nstar,theta1,theta2,b,nsim,
                                      d=NULL){
  
  # Computes the approximate composite likelihood considering all posible   
  # combinations of points with d distance
 
  # d<length(Nstar)

  curr=0
  Times=length(Nstar)
  mu=rep(theta1, 2)

  if(is.null(d)){d=Times}
  if(d>Times){cat('d must be less than the length of the time series')}
  
  for (t1 in 1:(Times-1)) {
    limit=min(t1+d,Times)    
    
    for (t2 in (t1+1):limit) {
      B = (1+b)^(abs(outer(c(1,abs(t2-t1+1)),c(1,abs(t2-t1+1)) , "-"))) 
      cov_matrix=theta2*B
      N.sample = exp(MASS::mvrnorm(nsim, mu = mu, Sigma = cov_matrix))
      logval = colSums(dpois(Nstar[c(t1,t2)],t(N.sample),log=TRUE))
      curr = curr-log(nsim) + log( sum(exp(logval)))
    }    
  
  }
 
  return(curr)
}
################################################################################



################################################################################
log_targetPoissonGauss = function(x,n,tauSq,mu){
  -exp(x) + x*n - 0.5/tauSq * (x-mu)^2
}
################################################################################



################################################################################
log_proposalGauss = function(x,tauSq,xi){
  - 0.5/tauSq * (x-xi)^2
}
################################################################################





