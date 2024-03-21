
mcmc.quasiGibbs_gprior=function(iters=3500, burn=1000, n.chains=2, #theta1.init, theta2.init,
                  b.init,
                  # Observations
                  Nt.obs,
                  # prior parameters
                  phi.1, phi.2, eta.1,c=1)
  
{
  keep.theta1 = array(0, dim = c(iters,n.chains))
  keep.theta2 = array(0, dim = c(iters,n.chains))
  keep.b = array(0, dim = c(iters,n.chains)) 
  
  start=proc.time()[3]
  cat(as.character(Sys.time())," MCMC go!.\n")
  V = pi^2 / 3 # initialize V equal to its prior mean
  TT = length(Nt.obs)
  shape.t2 = phi.1 + TT/2
  eta.1T=rep(eta.1,TT)

  d1 = -(2*phi.1+TT)/2
  d2 = 1/(2*phi.2)
  
  for (c in 1:n.chains)
  {
    # initialize parameter values
    Zt = log(Nt.obs)
    if(sum(is.infinite(Zt))>0){Zt[is.infinite(Zt)]= 0}
    b = b.init[c]
    B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
    beta = log( -b / (1+b) )
    
    for (i in 1:iters)
    {
      #########################################################
      ### 1. Sample beta through elliptical slice sampling ####
      #########################################################
      
      # update log-likelihood of current value of b
      eta.2 = c(1/TT * t(rep(1,TT))%*%B%*%(rep(1,TT)))
      ES1 = chol( eta.2*matrix(1, nrow = TT, ncol = TT) + B )
      loglike_Z = -sum( log(diag(ES1)) )
      ES1 = backsolve( ES1, diag(1, nrow = TT) )
      QF = tcrossprod( t(Zt-eta.1T)%*%ES1 )
      loglike_Z = loglike_Z + d1*log( 1+d2*QF )
      
      # update b, B and beta!
      #y0 = loglike_Z
      #uu = runif(n=1,min=0,max=1)
      y = loglike_Z + log(runif(n=1,min=0,max=1))
      delta = runif(n=1,min=0,max=2*pi)
      delta_min = delta-2*pi
      delta_max = delta
      betaStar = rnorm(1,0,sqrt(c*V))# added c here
     #cc = 0
      while(TRUE)
      {
       #cc = cc + 1
        betaProp = beta*cos(delta) + betaStar*sin(delta)
        b = -1/(1+exp(-betaProp))
        B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
        eta.2 = c(1/TT * t(rep(1,TT))%*%B%*%(rep(1,TT)))
        ES1 = chol( eta.2*matrix(1, nrow = TT, ncol = TT) + B )
        loglike_Z = -sum( log(diag(ES1)) )
        ES1 = backsolve( ES1, diag(1, nrow = TT) )
        QF = tcrossprod( t(Zt-eta.1T)%*%ES1 )
        loglike_Z = loglike_Z + d1*log( 1+d2*QF )
        #print(paste(betaProp, b, loglike_Z, sep = '  '))
        if(any(loglike_Z > y, delta_max - delta_min < sqrt(.Machine$double.eps)))
        {
          beta = betaProp
          keep.b[i,c] = b
          break
        }
        else
        {
          if(delta<=0){delta_min=delta}else{delta_max=delta}
          delta = runif(n=1,min=delta_min,max=delta_max)
        }
        #print(paste(cc, abs(delta_min-delta_max), sep=' '))
      }
      
      #########################################################  
      ###         2. Sample theta2 (Inverse Gamma)         ####  
      #########################################################
      
      #ES1 = chol(eta.2TT+B)
      #ES1 = backsolve(ES1, diag(1, nrow = TT)) 
      rate.t2=phi.2+1/2*QF
      Theta2 = 1/rgamma(1, shape=shape.t2, rate=rate.t2)
      keep.theta2[i,c] = Theta2
            
      #########################################################
      ###           3. Sample theta1 (Gaussian)             ###
      #########################################################
      
      inv.B = chol2inv(chol(B))
      pp = rep(eta.2, TT)%*%inv.B
      mean.t1 = pp%*%Zt+eta.1
      var.t1 = Theta2*eta.2
      pp = pp%*%rep(1,TT) + 1
      mean.t1 = mean.t1/pp
      var.t1 = var.t1/pp
      Theta1 = rnorm(1, mean.t1, sd=sqrt(var.t1))
      keep.theta1[i,c] = Theta1
      
      #########################################################
      ###  4. Sample V through acceptance-rejection method  ###
      #########################################################
      
      V = rpostlogiskolmo(n = 1, x = beta/sqrt(c))# added sqrt(c) here
      
      #########################################################
      ###     5. Sample Z=log(Nt) Acceptance-Rejection Al   ###
      #########################################################
      
      # update sigmaSq and a
      sigmaSq=-Theta2*b*(2+b)
      a=-b*Theta1

      for (t in 1:TT)
      {
      if(t==1)
      {
      mu=(Theta1*sigmaSq+(1+b)*(Zt[2]-a)*Theta2)/(sigmaSq+(1+b)^2*Theta2)
      tauSq=sigmaSq*Theta2/(sigmaSq+(1+b)^2*Theta2)
      }else if (t==TT) {
      mu=a+(1+b)*Zt[TT-1]
      tauSq=sigmaSq/(1+(1+b)^2)
      }else {
      mu=(a+(1+b)*(Zt[t+1]+Zt[t-1]-a))/(1+(1+b)^2)    
      tauSq=sigmaSq
      }

      xi = Nt.obs[t]*tauSq + mu - lambertW_expArg(log(tauSq) + Nt.obs[t]*tauSq + mu)
      logC = log_targetPoissonGauss(x=xi,n=Nt.obs[t],tauSq=tauSq,mu=mu) - log_proposalGauss(x=xi,tauSq=tauSq,xi=xi)
      while(TRUE)
      {
      x = rnorm(100, mean = xi, sd = sqrt(tauSq)) 
      u = runif(100, 0, 1)
      check = log(u) <= log_targetPoissonGauss(x=x,n=Nt.obs[t],tauSq=tauSq,mu=mu) - log_proposalGauss(x=x,tauSq=tauSq,xi=xi) - logC
      if(sum(check)>0)
      {
        Zt[t] = x[check][1]
        break
      }
      }
      
      }
      
      #print(i)
      if(i%%1000 == 0){
        cat(as.character(Sys.time()), " ", i, '\n')
      }
    }
    
    
  }
end=proc.time()[3]
   out=list(Keep.theta1=keep.theta1, Keep.theta2=keep.theta2,
           Keep.b=keep.b,time=end-start)
  return(out)
}


