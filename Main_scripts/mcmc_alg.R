#### MCMC
source('Extrafunctions.R')
#source('GompUnkGP_likelihood.R')
source('dkolmo.R')
source('rpostlogiskolmo.R')

###############################################################################################
#                                                                                             #
#                                                                                             #
###############################################################################################


mcmc.quasiGibbs_v2=function(iters=3500, burn=1000, n.chains=2, theta1.init, theta2.init,
                  b.init,
                  # Observations
                  Nt.obs,
                  # prior parameters
                  phi.1, phi.2, eta.1, eta.2,c=1)
  
{
  keep.theta1 = array(0, dim = c(iters,n.chains))
  keep.theta2 = array(0, dim = c(iters,n.chains))
  keep.b = array(0, dim = c(iters,n.chains)) 
  
  start=proc.time()[3]
  cat(as.character(Sys.time())," MCMC go!.\n")
  V = pi^2 / 3 # initialize V equal to its prior mean
  lambda = phi.2/phi.1 # initialize lambda equal to its prior mean
  TT = length(Nt.obs)
  #shape.t2 = phi.1 + TT/2
  shape.t2 = 1 + TT/2
  eta.1T=rep(eta.1,TT)
  eta.2TT =  eta.2*matrix(1, nrow = TT, ncol = TT)
  eta.2T = rep(eta.2, TT)
  # d1 = -(2*phi.1+TT)/2 # for inverse gamma prior theta2
  # d2 = 1/(2*phi.2)
  
  for (c in 1:n.chains)
  {
    # initialize parameter values
    Zt = log(Nt.obs)
    if(sum(is.infinite(Zt))>0){Zt[is.infinite(Zt)]= 0}
    b = b.init[c]
    theta1=theta1.init[c]
    theta2=theta2.init[c]
    B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
    beta = log( -b / (1+b) )
    
    for (i in 1:iters)
    {
      #########################################################
      ### 1. Sample beta through elliptical slice sampling ####
      #########################################################
      
      if(FALSE){
        # update log-likelihood of current value of b
      ES1=chol(B)
      loglike_Z=-sum(log(diag(ES1))) 
      ES1 = backsolve( ES1, diag(1, nrow = TT) )
      ES1 = tcrossprod( t(Zt-rep(theta1,TT))%*%ES1 )
      loglike_Z = loglike_Z -0.5*ES1/theta2
            # update b, B and beta!
      #y0 = loglike_Z
      #uu = runif(n=1,min=0,max=1)
      y = loglike_Z + log(runif(n=1,min=0,max=1))
      delta = runif(n=1,min=0,max=2*pi)
      delta_min = delta-2*pi
      delta_max = delta
      betaStar = rnorm(1,0,sqrt(c*V))
     #cc = 0
      while(TRUE)
      {
       #cc = cc + 1
        betaProp = beta*cos(delta) + betaStar*sin(delta)
        b = -1/(1+exp(-betaProp))
        B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
        ES1=chol(B)
        loglike_Z=-sum(log(diag(ES1))) 
        ES1 = backsolve( ES1, diag(1, nrow = TT) )
        ES1 = tcrossprod( t(Zt-rep(theta1,TT))%*%ES1 )
        loglike_Z = loglike_Z -0.5*ES1/theta2
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
      }
      
      #########################################################  
      ### 2. Sample theta2 (Conditional Inverse Gamma)     ####  
      #########################################################
      
      ES2 = backsolve(chol(eta.2TT+B), diag(1, nrow = TT))
      ES2 = tcrossprod( t(Zt-eta.1T)%*%ES2 )
      
      rate.t2=lambda+1/2*ES2      
      Theta2 = 1/rgamma(1, shape=shape.t2, rate=rate.t2)
      # original inverse gamma update
      #rate.t2=phi.2+1/2*QF
      #Theta2 = 1/rgamma(1, shape=shape.t2, rate=rate.t2)
      keep.theta2[i,c] = Theta2
            
      #########################################################
      ###           3. Sample theta1 (Gaussian)             ###
      #########################################################
      
      inv.B = chol2inv(chol(B))
      pp = eta.2T%*%inv.B
      mean.t1 = pp%*%Zt+eta.1
      var.t1 = Theta2*eta.2
      pp = pp%*%rep(1,TT) + 1
      mean.t1 = mean.t1/pp
      var.t1 = var.t1/pp
      Theta1 = rnorm(1, mean.t1, sd=sqrt(var.t1))
      keep.theta1[i,c] = Theta1
      
      #########################################################
      ###             4. Sample V and lambda                ###
      #########################################################
      
      V = rpostlogiskolmo(n = 1, x = beta/sqrt(c))
      lambda = rgamma(n=1,shape=phi.2+1,rate=1/theta2+phi.1)
      
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



###############################################################################################
#                                                                                             #
#                                                                                             #
###############################################################################################

mcmc.quasiGibbs=function(iters=3500, burn=1000, n.chains=2, #theta1.init, theta2.init,
                  b.init,
                  # Observations
                  Nt.obs,
                  # prior parameters
                  phi.1, phi.2, eta.1, eta.2,c=1)
  
{
  keep.theta1 = array(0, dim = c(iters,n.chains))
  keep.theta2 = array(0, dim = c(iters,n.chains))
  keep.b = array(0, dim = c(iters,n.chains)) 
  
  start=proc.time()[3]
  cat(as.character(Sys.time())," MCMC go!.\n")
  V = pi^2 / 3 # initialize V equal to its prior mean
  lambda = phi.2/phi.1 # initialize lambda equal to its prior mean
  TT = length(Nt.obs)
  #shape.t2 = phi.1 + TT/2
  shape.t2 = 1 + TT/2
  eta.1T=rep(eta.1,TT)
  eta.2TT =  eta.2*matrix(1, nrow = TT, ncol = TT)
  eta.2T = rep(eta.2, TT)
  # d1 = -(2*phi.1+TT)/2 # for inverse gamma prior theta2
  # d2 = 1/(2*phi.2)
  d1 = -(2+TT)/2
  d2 = 1/(2*lambda)
  
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
      
      for(ithin in 1:1e+4){

      # update log-likelihood of current value of b
      ES1 = chol(eta.2TT+B)
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
      betaStar = rnorm(1,0,sqrt(c*V))
     #cc = 0
      while(TRUE)
      {
       #cc = cc + 1
        betaProp = beta*cos(delta) + betaStar*sin(delta)
        b = -1/(1+exp(-betaProp))
        B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
        ES1 = chol(eta.2TT+B)
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

      }
      
      #########################################################  
      ###         2. Sample theta2 (Inverse Gamma)         ####  
      #########################################################
      
      #ES1 = chol(eta.2TT+B)
      #ES1 = backsolve(ES1, diag(1, nrow = TT)) 
      rate.t2=lambda+1/2*QF      
      Theta2 = 1/rgamma(1, shape=shape.t2, rate=rate.t2)
      # original inverse gamma update
      #rate.t2=phi.2+1/2*QF
      #Theta2 = 1/rgamma(1, shape=shape.t2, rate=rate.t2)
      keep.theta2[i,c] = Theta2
            
      #########################################################
      ###           3. Sample theta1 (Gaussian)             ###
      #########################################################
      
      inv.B = chol2inv(chol(B))
      pp = eta.2T%*%inv.B
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
      
      V = rpostlogiskolmo(n = 1, x = beta/sqrt(c))
      lambda = rgamma(n=1,shape=phi.2+1,rate=1/theta2+phi.1) 
      d2 = 1/(2*lambda)
      
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

mcmc.al1=function(iters=3500, burn=1000, n.chains=2, theta1t.init=1, theta2t.init=theta2,theta3t.init=theta3,
                  M=3500,
                 # Observations
                 Nt.obs,
                 # real parameters
                 theta1t=NULL,theta2t=NULL,theta3t=NULL,
                 # prior parameters
                 mu.theta1t=0, sd.theta1t=1, mu.theta2t=0, sd.theta2t=1,
                 # proposal parameters
                 p.theta1_1=0.1, p.theta2_1=0.05, 
		p.theta1_2=0.01, p.theta2_2=0.02)
{
  # Keep values
  keep.theta1t=array(0, dim = c(iters,n.chains))
  keep.theta2t=array(0, dim = c(iters,n.chains))
  keep.theta3t=array(0, dim = c(iters,n.chains))

  
  ac.theta1t=matrix(0,iters,n.chains)
  ac.theta2t=matrix(0,iters,n.chains)
  
  start=proc.time()[3]

  
  #### GO!
  for (c in 1:n.chains)
  {
    
    Theta1t=theta1t.init[c]
    Theta2t=theta2t.init[c]
    
    Theta1=exp(Theta1t)
    Theta2=exp(Theta2t)
    
    curll=bis.approx.CL1(Nt.obs,Theta1,Theta2,M)
    #curll=approx.CL1(Nt.obs,Theta1,Theta2,M)

    p.theta1=p.theta1_1
    p.theta2=p.theta2_1	
    for (i in 1:iters)
    {
       if(i<1000){
        p.theta1=p.theta1_1
        p.theta2=p.theta2_1
        m.theta1t=Theta1t
        m.theta2t=Theta2t
        
	
        }else if(i>=1000 & i<=1500){
	p.theta1=sd(keep.theta1t[600:(i-1),c])+0.001
	#m.theta1t=(mean(keep.theta1t[1:(i-1),c]))# log o exp here?
	m.theta1=Theta1t
	#if(is.na(p.theta1)){p.theta1=p.theta1_1}
	p.theta2=sd(keep.theta2t[100:(i-1),c])+0.001
        m.theta2t=Theta2t
	#m.theta2t=(mean(keep.theta2t[1:(i-1),c]))
	#if(is.na(p.theta2)){p.theta2=p.theta1_2}

	}else{
	
	 p.theta1=p.theta1_2
	 p.theta2=p.theta2_2

	m.theta1t=Theta1t
	m.theta2t=Theta2t

	}

    
      # Theta 1
      #cantheta1t = rnorm(1,m.theta1t,p.theta1)
      #cantheta1=exp(cantheta1t)
      #canll=bis.approx.CL1(Nt.obs,cantheta1,Theta2,M)
      ##canll= approx.CL1(Nt.obs,cantheta1,Theta2,M)
      
      #MH <- canll-curll+
      #  dnorm(cantheta1t,mu.theta1t,sd.theta1t,log=TRUE)-
      # dnorm(Theta1t,mu.theta1t,sd.theta1t,log=TRUE)
      
      #if(log(runif(1))<MH){
       # Theta1t  <- cantheta1t
	      #Theta1 <- exp(Theta1t)
        #curll  <- canll
        #ac.theta1t[i,c]=1
        #}
      
      # Theta 2
      cantheta2t = rnorm(1,m.theta2t,p.theta2)
      cantheta2=exp(cantheta2t)
      canll=bis.approx.CL1(Nt.obs,Theta1,cantheta2,M)
      #canll=approx.CL1(Nt.obs,Theta1,cantheta2,M)
      
      MH <- canll-curll+
      dnorm(cantheta2t,mu.theta2t,sd.theta2t,log=TRUE)-
      dnorm(Theta2t,mu.theta2t,sd.theta2t,log=TRUE)

      if(log(runif(1))<MH){
        Theta2t  <- cantheta2t
	      Theta2 <- exp(Theta2t)
        curll  <- canll
        ac.theta2t[i,c]=1
        }
      
      keep.theta1t[i,c]=Theta1t
      keep.theta2t[i,c]=Theta2t
    #cat('theta1t.\n')
     # print(Theta1t)
     # print(theta1t)
    cat('theta2t.\n') 
     print(Theta2t) 
     print(theta2t)
     print(i)		
    cat("acceptance ratio.\n")
    print(sum(ac.theta1t[,c])/(i))
    #print(sum(ac.theta2t)/(i))

	}
    
  }
  avaelapsed = proc.time()[3] - start  
  
  ## Acceptance ratio
  ar1=(sum(ac.theta1t))/(iters*n.chains)
  ar2=(sum(ac.theta2t))/(iters*n.chains)
  
  ## Check coverage 
  cover1=0
  if(!is.null(theta1t))
  {low=quantile(keep.theta1t[burn:iters,],probs = c(0.025))
  high=quantile(keep.theta1t[burn:iters,],probs = c(0.975))
  if (theta1t > low & theta1t < high) {cover1 = cover1 + 1}}
  
  # Theta2
  cover2=0
  if(!is.null(theta2t))
  {low=quantile(keep.theta2t[burn:iters,],probs = c(0.025))
  high=quantile(keep.theta2t[burn:iters,],probs = c(0.975))
  if (theta2t > low & theta2t < high) {cover2 = cover2 + 1}}
  
  out=list(time=avaelapsed,cover1=cover1,cover2=cover2,pseudo1=keep.theta1t,
                 pseudo2=keep.theta2t,ar1=ar1,ar2=ar2)
  return(out)
  
}

  

mcmc.al2=function(iters=3500, burn=1000, n.chains=2, theta1t.init=1, theta2t.init=theta2,
theta3t.init=theta3,
                  M=3500,
                 # Observations
                 Nt.obs,
                 # real parameters
                 theta1t=NULL,theta2t=NULL,theta3t=NULL,
                 # prior parameters
                 mu.theta1t=0, sd.theta1t=1, mu.theta2t=0, sd.theta2t=1, mu.theta3t=0, sd.theta3t=1,
                 # proposal parameters
                 p.theta1_1=0.1, p.theta2_1=0.05,p.theta3_1=0.05, 
		p.theta1_2=0.01, p.theta2_2=0.02,p.theta3_2=0.03)
{
  # Keep values
  keep.theta1t=array(0, dim = c(iters,n.chains))
  keep.theta2t=array(0, dim = c(iters,n.chains))
  keep.theta3t=array(0, dim = c(iters,n.chains))

  
  ac.theta1t=matrix(0,iters,n.chains)
  ac.theta2t=matrix(0,iters,n.chains)
  ac.theta3t=matrix(0,iters,n.chains)
  
  start=proc.time()[3]

  
  #### GO!
  for (c in 1:n.chains)
  {
    
    Theta1t=theta1t.init[c]
    Theta2t=theta2t.init[c]
    Theta3t=theta3t.init[c]
    
    Theta1=exp(Theta1t)
    Theta2=exp(Theta2t)
    Theta3=exp(Theta3t)
    B=-2*(Theta3)/(1+Theta3)


    #curll=bis.approx.CL2(Nt.obs,Theta1,Theta2,Theta3,M)
    curll=bis.approx.cpp.CL2(Nt.obs,Theta1,Theta2,B,M)

    p.theta1=p.theta1_1
    p.theta2=p.theta2_1	
   
    for (i in 1:iters)
    {
       if(i<500){
        p.theta1=p.theta1_1
        p.theta2=p.theta2_1
        m.theta1t=Theta1t
        m.theta2t=Theta2t
        }else if(i>=500 & i<=1500){
	p.theta1=sd(keep.theta1t[400:(i-1),c])+0.001
	.theta1t=(mean(keep.theta1t[1:(i-1),c]))# log o exp here?
	m.theta1=Theta1t
	if(is.na(p.theta1)){p.theta1=p.theta1_1}
	p.theta2=sd(keep.theta2t[100:(i-1),c])+0.001
        m.theta2t=Theta2t
	m.theta2t=(mean(keep.theta2t[1:(i-1),c]))
	if(is.na(p.theta2)){p.theta2=p.theta1_2}

	}else{
	
	 p.theta1=p.theta1_2
	 p.theta2=p.theta2_2

	m.theta1t=Theta1t
	m.theta2t=Theta2t

	}

    
      # Theta 1
     # cantheta1t = rnorm(1,m.theta1t,p.theta1)
     # cantheta1=exp(cantheta1t)
     # canll=bis.approx.cpp.CL2(Nt.obs,cantheta1,Theta2,B,M)
      
     # MH <- canll-curll+
     #   dnorm(cantheta1t,mu.theta1t,sd.theta1t,log=TRUE)-
     #  dnorm(Theta1t,mu.theta1t,sd.theta1t,log=TRUE)
     # 
     # if(log(runif(1))<MH){
     #  Theta1t  <- cantheta1t
	   #  Theta1 <- exp(Theta1t)
     #  curll  <- canll
     #  ac.theta1t[i,c]=1
     #   }
      
      # Theta 2
      cantheta2t = rnorm(1,m.theta2t,p.theta2)
      cantheta2=exp(cantheta2t)
      #canll=bis.approx.CL2(Nt.obs,Theta1,cantheta2,Theta3,M)
      canll=bis.approx.cpp.CL2(Nt.obs,Theta1,cantheta2,B,M)
      
      MH <- canll-curll+
      dnorm(cantheta2t,mu.theta2t,sd.theta2t,log=TRUE)-
      dnorm(Theta2t,mu.theta2t,sd.theta2t,log=TRUE)

      if(log(runif(1))<MH){
        Theta2t  <- cantheta2t
	      Theta2 <- exp(Theta2t)
        curll  <- canll
        ac.theta2t[i,c]=1
        }
      
      keep.theta1t[i,c]=Theta1t
      keep.theta2t[i,c]=Theta2t
    cat('theta1t.\n')
      print(Theta1t)
      print(theta1t)
    cat('theta2t.\n') 
     print(Theta2t) 
     print(theta2t)
     print(i)		
    #cat("acceptance ratio.\n")
    #print(sum(ac.theta1t[,c])/(i))
    #print(sum(ac.theta2t)/(i))

	}
    
  }
  avaelapsed = proc.time()[3] - start  
  
  ## Acceptance ratio
  ar1=(sum(ac.theta1t))/(iters*n.chains)
  ar2=(sum(ac.theta2t))/(iters*n.chains)
  
  ## Check coverage 
  cover1=0
  if(!is.null(theta1t))
  {low=quantile(keep.theta1t[burn:iters,],probs = c(0.025))
  high=quantile(keep.theta1t[burn:iters,],probs = c(0.975))
  if (theta1t > low & theta1t < high) {cover1 = cover1 + 1}}
  
  # Theta2
  cover2=0
  if(!is.null(theta2t))
  {low=quantile(keep.theta2t[burn:iters,],probs = c(0.025))
  high=quantile(keep.theta2t[burn:iters,],probs = c(0.975))
  if (theta2t > low & theta2t < high) {cover2 = cover2 + 1}}
  
  out=list(time=avaelapsed,cover1=cover1,cover2=cover2,pseudo1=keep.theta1t,
                 pseudo2=keep.theta2t,ar1=ar1,ar2=ar2)
  return(out)
  
}

  

