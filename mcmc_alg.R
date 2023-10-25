#### MCMC
source('Extrafunctions.R')
source('GompUnkGP_likelihood.R')

###############################################################################################
## This function computes the alg considering the likelihood computed via important sampling.##
## parameters: theta1,theta2,b                                                               ##
###############################################################################################

mcmc.new=function(iters=3500, burn=1000, n.chains=2, theta1.init, theta2.init,
beta.init,
 # Observations
 Nt.obs,
 # prior parameters
 phi.1, phi.2, eta.1, eta.2)

{
 keep.theta1=array(0, dim = c(iters,n.chains))
 keep.theta2=array(0, dim = c(iters,n.chains))
 keep.b=array(0, dim = c(iters,n.chains)) 

 start=proc.time()[3]
 cat(" MCMC go!.\n")
   for (c in 1:n.chains)
  {
    # initialize parameter values
    Zt=Nt.obs
    TT=length(Nt.obs)
    Theta1=theta1.init[c]
    Theta2=theta2.init[c]
    Beta=beta.init[c]
    b=1/(1+exp(Beta))
    B= (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
   
   for (i in 1:iters)
   {
    #########################################################
    ### 1. Sample beta through elliptical slice sampling ####
    #########################################################
    
    # update b and B!
    
    #########################################################  
    ###         2. Sample theta2 (Inverse Gamma)         ####  
    #########################################################
    shape.t2 = shape.theta2 + TT/2
    #inv.m= solve(eta.2*diag(TT)+B) 
    ES1 = eigen(eta.2*diag(TT)+B)
    S1_G = ES1$vectors
    S1_D = ES1$values
    inv.m = S1_G%*%diag(1/S1_D)%*%t(S1_G) 
    rate.t2 = rate.theta2 +1/2 * t(Zt-eta.1)%*%inv.m%*%(Zt-eta.1)
    Theta2=1/rgamma(1, shape.t2, rate.t2)
    keep.theta2[i,c]=Theta2

    #########################################################
    ###           3. Sample theta1 (Gaussian)             ###
    #########################################################
    
    ES = eigen(B)
    S1_E = ES$vectors
    S1_V = ES$values
    inv.B = S1_E%*%diag(1/S1_V)%*%t(S1_E) 
    OnesV=rep(1,TT)
    pp= eta.2*t(OnesV)%*%inv.B
    mean=(pp%*%Zt+eta.1)/(pp%*%OnesV+1)
    var=Theta2*eta.2/(eta.2*t(OnesV)%*%inv.B%*%OnesV)
    Theta1=rnorm(1, mean, sd=sqrt(var))
    keep.theta1[i,c]=Theta1

    #########################################################
    ###        4.Sample V acceptance-rejection method     ###
    #########################################################

    #curll=GompUnkGP_likelihood(Theta1, Theta2, b, Nt.obs, nsim, log =TRUE)  
    
    #########################################################
    ###     5. Sample Z=log(Nt) Acceptance-Rejection Al   ###
    #########################################################





     
   
   
   }
  
  
  }





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

  



