load('./compute_likelihoodV0.RData')

GompUnkGP_likelihood0 = function(Zstar, theta1, theta2, b, M = 1e+4, log = FALSE){
    eZstar = exp(Zstar)
    a = -b*theta1
    sigma =sqrt(-b*(2+b))*sqrt(theta2)
    TT = length(Zstar)

    Zs = matrix(data = NA, nrow = TT, ncol = M)
    Zs[1,] = rnorm(n = M, mean = theta1, sd = sqrt.theta2)
    for(t in 2:TT){
     Zs[t,] = a + (1+b)*Zs[t-1,] + rnorm(M, mean = 0, sd = sigma)
    }
    
    logweights = colSums(dpois(eZstar,exp(Zs),log=TRUE))
    maxLogweights = max(logweights)

    logres = -log(M) + maxLogweights + log( sum(exp(logweights-maxLogweights)) )
    if(log){return(logres)}else{return(exp(logres))}
}

theta1=1.9244
sqrt.theta2=0.4726
b=-.8
TT = 1000
a = -b*theta1
sigma =sqrt(-b*(2+b))*sqrt.theta2


B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt.theta2*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )

bs=-c(1:1000)/(1001) 
ys=numeric(length(bs))

for(i in 1:length(bs)){
  ys[i] = GompUnkGP_likelihood0(Zstar, theta1, theta2 = sqrt.theta2^2, b = bs[i], M = 1e+4, log = TRUE)
}

#######################################################
theta1=1.9244
sqrt.theta2=0.4726
b=-.24
TT = 1000
a = -b*theta1
sigma =sqrt(-b*(2+b))*sqrt.theta2


B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt.theta2*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )

bs=-c(1:1000)/(1001) 
ys_2=numeric(length(bs))

for(i in 1:length(bs)){
  ys_2[i] = GompUnkGP_likelihood0(Zstar, theta1, theta2 = sqrt.theta2^2, b = bs[i], M = 1e+4, log = TRUE)
}


#######################################################
theta1=1.9244
sqrt.theta2=0.4726
b=-.5
TT = 1000
a = -b*theta1
sigma =sqrt(-b*(2+b))*sqrt.theta2


B = (1+b)^(abs(outer(1:TT, 1:TT, "-"))) # get B matrix
Z=theta1+sqrt.theta2*t(chol(B))%*%rnorm(TT)
Zstar = log( rpois(TT, exp(Z)) )

bs=-c(1:1000)/(1001) 
ys_3=numeric(length(bs))

for(i in 1:length(bs)){
  ys_3[i] = GompUnkGP_likelihood0(Zstar, theta1, theta2 = sqrt.theta2^2, b = bs[i], M = 1e+4, log = TRUE)
}



# lines(bs,ys,col='blue')
# abline(v=-0.8, col = "blue", lty  = 2)
# lines(bs,ys_2,col='red')
# abline(v=-0.24, col = "red", lty  = 2)
# lines(bs,ys_3,col='green3')
# abline(v=-0.5, col = "green3", lty  = 2)


save.image('./compute_likelihoodV0.RData')

