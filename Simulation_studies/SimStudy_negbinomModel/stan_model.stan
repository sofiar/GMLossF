data {
  int<lower=0> M;
  int Nstar[M];
}

parameters {
  real theta1;
  real<lower=-2, upper=0> b;
  real<lower=0> theta2;
  vector<lower=0>[M] N;
}
  
model {
  // Priors
  theta2 ~ inv_gamma(0.1, 0.1); 
  theta1 ~ normal(0, 10*sqrt(theta2));
  b ~ uniform(-2, 0);
    
  // Likelihood
  N[1] ~ lognormal(theta1, sqrt(theta2));
  
  for (t in 2:M) {
    N[t] ~ lognormal(-b*theta1 + (1+b)*log(N[t-1]), sqrt(-b*(2+b)*theta2));
  }
  
  for(i in 1:M){
    Nstar[i] ~ poisson(N[i]);
  }
  
  }
