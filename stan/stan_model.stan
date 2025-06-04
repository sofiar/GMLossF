data {
  int<lower=0> M;
  int Nstar[M];
}

parameters {
  real theta1;
  real b;
  real<lower=0> theta2;
  vector<lower=0>[M] N;
}

transformed parameters {
  vector[M] mu; // Mean vector
  cov_matrix[M] Sigma; // Covariance matrix
  
  mu = rep_vector(theta1, M);
 
  // Compute Sigma
  for (i in 1:M) {
    for (j in 1:M) {
      Sigma[i, j] = (theta2) * pow(1 + b, abs(i - j));
    }
  }
  
}
  
model {
  // Priors
  theta2 ~ inv_gamma(0.1, 0.1); 
  theta1 ~ normal(0, 10*sqrt(theta2));
  b ~ uniform(-2, 0);
    
  // Likelihood
  log(N[1]) ~ normal(theta1, sqrt(theta2));
  
  for (t in 2:M) {
    log(N[t]) ~ normal(-b*theta1 + (1+b)*log(N[t-1]), sqrt(-b*(2+b)*theta2));
  }
  
  for(i in 1:M){
    Nstar[i] ~ poisson(N[i]);
  }
  
  }
