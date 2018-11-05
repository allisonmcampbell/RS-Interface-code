data{
  int<lower=1> N;
  int<lower = 1> N_test;
  int<lower=1> ncov;
  int<lower = 1> K;
  int y[N];
  matrix[N,(ncov*K)] PHI;
  matrix[N_test,(ncov*K)] phi_test;
  real<lower = 0> sigma_indic;
  real mu_indic;
  real<lower = 0> tau;
}

parameters{
  vector[ncov*K] beta_raw;
  vector[ncov*K] indic_raw;
  real beta_0;
}

transformed parameters{
  vector[(ncov*K)] beta;
  vector<lower=0,upper=1>[(ncov*K)] indic;
  
  indic = inv_logit(mu_indic + sigma_indic*indic_raw);
  
  for (j in 1:ncov){
  beta[(1 + K*(j-1))] = tau*indic[(1 + K*(j-1))] .* beta_raw[(1 + K*(j-1))];
  }
  
  for (j in 1:ncov){
    for (i in 2:K){
      beta[(K*(j-1) + i)] = tau * indic[(1 + K*(j-1))] * indic[(K*(j-1) + i)] * beta_raw[(K*(j-1) + i)];
    }
  }
  
}

model{
  beta_raw ~ normal(0,1);
  
  indic_raw ~ normal(0,1);
  
  y ~ bernoulli_logit(beta_0 + PHI*beta);
  
}

generated quantities{
  int y_rep[N_test];
  
  for (i in 1:N_test){
    y_rep[i] = bernoulli_logit_rng(beta_0 + phi_test[i,]*beta);
  }
}
