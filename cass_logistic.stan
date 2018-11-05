data{
  int<lower=1> N;
  int<lower=1> ncov;
  int y[N];
  matrix[N,ncov] x;
  real<lower = 0> sigma_indic;
  real mu_indic;
  real<lower = 0> tau;
}

parameters{
  vector[ncov] beta_raw;
  vector[ncov] indic_raw;
  real beta_0;
}

transformed parameters{
  vector[ncov] beta;
  vector<lower=0,upper=1>[ncov] indic;
  
  indic = inv_logit(mu_indic + sigma_indic*indic_raw);
  
  beta = tau * indic .* beta_raw;
  
}

model{
  beta_raw ~ normal(0,1);
  
  indic_raw ~ normal(0,1);
  
  y ~ bernoulli_logit(beta_0 + x*beta);
  
}
