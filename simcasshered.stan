data{
  int<lower=1> N;
  int<lower=1> ncov;
  vector[N] y;
  matrix[N,ncov] X;
  real<lower = 0> sigma_indic;
  real mu_indic;
  int<lower = 0> ngroup;
}

parameters{
  vector[ncov] beta_raw;
  vector[ncov] indic_raw;
  real beta_0;
  real<lower = 0> sigma;
  vector [ngroup] beta_pop_raw;
  vector[ngroup] indic_pop_raw;
}

transformed parameters{
  vector[ncov] beta;
  vector<lower=0,upper=1>[ncov] indic;
  vector[ngroup] beta_pop;
  vector[ngroup] indic_pop;
  
  indic = inv_logit(mu_indic + sigma_indic*indic_raw);
  indic_pop = inv_logit(mu_indic + sigma_indic * indic_pop_raw);
  
  beta_pop = 5 * indic_pop .* beta_pop_raw;
  
  for (i in 1:ngroup){
  
  for (j in 1:5){
  beta[5*(i-1) + j] = beta_pop[i] + 5 * indic [5*(i-1) + j] * indic_pop[i] * beta_raw[5*(i-1) + j];
  }
  
  }
  
}

model{
  
  beta_raw ~ normal(0,1);
  beta_pop_raw ~ normal(0,1);
  
  indic_raw ~ normal(0,1);
  indic_pop_raw ~ normal(0,1);
  
  sigma ~ cauchy(0,1);
  
  y ~ normal(beta_0 + X*beta,sigma);
  
}
