### linear model with interactions to check robustness to multi-modality

set.seed(42)

require(rstan)

d <- 15
n <- 100

X <- data.frame(matrix(rnorm(d*n),n,d))

designMat <- model.matrix(~ . + .^2 + 0,data = X)

p <- ncol(designMat)

beta0_true <- rnorm(1)

beta_true <- rnorm(p)

beta_true[sample.int(p,p-10)] <- 0

y <- as.numeric(beta0_true + designMat%*%beta_true + 0.1*rnorm(n))

model_code <- "data{
int<lower=1> N;
int<lower=1> ncov;
real y[N];
matrix[N,ncov] x;
real<lower = 0> sigma_indic;
real mu_indic;
real<lower = 0> tau;
}

parameters{
vector[ncov] beta_raw;
vector[ncov] indic_raw;
real beta_0;
real<lower = 0> sigma_noise;
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

y ~ normal(beta_0 + x*beta,sigma_noise);

}
"

sfit <- stan(model_code = model_code,data = list(N = 100,
                                         ncov = 120,
                                         y = y,
                                         x = designMat,
                                         sigma_indic = 10,
                                         mu_indic = 0,
                                         tau = 5),
              chains = 4,iter = 4000, control = list(adapt_delta = 0.99))

require(shinystan)

launch_shinystan(sfit) ## see parameters beta[13],[28],[70],[87],[106],[113],[115],[119]

