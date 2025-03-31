data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> J;
  matrix[N, K] X;
  matrix[N, J] Y;
}

parameters {
  vector[J] beta0;
  matrix[K, J] beta;
  vector<lower=0>[J] sigma;
  cholesky_factor_corr[J] L_Omega;
}

model {
  beta0 ~ normal(0, 10);
  to_vector(beta) ~ normal(0, 5);
  sigma ~ cauchy(0, 2.5);
  L_Omega ~ lkj_corr_cholesky(2);

  for (n in 1:N) {
    vector[J] mu = beta0 + to_vector(X[n] * beta);
    Y[n] ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L_Omega));
  }
}

generated quantities {
  matrix[N, J] y_rep;

  for (n in 1:N) {
    vector[J] mu = beta0 + to_vector(X[n] * beta);
    y_rep[n] = multi_normal_cholesky_rng(mu, diag_pre_multiply(sigma, L_Omega));
  }

  vector[N] log_lik;
  for (n in 1:N) {
    vector[J] mu = beta0 + to_vector(X[n] * beta);
    log_lik[n] = multi_normal_cholesky_lpdf(Y[n] | mu, diag_pre_multiply(sigma, L_Omega));
  }

}
