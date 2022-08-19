data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
// switch parameters and transformed parameters
parameters {
  vector[K] b;  // population-level effects
  real<lower=0,upper=1> p0;
}
transformed parameters {
  real b0;
  b0 = logit(p0);
}
model {
  if (!prior_only) {
    vector[N] mu = X * b;
    for (i in 1:N) {
      target += bernoulli_logit_lpmf(Y[i] | b0 + mu[i]);
    }
  }
  target += uniform_lpdf(p0 | 0, 1);
  target += normal_lpdf(b | 0, 10);
}
