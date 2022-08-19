data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
  real b0;
}
transformed parameters {
  real p0;
  p0 = inv_logit(b0);
}
model {
  if (!prior_only) {
    vector[N] mu = X * b;
    for (i in 1:N) {
      target += bernoulli_logit_lpmf(Y[i] | b0 + mu[i]);
    }
  }
  // prior on non-linear transformed parameter, needs Jacobian
  target += uniform_lpdf(p0 | 0, 1);
  target += normal_lpdf(b | 0, 10);
}
