// survreg_ph_mspline.stan

data {
  int<lower=1> N; // number of observations
  int<lower=1> M; // number of mSpline knots
  int<lower=1> Mb;
  array[N] real y;      // event times
  array[N] int v;       // censoring indicator
  vector[N] x;       // treatment - 1 or 0
  matrix[N, M] Z; // mSpline matrix (with intercept included)
  matrix[N, M] W; // iSpline matrix
  matrix[N, Mb] B; // bSpline matrix
  vector[M] alpha; // Dirichlet hyperparameter
  real mu0;
  int<lower=0, upper=1> estimate; // estimate the model or just simulate data
}

parameters {
  // lp intercept
  real beta0;
  // b-spline intercept
  real theta0;
  simplex[M] gamma;
  vector[Mb] epsilon;
  real<lower=0> tau;
}

transformed parameters {
  vector[Mb] theta;
  vector[N] beta;
  theta[1] = epsilon[1];
  for(m in 2:Mb) {
    theta[m] = theta[m-1] + tau*epsilon[m];
  }
  beta = theta0 + B * theta;
}

model {
  vector[N] mu;
  vector[N] z;
  vector[N] w;

  // baseline hazard/cumulative hazard
  z = log(Z*gamma);
  w = W*gamma;

  // prior
  // random walk prior on B-spline coefficients
  epsilon[1] ~ normal(0, 10);
  epsilon[2:Mb] ~ normal(0, 1);
  tau ~ cauchy(0, 1);
  // priors on intercepts
  theta0 ~ normal(0, 10);
  beta0 ~ normal(mu0, 10);
  // prior on M-spline coefficients
  gamma ~ dirichlet(alpha);

  // likelihood
  if(estimate) {
    for(n in 1:N) {
      mu[n] = beta0 + x[n] * beta[n];
      if (v[n] == 1)
        target += mu[n] + z[n];
      target += -exp(mu[n]) * w[n];
    }
  }
}
