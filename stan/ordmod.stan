functions {
  vector make_cutpoints(vector p, real scale) {
    int C = rows(p) - 1;
    vector[C] cutpoints;
    real running_sum = 0;
    for(c in 1:C) {
      running_sum += p[c];
      cutpoints[c] = logit(running_sum);
    }
    return scale * cutpoints;
  }
}

data {
  int N; // number of records
  int K; // number of response levels
  int P; // number of covariates
  int<lower=1,upper=K> y[N];
  matrix[N, P] X; // design matrix
  vector[K] prior_counts; // prior for Dirichlet on cuts
  vector[P] prior_sd; // prior SD for beta coefficients
}

parameters {
  simplex[K] pi;    // category probabilities
  vector[P] beta_raw;   // covariate coefficients
}

transformed parameters {
  vector[K-1] c; // cut-points on unconstrained scale
  vector[P] beta = prior_sd .* beta_raw;
  c = make_cutpoints(pi, 1);
}

model {
  // Linear predictor
  vector[N] eta = X * beta;

  // Prior model
  pi ~ dirichlet(prior_counts);
  beta_raw ~ normal(0, 1);

  // Observational model
  y ~ ordered_logistic(eta, c);
}

generated quantities {
  int<lower=1, upper=K> y_ppc[N];
  // posterior predictive draws
  for (n in 1:N)
    y_ppc[n] = ordered_logistic_rng(X[n]*beta, c);
}
