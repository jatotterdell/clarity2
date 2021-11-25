// James Totterdell
// Date: 2021-10-13
//
// if undertaking simulations we would generally not have
// any continuous covariates, but instead just discrete treatment groups
// Therefore, to save computation time we can aggregate over
// the covariate patterns and weight by counts.
// E.g. rather than analyse N = 2,100 participants 1:1:1 to 3 arms,
// We can analyse the 3 covariate patterns weighted by sample size.

functions {
  // calculate the multinomial coefficient
  // required if want to include constants in log likelihood evaluation
  // - y: observed counts
  real multinomial_coef(vector y) {
    return lgamma(sum(y) + 1) - sum(lgamma(y + 1));
  }

  // make_cutpoints
  // - p: outcome level probabilities
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

  // Pr(y == k) for k=1,...,K
  // - c: cut-points for outcome levels
  // - eta: linear predictor for each pattern
  matrix Pr(vector c, vector eta) {
    int N = num_elements(eta);
    int K = num_elements(c) + 1;
    matrix[N, K] out;
    // for stability, work on log-scale
    for (n in 1:N) {
      out[n, 1] = log1m_exp(-log1p_exp(-(eta[n] - c[1]))); // ln(1 - inv_logit(eta[n] - c[1]))
      out[n, K] = -log1p_exp(-(eta[n] - c[K-1]));          // ln(inv_logit(eta[n] - c[K-1]))
      for (k in 2:(K - 1)) {
        // ln(inv_logit(eta[n] - c[k-1]) - inv_logit(eta[n] - c[k]))
        out[n, k] = log_diff_exp(-log1p_exp(-(eta[n] - c[k-1])), -log1p_exp(-(eta[n] - c[k])));
      }
    }
    return exp(out);
  }

  // log-likelihood (multinomial)
  // - p: level probabilities for each pattern
  // - y: observed count for each level for each pattern
  vector log_lik(matrix p, matrix y) {
    int N = rows(y);
    int K = cols(y);
    vector[N] out;
    for(n in 1:N) {
      out[n] = 0.0;
      for(k in 1:K) {
        out[n] += y[n, k] * log(p[n, k]);
      }
    }
    return out;
  }
}

data {
  int N; // number of records
  int K; // number of response levels
  int P; // number of covariates
  matrix[N, K] y; // response record x level
  matrix[N, P] X; // design matrix
  vector[K] prior_counts; // prior for Dirichlet on cuts
  vector[P] prior_sd; // prior SD for beta coefficients
}

transformed data {
  vector[N] multinom_coef;
  for(i in 1:N)
    multinom_coef[i] = multinomial_coef(to_vector(y[i]));
}

parameters {
  simplex[K] pi;      // outcome level probabilities for reference level
  vector[P] beta_raw; // covariate coefficients
}

transformed parameters {
  vector[K-1] alpha;    // outcome level cuts for reference pattern
  matrix[N, K] p;   // matrix of level probabilities for covariate pattern and outcome level
  vector[N] loglik; // store covariate pattern loglik contribution
  vector[P] beta = prior_sd .* beta_raw;
  vector[N] eta = X * beta;
  alpha = make_cutpoints(pi, 1);
  p = Pr(alpha, eta);
  loglik = multinom_coef + log_lik(p, y);
}

model {
  // Prior model
  target += normal_lpdf(beta_raw | 0, 1);
  target += dirichlet_lpdf(pi | prior_counts);
  // Observational model
  target += loglik;
}
