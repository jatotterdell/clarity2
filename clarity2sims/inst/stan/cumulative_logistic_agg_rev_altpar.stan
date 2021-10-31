// James Totterdell
// Date: 2021-10-13
//
// parameterised in terms of cut-points rather than probabilities,
// but still using Dirichlet prior

functions {
  // Induce Dirichlet prior on outcome level probabilities
  // - c: The cut-points
  // - alpha: The prior weights to each level
  // - phi: optional centering variable
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
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
      out[n, 1] = log1m_exp(-log1p_exp(eta[n] + c[1])); // ln(1 - inv_logit(-eta[n] - c[1]))
      out[n, K] = -log1p_exp(eta[n] + c[K-1]);          // ln(inv_logit(-eta[n] - c[K-1]))
      for (k in 2:(K - 1)) {
        // ln(inv_logit(-eta[n] - c[k-1]) - inv_logit(-eta[n] - c[k]))
        out[n, k] = log_diff_exp(-log1p_exp(eta[n] + c[k-1]), -log1p_exp(eta[n] + c[k]));
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

parameters {
  ordered[K-1] alpha;      // outcome level probabilities for reference level
  vector[P] beta_raw; // covariate coefficients
}

transformed parameters {
  matrix[N, K] p;   // matrix of level probabilities for covariate pattern and outcome level
  vector[N] loglik; // store covariate pattern loglik contribution
  vector[P] beta = prior_sd .* beta_raw;
  vector[N] eta = X * beta;
  p = Pr(alpha, eta);
  loglik = log_lik(p, y);
}

model {
  // Prior model
  target += normal_lpdf(beta_raw | 0, 1);
  target += induced_dirichlet_lpdf(alpha | prior_counts, 0);
  // Observational model
  target += loglik;
}
