#' @title impute_data
#' @description Impute future data using posterior predictive
#' @param npred A vector of length N giving the number of future
#' participants assigned to treatment n
#' @param post_p An array of posterior draws of dimension N x K x B,
#' where N is treatment arms, K is outcome levels, and B is repetitions
#' @return An array of same dimension as post_p, but with posterior-predicted
#' aggregated data.
impute_data <- function(npred, post_p) {
  N <- dim(post_p)[1]
  K <- dim(post_p)[2]
  B <- dim(post_p)[3]
  yagg <- array(0, dim = list(N, K, B))
  for (n in 1:N) {
    yagg[n, , ] <- apply(post_p[n, , ], 2, function(j) stats::rmultinom(1, npred[n], prob = j))
  }
  return (yagg)
}


#' @title approximate_posterior
#' @description Use Laplace approximation to posterior
#' @param mod Stan model
#' @param moddat Model data
#' @param collapse_levels Optimisation likely to fail if any outcome levels have zero observations.
#' If this parameter is TRUE (default) then drops outcome levels with zero observations and shares
#' the prior weight equally across all remaining levels
#' @return A list with posterior mean and covariance
approximate_posterior <- function(mod, moddat, collapse_levels = TRUE) {
  level_totals <- colSums(moddat$y)
  if(collapse_levels & any(level_totals == 0)) {
    moddat$y <- moddat$y[, level_totals != 0]
    moddat$K <- sum(level_totals != 0)
    moddat$prior_counts <- moddat$prior_counts[level_totals != 0] * sum(moddat$prior_counts) / sum(moddat$prior_counts[level_totals != 0])
  }
  fit <- rstan::optimizing(mod, data = moddat, hessian = TRUE, draws = 0)
  conv <- fit$return_code # Should check for convergence before summarising
  id <- grepl("beta_raw", colnames(fit$hessian))
  M <- fit$par[grepl("beta\\[", names(fit$par))]
  # If there's issues with covariance matrix just return NA
  V <- tryCatch ({
    H <- fit$hessian
    invH <- solve(-H, tol = sqrt(.Machine$double.eps))[id, id]
    v <- diag(moddat$prior_sd) %*% invH  %*% diag(moddat$prior_sd)
    if (!isSymmetric(v)) {
      crossprod(chol(v))
    } else {
      v
    }
  }, error = function(e) e)
  if (inherits(V, "error")) V <- NA
  return (list(M = M, V = V))
}


#' @title calc_ppos
#' @description Calculate posterior probability of success.
#' @details Calculates PPoS using Laplace approximation
#' @param mod An rstan model
#' @param dat The current data for the model
#' @param yimp The imputed sets (impute_data) of data to be combined with dat
#' @param epsilon Threshold for deciding superiority at maximum sample size
#' @param n_draws Number of posterior draws to use in calculating PPoS
#' @importFrom rstan optimizing
calc_ppos <- function(mod, dat, yimp, epsilon = 0.975, n_draws = 0) {
  B <- dim(yimp)[3]
  yobs <- dat$y
  pr_clt0 <- matrix(0, B, 3)
  Cmat <- rbind(c(0, 1), c(-1, 1))  # Or whatever contrast
  for (b in 1:B) {
    dat$y <- yobs + yimp[, , b]
    fit <- approximate_posterior(mod, dat)
    M <- fit$M
    V <- fit$V
    if (any(is.na(V))) {
      pr_clt0[b, ] <- NA
    } else {
      CM <- drop(Cmat %*% M)
      CVCt <- Cmat %*% V %*% t(Cmat)
      pr_clt0[b, -3] <- 1 - stats::pnorm(0, CM, sqrt(diag(CVCt)))
      pr_clt0[b, 3] <- mvtnorm::pmvnorm(lower = 0, mean = CM, sigma = CVCt)
    }
  }
  return (pr_clt0)
}


#' @title assess_futility
#' @description Assess futility rule at given sample size.
#' @details The observation model is calculated using dynamic HMC.
#' Posterior predictive draws are generated from this model.
#' The PPoS calculation uses Laplace approximations for the imputed data posteriors.
#' @param mod The assumed model
#' @param moddat The data required for the model (excluding outcomes which will be simulated)
#' @param n The observed sample size in each arm
#' @param p The probability of each outcome level in each arm
#' @param n_max The maximum sample size in each arm
#' @param n_sim The number of simulations to run
#' @param B The number of imputations for PPoS
#' @return A list with two elements:
#' - the first element is an array of dimension B x nrow(p) - 1 x n_sim which each
#'   element a predicted posterior probability that beta is less than zero
#' - the second element is a matrix of posterior probabilities
assess_futility <- function(mod, moddat, n, p, n_max, n_sim = 1, B = 500) {
  res <- parallel::mclapply(seq_len(n_sim), function(i) {
    yobs <- generate_data_p(n, p)
    moddat$y <- yobs
    fitobs <- sampling(mod, data = moddat, chains = 4, iter = 3000, warmup = 500, refresh = 0)
    post_draws <- as_draws_matrix(fitobs)
    post_p <- post_draws[, grepl("p\\[", colnames(post_draws))]
    P <- aperm(array(c(post_p), dim = list(nrow(post_p), nrow(p), ncol(p))), c(2, 3, 1))
    yimp <- impute_data(n_max - n, P[, , sample.int(dim(P)[3], size = B)])
    ppos <- calc_ppos(mod, moddat, yimp)
    pprob <- colMeans(post_draws[, grepl("beta\\[", colnames(post_draws))] < 0)
    ebeta <- colMeans(post_draws[, grepl("beta\\[", colnames(post_draws))])
    return (list(ppos = ppos, pprob = pprob, ebeta = ebeta))
  }, mc.cores = 15)
  return (list(
    ppos = simplify2array(lapply(res, function(x) x$ppos)),
    pprob = sapply(res, function(x) x$pprob),
    ebeta = sapply(res, function(x) x$ebeta)
  ))
}



#' @title sim_clarity2_ppos_trial
#' @description Simulate Clarity 2 with futility check
#' @details Note, assume the model is reverse coded, so OR > 1 => Pr(y < k) increases.
#' @param mod A list including the stan model as first element, and standata as second, and approximation stan model as third element
#' @param n_seq Sequence of interim analysis sample sizes
#' @param p_assign Assignment probabilities
#' @param alpha Intercept parameter
#' @param eta Effect parameter
#' @param eff_eps Effectiveness threshold
#' @param fut_eps Futility threshold
#' @param B_ppos Number of draws to use for calculating PPoS
#' @param ... Additional arguments to rstan::sampling
#' @return A list
#' @importFrom rstan sampling
#' @export
sim_clarity2_ppos_trial <- function(
  mod,
  n_seq = seq(600, 2100, 300),
  p_assign = rep(1 / 3, 3),
  alpha = stats::qlogis(cumsum(c(16, 28, 32, 12, 2, 2, 2, 6) / 100)[1:7]),
  eta = c(0, -0.5, 0.5),
  eff_eps = 0.975,
  fut_eps = 0.02,
  B_ppos = 500,
  ...) {

  N <- length(p_assign)
  K <- length(alpha) + 1
  X <- rbind(0, diag(1, N - 1))
  n_max <- max(n_seq)
  n_new <- diff(c(0, n_seq))
  n_int <- length(n_seq)
  n_left <- n_max - cumsum(n_new)
  p <- ordinal_probs_mats(alpha, eta)

  # Storage
  labs <- list("analysis" = 1:n_int, "variable" = 1:3)
  n_obs <- matrix(0, n_int, N, dimnames = labs)
  e_alpha <- matrix(0, n_int, K - 1, dimnames = list("analysis" = 1:n_int, "variable" = 1:(K - 1)))
  v_alpha <- e_alpha
  e_beta <- matrix(0, n_int, N - 1, dimnames = list("analysis" = 1:n_int, "variable" = 2:3))
  v_beta <- matrix(0, n_int, N - 1, dimnames = list("analysis" = 1:n_int, "variable" = 2:3))
  pr_eff <- matrix(0, n_int, N - 1, dimnames = list("analysis" = 1:n_int, "variable" = 2:3))
  pr_ctr <- matrix(0, n_int, 3, dimnames = list("analysis" = 1:n_int, "contrast" = c("3v1", "3v2", "3v1and2")))
  ppos <- matrix(0, n_int, 3, dimnames = list("analysis" = 1:n_int, "contrast" = c("3v1", "3v2", "3v1and2")))
  i_ctr <- pr_ctr
  i_eff <- pr_eff
  i_ppos <- ppos
  y <- matrix(0, N, K)
  y_obs <- array(0,
                 dim = c(n_int, N, K),
                 dimnames = list("analysis" = seq_len(n_int),
                                 "variable" = seq_len(N),
                                 "level" = seq_len(K)))

  for (i in 1:n_int) {
    # Treatment assignment
    x_new <- permuted_block_rand(p_assign, n_new[i], 2 * N)[["trt"]]
    n_x <- table(factor(x_new, levels = seq_len(length(p_assign))))
    y_new <- generate_data_p(n_x, p)
    y <- y + y_new
    y_obs[i, , ] <- y
    mod[[2]]$y <- y
    fit <- rstan::sampling(mod[[1]], data = mod[[2]], refresh = 0, chains = 4, iter = 3000, warmup = 500, ...)
    post_draws <- as_draws_matrix(fit)
    post_p <- post_draws[, grepl("p\\[", colnames(post_draws))]
    post_alpha <- post_draws[, grepl("alpha", colnames(post_draws))]
    post_beta <- post_draws[, grepl("beta\\[", colnames(post_draws))]
    if (i < n_int) {
      P <- aperm(array(c(post_p), dim = list(nrow(post_p), nrow(p), ncol(p))), c(2, 3, 1))
      x_imp <- table(factor(permuted_block_rand(p_assign, n_left[i], 2 * N)[["trt"]], levels = seq_len(N)))
      y_imp <- impute_data(x_imp, P[, , sample.int(dim(P)[3], size = B_ppos)])
      ppos[i, ] <- matrixStats::colMeans2(calc_ppos(mod[[3]], mod[[2]], y_imp, n_draws = 0) > eff_eps)
    } else {
      ppos[i, ] <- NA
    }
    n_obs[i, ] <- rowSums(y)
    e_alpha[i, ] <- matrixStats::colMeans2(post_alpha)
    v_alpha[i, ] <- diag(var(post_alpha))
    e_beta[i, ] <- matrixStats::colMeans2(post_beta)
    v_beta[i, ] <- diag(var(post_beta))
    # is arm 2,3 better than arm 1?
    pr_eff[i, ] <- matrixStats::colMeans2(post_beta > 0)
    # is arm 3 better than arm 2, and arm 1 and 2?
    pr_ctr[i, ] <- c(pr_eff[i, 2],
                     matrixStats::colMeans2(post_beta[, 2] - post_beta[, 1] > 0),
                     matrixStats::colMeans2(post_beta[, 2] > 0 & post_beta[, 2] - post_beta[, 1] > 0))
    i_eff[i, ] <- as.integer(pr_eff[i, ] > eff_eps)
    i_ctr[i, ] <- as.integer(pr_ctr[i, ] > eff_eps)
    i_ppos[i, ] <- as.integer(ppos[i, ] < fut_eps)
    # if (i < n_int && i_ppos[i, 1] == 1) break
  }
  ind <- 1:i

  # Summarise outputs by relevant groupings
  out_alpha <- list(
    e_alpha = e_alpha[ind, , drop = F],
    v_alpha = v_alpha[ind, , drop = F]
  )

  out_arm <- list(
    n_obs = n_obs[ind, , drop = F],
    e_beta = e_beta[ind, , drop = F],
    v_beta = v_beta[ind, , drop = F],
    pr_eff = pr_eff[ind, , drop = F]
  )

  out_ctr <- list(
    pr_ctr = pr_ctr[ind, , drop = F],
    ppos   = ppos[ind, , drop = F]
  )

  y_obs <- y_obs[ind, , , drop = FALSE]

  return (list(
    alpha = list_as_dt(out_alpha),
    trial = list_as_dt(out_arm),
    contr = list_as_dt(out_ctr),
    yobs = dcast(array_as_dt(y_obs), ... ~ paste0("y", level), value.var = "y_obs")
  ))
}
