#' @title update_model
#' @description Update the model used in clarity2sims
#' @param mod A compiled cmdstanr stanmodel
#' @param moddat Model data for use with mod
#' @param ... Additional arguments to sample
#' @return Posterior draws
#' @import posterior
update_model <- function(mod, moddat, ...) {
  drp <- utils::capture.output(fit <- mod$sample(
    data = moddat, ...
  ))
  draws <- posterior::as_draws_rvars(fit$draws(c("alpha", "beta")))
  return(draws)
}


#' @title sim_clarity2_trial
#' @description Simulate a clarity2 trial
#' @details
#' Simulates a Clarity 2.0 trial.
#' @param mod A list (model, model data)
#' @param n_seq Sequence of interim analysis sample sizes
#' @param p_assign Assignment probabilities to treatment arms
#' @param alpha True intercept parameter
#' @param eta True eta parameter
#' @param eff_eps Effectiveness threshold
#' @param ... Other arguments to cmdstanr::sample, e.g. adapt_delta, chains, etc.
#' @return A list of trial related data.tables
#' @export
sim_clarity2_trial <- function(mod,
                               n_seq = seq(600, 2100, 300),
                               p_assign = rep(1 / 3, 3),
                               alpha = stats::qlogis(cumsum(c(16, 29, 32, 13, 2, 1, 1, 6) / 100)[1:7]),
                               eta = c(0, -0.5, 0.5),
                               eff_eps = 0.975,
                               ...) {
  N <- length(p_assign)
  K <- length(alpha) + 1
  n_max <- max(n_seq)
  n_new <- diff(c(0, n_seq))
  n_int <- length(n_seq)
  y <- matrix(0, N, K)

  # Storage
  labs <- list(c("analysis" = 1:n_int), c("arm" = 1:3))
  n_obs <- matrix(0, n_int, N, dimnames = labs)
  e_alpha <- matrix(0, n_int, K - 1, dimnames = list(labs[[1]], c("cut" = 1:(K - 1))))
  v_alpha <- e_alpha
  e_beta <- matrix(0, n_int, N - 1, dimnames = list(labs[[1]], labs[[2]][-1]))
  v_beta <- matrix(0, n_int, N - 1, dimnames = list(labs[[1]], labs[[2]][-1]))
  pr_eff <- matrix(0, n_int, N - 1, dimnames = list(labs[[1]], labs[[2]][-1]))
  pr_eff2 <- matrix(0, n_int, 1, dimnames = list(labs[[1]], c("arm" = 3)))
  i_eff2 <- pr_eff2

  for (i in 1:n_int) {
    # Treatment assignment
    x_new <- permuted_block_rand(p_assign, n_new[i], 2 * N)[["trt"]]
    n_x <- table(factor(x_new, levels = seq_len(length(p_assign))))
    y_new <- generate_data_agg(n_x, alpha, eta)
    y <- y + y_new
    mod[[2]]$y <- y
    fit <- update_model(mod[[1]], mod[[2]], ...)
    n_obs[i, ] <- rowSums(y)
    e_alpha[i, ] <- posterior::E(fit$alpha)
    v_alpha[i, ] <- posterior::var(fit$alpha)
    e_beta[i, ] <- posterior::E(fit$beta)
    v_beta[i, ] <- posterior::var(fit$beta)
    # is arm 2,3 better than arm 1?
    pr_eff[i, ] <- posterior::Pr(fit$beta > 0)
    # is arm 3 better than arm 1 and 2?
    pr_eff2[i, ] <- posterior::Pr(fit$beta[2] > 0 & fit$beta[2] - fit$beta[1] > 0)
    i_eff2[i, ] <- as.integer(pr_eff2[i, ] > eff_eps)

    if (i_eff2[i, ] == 1) break
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
    pr_eff = pr_eff[ind, , drop = F],
    pr_eff2 = pr_eff2[ind, , drop = F],
    i_eff2 = i_eff2[ind, , drop = F]
  )

  return(list(
    alpha = list_as_dt(out_alpha),
    trial = list_as_dt(out_arm)
  ))
}
