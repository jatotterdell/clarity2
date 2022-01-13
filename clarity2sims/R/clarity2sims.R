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
  draws <- posterior::as_draws_rvars(fit$draws(c("alpha", "beta", "eta")))
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
#' @param stage2 Which arms to take into Stage 2? "all" takes both controls, "first" always takes control 1
#' and "best" always takes the one with highest posterior mean log-odds.
#' We can equate "first" to a random choice if both controls are equivalent.
#' @param ... Other arguments to cmdstanr::sample, e.g. adapt_delta, chains, etc.
#' @return A list of trial related data.tables
#' @export
sim_clarity2_trial <- function(mod,
                               n_seq = seq(600, 2100, 300),
                               p_assign = rep(1 / 3, 3),
                               alpha = stats::qlogis(cumsum(c(16, 29, 32, 13, 2, 1, 1, 6) / 100)[1:7]),
                               eta = c(0, 0, 0),
                               eff_eps = 0.975,
                               stage2 = "all",
                               ...) {

  N <- length(p_assign)
  K <- length(alpha) + 1
  X <- mod[[2]]$X
  C <- mod[[2]]$C
  n_max <- max(n_seq)
  n_new <- diff(c(0, n_seq))
  n_int <- length(n_seq)
  y <- matrix(0, N, K)
  y_obs <- array(0,
                 dim = c(n_int, N, K),
                 dimnames = list("analysis" = seq_len(n_int),
                                 "variable" = seq_len(N),
                                 "level" = seq_len(K)))

  # Storage
  labs <- list("analysis" = 1:n_int, "arm" = 1:3)
  n_obs <- matrix(0, n_int, N, dimnames = labs)
  e_alpha <- matrix(0, n_int, K - 1, dimnames = list(labs[[1]], c("cut" = 1:(K - 1))))
  v_alpha <- e_alpha
  e_beta <- matrix(0, n_int, N - 1, dimnames = list(labs[[1]], labs[[2]][-1]))
  v_beta <- e_beta
  e_ctr <- matrix(0, n_int, 3, dimnames  = list("analysis" = seq_len(n_int), "contrast" = c("3v1", "3v2", "3v1and2")))
  pr_ctr <- hi_ctr <- lo_ctr <- v_ctr <- e_ctr

  for (i in 1:n_int) {
    # After stage 1, if not "all" then drop one of the control arms
    if (i == 2 & stage2 != "all" & N == 3) {
      if (stage2 == "first") {
        p_assign <- c(0.5, 0, 0.5)
      } else if (stage2 == "best" & e_beta[i - 1, 1] > 0) {
        p_assign <- c(0, 0.5, 0.5)
      } else {
        p_assign <- c(0.5, 0, 0.5)
      }
    }
    # Treatment assignment
    x_new <- permuted_block_rand(p_assign, n_new[i], 2 * N)[["trt"]]
    n_x <- table(factor(x_new, levels = seq_len(length(p_assign))))
    y_new <- generate_data_agg(n_x, alpha, eta)
    y <- y + y_new
    mod[[2]]$y <- y
    n_obs[i, ] <- rowSums(y)
    y_obs[i, , ] <- y

    fit <- update_model(mod[[1]], mod[[2]], ...)
    # Define contrasts
    ctr <- c(fit$eta[3] - fit$eta[1], fit$eta[3] - fit$eta[2], fit$eta[3] - posterior::rvar_max(c(fit$eta[1], fit$eta[2])))
    # Summarise contrasts
    e_ctr[i, ] <- posterior::E(ctr)
    v_ctr[i, ] <- posterior::var(ctr)
    tmp        <- posterior::quantile2(ctr, c(0.025, 0.975))
    lo_ctr[i, ]<- tmp[1, ]
    hi_ctr[i, ]<- tmp[2, ]
    pr_ctr[i, ]<- posterior::Pr(ctr > 0)
    # Other summaries
    e_alpha[i, ] <- posterior::E(fit$alpha)
    v_alpha[i, ] <- posterior::var(fit$alpha)
    e_beta[i, ]  <- posterior::E(fit$beta)
    v_beta[i, ]  <- posterior::var(fit$beta)

  }
  ind <- 1:i

  out_alpha <- list(
    e_alpha = e_alpha[ind, , drop = F],
    v_alpha = v_alpha[ind, , drop = F]
  )
  out_arm <- list(
    n_obs = n_obs[ind, , drop = F],
    e_beta = e_beta[ind, , drop = F],
    v_beta = v_beta[ind, , drop = F]
  )
  out_ctr <- list(
    e_ctr = e_ctr[ind, , drop = F],
    v_ctr = v_ctr[ind, , drop = F],
    lo_ctr = lo_ctr[ind, , drop = F],
    hi_ctr = hi_ctr[ind, , drop = F],
    pr_ctr = pr_ctr[ind, , drop = F]
  )
  y_obs <- y_obs[ind, , , drop = FALSE]

  return(list(
    alpha = list_as_dt(out_alpha),
    trial = list_as_dt(out_arm),
    contr = list_as_dt(out_ctr),
    yobs  = dcast(array_as_dt(y_obs), ... ~ paste0("y", level), value.var = "y_obs")
  ))
}
