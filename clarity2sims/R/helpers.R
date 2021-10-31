# %%%%% Miscellaneous functions %%%%%


#' @title var_to_string
#' @description Convert a variable name to a string
#' @param var The parameter
#' @return A string
var_to_string <- function(var) {
  deparse(substitute(var))
}

#' @title list_as_dt
#' @description Convert a trial list into a trial data.table
#' @param res A list which is the result of simulate_longitudinal_trial
#' @return The list converted to a data.table with one row being an analysis by treatment result
#' @import data.table
list_as_dt <- function(res) {
  tmp <- lapply(
    1:length(res), function(a) {
      # if(length(dim(res[[a]])) > 2) {
      #   return(as.data.table(res[[a]], value.name = names(res)[a]))
      # } else {
      melt(as.data.table(res[[a]], keep.rownames = "analysis"),
           id.vars = "analysis",
           value.name = names(res)[a])
      # }
    })
  tmp <- Reduce(function(x, y) merge(x, y, on = .(analysis, variable), all = TRUE), tmp)
  setkey(tmp, analysis, variable)
  return(tmp)
}


#' @title array_as_dt
#' @description Convert a trial array into a data.table
#' @param arr An array
#' @return The array converted to a data.table=
#' @import data.table
array_as_dt <- function(arr) {
  if(length(dim(arr)) > 2) {
    return(as.data.table(arr, value.name = deparse(substitute(arr))))
  } else {
    return(melt(
      data.table::as.data.table(arr, keep.rownames = "analysis"),
      id.vars = "analysis",
      value.name = deparse(substitute(arr))))
  }
}


#' @title log1m
#' @description Calculate log(1 - x)
#' See https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
#' @param x A real valued input
#' @return The result of log(1 - x)
log1m <- function(x) {
  if(x > 1) return(NA_real_)
  return(log1p(-x))
}


#' @title log1m_exp
#' @description Calculates log(1 - exp(x))
#' @param x A real valued input
#' @return The result of log(1 - exp(x))
log1m_exp <- function(x) {
  if (x > 0) {
    return(NA_real_);
  } else if (x > -0.6931472) {
    return(log(-expm1(x)))  # 0.693147 ~= log(2) = a0
  } else {
    return(log1m(exp(x)))
  }
}


#' @title log1p_exp
#' @description   Calculates log(1 + exp(x))
#' @param x A real valued input
#' @return The result of log(1 + exp(x))
log1p_exp <- function(x) {
  # like log_sum_exp with y=0.0; prevents underflow
  if (x > 0.0) {
    return(x + log1p(exp(-x)))
  }
  return(log1p(exp(x)))
}


#' @title log_sum_exp
#' @description Calculate log(exp(x) + exp(y))
#' @param x A real valued input
#' @param y A real valued input
#' @return The result of log(exp(x) + exp(y))
log_sum_exp <- function(x, y) {
  if(x == -Inf) return(y)
  if(x == Inf & y == Inf) return(Inf)
  if(x > y) return(x + log1p_exp(y - x))
  return(y + log1p_exp(x - y))
}


#' @title log_diff_exp
#' @description Calculate log(exp(x) - exp(y))
#' @param x A real valued input
#' @param y A real valued input
#' @return The result of log(exp(x) - exp(y))
log_diff_exp <- function(x, y) {
  if(x <= y) {
    if(x < Inf && x  == y) {
      return(-Inf)
    } else {
      return(NA_real_)
    }
  }
  return(x + log1m_exp(y - x))
}


#' @title ordinal_probs
#' @description Calculate marginal probabilities from cut-points
#' @param alpha A vector of cut-points (increasing)
#' @return A vector of level probabilities
ordinal_probs <- function(alpha) {
  K <- length(alpha) + 1
  out <- numeric(K)
  out[1] <- log1m_exp(-log1p_exp(alpha[1])) # ln(1 - inv_logit(-alpha[1]))
  out[K] <- -log1p_exp(alpha[K-1])          # ln(inv_logit(-alpha[K-1]))
  if(K == 2) return(exp(out))
  for(k in 2:(K-1)) {
    # ln(inv_logit(-alpha[k-1]) - inv_logit(-alpha[k]))
    out[k] = log_diff_exp(-log1p_exp(alpha[k-1]), -log1p_exp(alpha[k]))
  }
  return(exp(out))
}


#' @title ordinal_probs_mat
#' @description Calculate marginal probabilities matrix from cut-points
#' @param alpha A vector of cut-points (increasing)
#' @param eta A vector of linear predictors used to shift alpha (proportional odds)
#' @return A vector of level probabilities
ordinal_probs_mats <- function(alpha, eta) {
  K <- length(alpha) + 1
  N <- length(eta)
  out <- matrix(0, N, K)
  for(n in 1:N) {
    out[n, 1] <- log1m_exp(-log1p_exp(alpha[1] + eta[n])) # ln(1 - inv_logit(-eta[n]-alpha[1]))
    out[n, K] <- -log1p_exp(alpha[K-1] + eta[n])          # ln(inv_logit(-eta[n]-alpha[K-1]))
    if(K > 2) {
      for(k in 2:(K-1)) {
        # This is ln(inv_logit(-eta[n]-alpha[k-1]) - inv_logit(-eta[n]-alpha[k]))
        out[n, k] <- log_diff_exp(-log1p_exp(alpha[k - 1] + eta[n]), -log1p_exp(alpha[k] + eta[n]))
      }
    }
  }
  return(exp(out))
}


#' @title compile_rstan_mod
#' Compile provided cumulative logistic Stan model for rstan.
#'
#' @return A `stanmodel` object
#' @importFrom rstan stan_model
compile_rstan_mod <- function() {
  rstan::stan_model(file = system.file("stan", "cumulative_logistic_agg_rev.stan", package = "clarity2sims"))
}


#' @title compile_rstan_mod
#' Compile provided cumulative logistic Stan model using alternative parameterisation for rstan.
#'
#' @return A `stanmodel` object
#' @importFrom rstan stan_model
compile_rstan_approx_mod <- function() {
  rstan::stan_model(file = system.file("stan", "cumulative_logistic_agg_rev_altpar.stan", package = "clarity2sims"))
}


#' @title compile_cmdstanr_mod
#' Compile provided cumulative logistic Stan model for cmdstanr.
#'
#' @return A `CmdStanModel` object
#' @importFrom cmdstanr cmdstan_model
compile_cmdstanr_mod <- function() {
  cmdstanr::cmdstan_model(stan_file = system.file("stan", "cumulative_logistic_agg_rev.stan", package = "clarity2sims"))
}
