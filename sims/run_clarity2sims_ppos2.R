#!/usr/bin/env Rscript

suppressMessages({
  library(clarity2sims)
  library(rstan)
  library(posterior)
  library(data.table)
  library(parallel)
  library(matrixStats)
  library(mvtnorm)
  library(optparse)
  library(bayestestR)
})

# ----- Command line arguments -----
option_list <- list(
  make_option(c("-c", "--cores"),
              type = "integer", default = 10,
              help = "number of cores to use [default %default]", metavar = "number"
  ),
  make_option(c("-n", "--nsim"),
              type = "integer", default = 10,
              help = "number of simulations to run under each configuration [default %default]", metavar = "number"
  ),
  make_option(c("-f", "--filename"),
              type = "character", default = "ppos_sims_",
              help = "the output file name for the simulations [default %default]", metavar = "character"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
num_cores <- opt$cores
num_sims <- opt$nsim
file_name <- opt$filename


# ----- Compile the model -----
ordmod <- clarity2sims:::compile_rstan_mod()
ordmodat <- list(
  N = 3,
  K = 8,
  P = 2,
  y = matrix(0, 3, 8),
  X = bayestestR::contr.orthonorm(3),
  C = rbind(c(0, 1), c(-1, 1)),
  prior_counts = 2 * rep(1 / 8, 8),
  prior_sd = rep(1, 2)
)
approxmod <- clarity2sims:::compile_rstan_approx_mod()
mod <- list(ordmod, ordmodat, approxmod)

# ----- Specify configurations to explore -----
cfg <- CJ(
  sims = num_sims,
  n_seq = list(
    seq(600, 2100, 300),
    c(600, 1050, 1500, 1950, 2100),
    c(600, 1200, 1800, 2100)
  ),
  eff_eps = c(0.93, 0.94),
  fut_eps = 0.025,
  alpha = list(qlogis(cumsum(c(67, 24, 1, 1, 1, 1, 1, 4) / 100)[1:7])),
  eta = list(
    rep(0, 3),
    c(0, 0, log(1 / 1.1)),
    c(0, 0, log(1.1)),
    c(0, 0, log(1.2)),
    c(0, 0, log(1.3)),
    c(0, 0, log(1.5)),
    c(0, log(1.1), log(1.3))
  ),
  sorted = FALSE
)


# ----- Which configurations do we want to run? -----
run_row <- seq_len(nrow(cfg))

# ----- PRNGs -----
RNGkind("L'Ecuyer-CMRG")
seed <- 791351
set.seed(seed)
mc.reset.stream()

# ----- Loop over configurations and save results -----
for (z in run_row) {
  start_time <- Sys.time()

  res <- mclapply(1:cfg[z][["sims"]], function(j) {
    sim_clarity2_ppos_trial(
      mod,
      n_seq = cfg[z][["n_seq"]][[1]],
      eff_eps = cfg[z][["eff_eps"]][[1]],
      fut_eps = cfg[z][["fut_eps"]][[1]],
      alpha = cfg[z][["alpha"]][[1]],
      eta = cfg[z][["eta"]][[1]],
      stage2 = "all",
      drop_controls = TRUE,
      refresh = 0
    )
  }, mc.cores = num_cores)
  # invisible(runif(1)) # step Global RNG stream forward
  set.seed(seed) # use same seed for every configuration

  resl_alpha <- rbindlist(lapply(res, \(x) x[["alpha"]]), idcol = "trial")
  resl_contr <- rbindlist(lapply(res, \(x) x[["contr"]]), idcol = "trial")
  resl_trial <- rbindlist(lapply(res, \(x) x[["trial"]]), idcol = "trial")
  resl_yobs <- rbindlist(lapply(res, \(x) x[["yobs"]]), idcol = "trial")
  resl_alpha[, analysis := as.numeric(analysis)]
  resl_contr[, analysis := as.numeric(analysis)]
  resl_trial[, analysis := as.numeric(analysis)]
  resl_yobs[, analysis := as.numeric(analysis)]

  end_time <- Sys.time()

  saveRDS(
    list(
      cfg = cfg[z],
      alpha = resl_alpha,
      contr = resl_contr,
      trial = resl_trial,
      yobs = resl_yobs,
      runtime = end_time - start_time
    ),
    paste0(
      "~/out_files/clarity2_sims/", file_name, formatC(z, width = 2, flag = "0"), ".rds"
    )
  )
}
