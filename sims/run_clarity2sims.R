#!/usr/bin/env Rscript

library(clarity2sims)
library(cmdstanr)
library(posterior)
library(data.table)
library(parallel)
library(optparse)

#  ----- Command line arguments -----
option_list <- list(
  make_option(c("-c", "--cores"), type = "integer", default = 15,
              help = "number of cores to use",
              metavar = "character"),
  make_option(c("-n", "--nsim"), type = "integer", default = 1000,
              help = "number of simulations to run under each configuration",
              metavar = "character")
);
opt <- parse_args(OptionParser(option_list = option_list))
num_cores <- opt$cores
num_sims  <- opt$nsim

#  ----- Compile the required model -----
ordmod <- clarity2sims:::compile_cmdstanr_mod()
ordmoddat <- list(
  N = 3, P = 2, K = 8,
  # X = rbind(0, c(1, 0), c(1, 1)),
  X = bayestestR::contr.orthonorm(3),
  prior_counts = rep(2 / 8, 8),
  prior_sd = rep(1, 2))
mod <- list(ordmod, ordmoddat)

#  ----- Randomisation -----
RNGkind("L'Ecuyer-CMRG")
set.seed(671359)
mc.reset.stream()

# ----- Specify configurations to explore -----
cfg <- CJ(
  sims = num_sims,
  n_seq = list(seq(600, 2100, 300)),
  eff_eps = 0.975,
  alpha = list(qlogis(cumsum(c(16, 28, 32, 12, 2, 2, 2, 6) / 100)[1:7])),
  eta = list(rep(0, 3),
             c(0, 0, log(1/1.1)),
             c(0, 0, log(1.1)),
             c(0, 0, log(1.2)),
             c(0, 0, log(1.3)),
             c(0, 0, log(1.5)),
             c(0, 0, log(2.0))),
  stage2 = c("all", "first", "best"),
  sorted = FALSE
)

# ----- Which configurations do we want to run? -----
run_row <- seq_len(nrow(cfg))

# ----- Loop over configurations and save results -----
for (z in run_row) {
  start_time <- Sys.time()

  res <- mclapply(1:cfg[z][["sims"]], function(j) {
    sim_clarity2_trial(
      mod,
      n_seq = unlist(cfg[z][["n_seq"]]),
      eff_eps = cfg[z][["eff_eps"]][[1]],
      alpha = cfg[z][["alpha"]][[1]],
      eta = cfg[z][["eta"]][[1]],
      stage2 = cfg[z][["stage2"]],
      refresh = 0)
  }, mc.cores = num_cores)

  resl_alpha <- rbindlist(lapply(res, \(x) x[["alpha"]]), idcol = "trial")
  resl_trial <- rbindlist(lapply(res, \(x) x[["trial"]]), idcol = "trial")
  resl_contr <- rbindlist(lapply(res, \(x) x[["contr"]]), idcol = "trial")
  resl_alpha[, analysis := as.numeric(analysis)]
  resl_trial[, analysis := as.numeric(analysis)]
  resl_contr[, analysis := as.numeric(analysis)]

  end_time <- Sys.time()

  saveRDS(list(cfg = cfg[z],
               alpha = resl_alpha,
               trial = resl_trial,
               contr = resl_contr,
               runtime = end_time - start_time),
          paste0("~/out_files/clarity2_sims/test_centered",
                 formatC(z, width = 2, flag = "0"), ".rds"))
}
