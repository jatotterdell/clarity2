library(clarity2sims)
library(rstan)
library(posterior)
library(data.table)
library(parallel)
library(matrixStats)
library(mvtnorm)
library(optparse)

# ----- Command line arguments -----
option_list = list(
  make_option(c("-c", "--cores"), type="integer", default=15,
              help="number of cores to use", metavar="character"),
  make_option(c("-n", "--nsim"), type="integer", default=100,
              help="number of simulations to run under each configuration", metavar="character")
);
opt <- parse_args(OptionParser(option_list = option_list))
num_cores <- opt$cores
num_sims  <- opt$nsim

# ----- Compile the model -----
ordmod <- clarity2sims:::compile_rstan_mod()
ordmodat <- list(
  N = 3,
  K = 8,
  P = 2,
  y = matrix(0, 3, 8),
  X = rbind(0, c(1, 0), c(0, 1)),
  prior_counts = 2*rep(1/8, 8),
  prior_sd = rep(1, 2)
)
mod <- list(ordmod, ordmodat)


# ----- PRNGs -----
RNGkind("L'Ecuyer-CMRG")
set.seed(613570)

# ----- Specify configurations to explore -----
cfg <- CJ(
  sims = num_sims,
  n_seq = list(seq(600, 2100, 300)),
  eff_eps = 0.975,
  fut_eps = 0.025,
  eta = list(rep(0, 3),
             c(0, 0, log(1.1)),
             c(0, 0, log(1.2))),
  sorted = FALSE
)

# ----- Which configurations do we want to run? -----
run_row <- seq_len(nrow(cfg))

# ----- Loop over configurations and save results -----

for(z in run_row) {
  start_time <- Sys.time()

  res <- mclapply(1:cfg[z][["sims"]], function(j) {
    sim_clarity2_ppos_trial(
      mod,
      n_seq = cfg[z][["n_seq"]][[1]],
      eff_eps = cfg[z][["eff_eps"]][[1]],
      fut_eps = cfg[z][["fut_eps"]][[1]],
      eta = cfg[z][["eta"]][[1]],
      refresh = 0)
  }, mc.cores = num_cores)

  resl_alpha <- rbindlist(lapply(res, \(x) x[["alpha"]]), idcol = "trial")
  resl_contr <- rbindlist(lapply(res, \(x) x[["contr"]]), idcol = "trial")
  resl_trial <- rbindlist(lapply(res, \(x) x[["trial"]]), idcol = "trial")
  resl_yobs  <- rbindlist(lapply(res, \(x) x[["yobs"]]), idcol = "trial")
  resl_alpha[, analysis := as.numeric(analysis)]
  resl_contr[, analysis := as.numeric(analysis)]
  resl_trial[, analysis := as.numeric(analysis)]
  resl_yobs[, analysis := as.numeric(analysis)]

  end_time <- Sys.time()

  saveRDS(list(
    cfg = cfg[z],
    alpha = resl_alpha,
    contr = resl_contr,
    trial = resl_trial,
    yobs = resl_yobs,
    runtime = end_time - start_time),
    paste0("~/out_files/clarity2_sims/ppos_",
           formatC(z, width = 2, flag = "0"), ".rds"))
}
