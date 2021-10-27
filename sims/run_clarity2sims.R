library(clarity2sims)
library(cmdstanr)
library(posterior)
library(data.table)
library(parallel)
library(pbapply)
library(optparse)

#  ----- Command line arguments -----
option_list <- list(
  make_option(c("-c", "--cores"), type = "integer", default = 15,
              help = "number of cores to use",
              metavar = "character"),
  make_option(c("-n", "--nsim"), type = "integer", default = 10,
              help = "number of simulations to run under each configuration",
              metavar = "character")
);
opt <- parse_args(OptionParser(option_list = option_list))
num_cores <- opt$cores
num_sims  <- opt$nsim

#  ----- Compile the required model -----
ordmod <- clarity2sims:::compile_cmdstanr_mod()
ordmoddat <- list(N = 3, P = 2, K = 8, X = rbind(0, c(1, 0), c(0, 1)),
                  prior_counts = rep(2 / 8, 8), prior_sd = rep(1, 2))
mod <- list(ordmod, ordmoddat)

#  ----- Randomisation -----
RNGkind("L'Ecuyer-CMRG")
set.seed(671359)

# ----- Specify configurations to explore -----
cfg <- CJ(
  sims = num_sims,
  n_seq = list(600, 2100),
  eff_eps = 0.975,
  eta = list(rep(0, 3)),
  prior_counts = list(
    5 * c(16, 29, 32, 13, 2, 1, 1, 6) / 100,
    2 * rep(1 / 8, 8)),
  prior_sd = list(rep(1, 2), rep(0.25, 2)),
  sorted = FALSE
)

# ----- Which configurations do we want to run? -----
run_row <- seq_len(nrow(cfg))

# ----- Loop over configurations and save results -----
pboptions(type = "txt")
for (z in run_row) {
  start_time <- Sys.time()

  res <- pblapply(1:cfg[z][["sims"]], function(j) {
    sim_clarity2_trial(
      mod,
      n_seq = unlist(cfg[z][["n_seq"]]),
      eff_eps = cfg[z][["eff_eps"]][[1]],
      eta = cfg[z][["eta"]][[1]],
      refresh = 0)
  }, cl = num_cores)

  resl_alpha <- rbindlist(lapply(res, \(x) x[["alpha"]]), idcol = "trial")
  resl_trial <- rbindlist(lapply(res, \(x) x[["trial"]]), idcol = "trial")
  resl_alpha[, analysis := as.numeric(analysis)]
  resl_trial[, analysis := as.numeric(analysis)]

  end_time <- Sys.time()

  saveRDS(list(cfg = cfg[z],
               alpha = resl_alpha,
               trial = resl_trial,
               runtime = end_time - start_time),
          paste0("~/out_files/clarity2_sims/fixed_ss_",
                 formatC(z, width = 2, flag = "0"), ".rds"))
}
