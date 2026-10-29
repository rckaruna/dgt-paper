# study_A_lognormal.R -- Study A of preregistration.md (frozen 2026-09-25)
# Usage:  Rscript study_A_lognormal.R [--pilot] [--workers N] [--cells 1,3,9]
#   --pilot   : 3 replications per cell, for timing and sanity
#   --workers : parallel workers (default 4), one chain per core
# Output: results/study_A.csv (one row per replicate, appended; safe to
#         restart, completed (cell, rep) pairs are skipped) and
#         results/study_A_targets.csv.

source("common.R")
args <- commandArgs(trailingOnly = TRUE)
pilot <- "--pilot" %in% args
arg_value <- function(flag, default) {
  # accepts "--flag=VALUE" and "--flag VALUE"
  i <- grep(paste0("^", flag, "(=|$)"), args)
  if (length(i) == 0) return(default)
  v <- sub(paste0("^", flag, "=?"), "", args[i[1]])
  if (!nzchar(v) && length(args) > i[1]) v <- args[i[1] + 1]
  v
}
workers <- as.integer(arg_value("--workers", "4"))
cells_sel <- { v <- arg_value("--cells", ""); if (nzchar(v)) as.integer(strsplit(v, ",")[[1]]) else 1:9 }
R_reps <- if (pilot) 3L else 200L
dir.create("results", showWarnings = FALSE); dir.create("logs", showWarnings = FALSE)
out_csv <- if (pilot) "results/study_A_pilot.csv" else "results/study_A.csv"

# ---- design ---------------------------------------------------------------
design <- expand.grid(icc_eta = c(0.10, 0.25, 0.50), s2eta = c(0.25, 1.00, 2.00))
design$cell  <- seq_len(nrow(design))
design$s2p   <- design$icc_eta * design$s2eta
design$s2eps <- design$s2eta - design$s2p
design$icc_Y <- (exp(design$s2p) - 1) / (exp(design$s2eta) - 1)
design$O     <- design$icc_eta / design$icc_Y
utils::write.csv(design, "results/study_A_targets.csv", row.names = FALSE)
n_p <- 50L; n_j <- 30L; mu <- 6.5

gen_A <- function(cell, rep) {
  d <- design[cell, ]; set.seed(20260925L + 1000L * cell + rep)
  nu <- rnorm(n_p, 0, sqrt(d$s2p))
  dat <- expand.grid(person = factor(seq_len(n_p)), cond = seq_len(n_j))
  dat$y <- exp(mu + nu[dat$person] + rnorm(nrow(dat), 0, sqrt(d$s2eps)))
  dat
}

# ---- compile the model once ----------------------------------------------
message("Compiling template (Study A)...")
template <- brm(y ~ 1 + (1 | person), data = gen_A(1L, 1L), family = lognormal(),
                prior = c(prior(normal(6.5, 5), class = Intercept),
                          prior(student_t(3, 0, 1), class = sd),
                          prior(student_t(3, 0, 1), class = sigma)),
                chains = 1, iter = 200, refresh = 0, backend = "cmdstanr", silent = 2)

one_rep <- function(cell, rep) {
  t0 <- Sys.time()
  dat <- gen_A(cell, rep); seed <- 20260925L + 1000L * cell + rep
  f <- fit_with_refit(template, dat, seed)
  res <- dgt::dgt_icc(f$fit, person_group = "person")
  e <- post_summary(res$draws$icc_eta); y <- post_summary(res$draws$icc_Y)
  data.frame(cell = cell, rep = rep, seed = seed,
             icc_eta_med = e[["med"]], icc_eta_lo = e[["lo"]], icc_eta_hi = e[["hi"]],
             icc_Y_med = y[["med"]], icc_Y_lo = y[["lo"]], icc_Y_hi = y[["hi"]],
             O_med = stats::median(res$draws$overestimation),
             max_rhat = f$max_rhat, n_div = f$n_div, refit = f$refit, ok = f$ok,
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# ---- run ------------------------------------------------------------------
cl <- parallel::makeCluster(workers)
invisible(parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({library(brms); library(dgt)}); NULL }))
parallel::clusterExport(cl, c("template", "design", "gen_A", "n_p", "n_j", "mu",
                              "fit_with_refit", "fit_diagnostics", "post_summary", "one_rep"))
for (cell in cells_sel) {
  done <- done_keys(out_csv)
  todo <- setdiff(seq_len(R_reps), as.integer(sub("^\\S+ ", "", done[grepl(paste0("^", cell, " "), done)])))
  if (length(todo) == 0) { message("cell ", cell, " complete"); next }
  message(sprintf("cell %d: %d replicates to run", cell, length(todo)))
  for (chunk in split(todo, ceiling(seq_along(todo) / workers))) {
    rows <- parallel::parLapply(cl, chunk, function(r, cell) tryCatch(one_rep(cell, r),
                                error = function(e) data.frame(cell = cell, rep = r, error = conditionMessage(e))),
                                cell = cell)
    for (row in rows) {
      if ("error" %in% names(row)) { message("  rep ", row$rep, " ERROR: ", row$error); next }
      append_row(out_csv, row)
      message(sprintf("  cell %d rep %3d  ICC_eta %.3f  ICC_Y %.3f  rhat %.3f div %d  %.0fs",
                      cell, row$rep, row$icc_eta_med, row$icc_Y_med, row$max_rhat, row$n_div, row$secs))
    }
  }
}
parallel::stopCluster(cl)
message("Study A ", if (pilot) "pilot " else "", "finished: ", out_csv)
