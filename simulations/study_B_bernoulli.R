# study_B_bernoulli.R -- Study B of preregistration.md (frozen 2026-09-25)
# Usage:  Rscript study_B_bernoulli.R [--pilot] [--workers N] [--cells 1,5]
# Output: results/study_B.csv (one row per replicate, appended, resume-safe)
#         results/study_B_targets.csv (population targets per cell)

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
cells_sel <- { v <- arg_value("--cells", ""); if (nzchar(v)) as.integer(strsplit(v, ",")[[1]]) else 1:12 }
R_reps <- if (pilot) 3L else 100L
dir.create("results", showWarnings = FALSE); dir.create("logs", showWarnings = FALSE)
out_csv <- if (pilot) "results/study_B_pilot.csv" else "results/study_B.csv"

# ---- design and population targets ---------------------------------------
design <- expand.grid(link = c("probit", "logit"), icc_eta = c(0.10, 0.25, 0.50), P = c(0.50, 0.20),
                      stringsAsFactors = FALSE)
design$cell <- seq_len(nrow(design))
targets_csv <- "results/study_B_targets.csv"
if (file.exists(targets_csv)) {
  targets <- utils::read.csv(targets_csv, stringsAsFactors = FALSE)
} else {
  message("Computing population targets (quadrature, exact probit check, ICC_I reference x5)...")
  targets <- do.call(rbind, lapply(seq_len(nrow(design)), function(c) {
    d <- design[c, ]; t <- bernoulli_targets(d$link, d$icc_eta, d$P, s2o = 0.3, seed = 100 + c)
    data.frame(cell = c, t, stringsAsFactors = FALSE)
  }))
  utils::write.csv(targets, targets_csv, row.names = FALSE)
}
print(targets[, c("cell", "link", "icc_eta", "P", "s2p", "mu", "icc_Y_abs", "icc_Y_rel", "icc_I", "icc_I_mcse")])
n_p <- 300L; n_j <- 20L

gen_B <- function(cell, rep) {
  t <- targets[cell, ]; set.seed(20260926L + 1000L * cell + rep)
  invlink <- dgt:::.inv_link_fun(t$link)
  nu <- rnorm(n_p, 0, sqrt(t$s2p)); om <- rnorm(n_j, 0, sqrt(t$s2o))
  dat <- expand.grid(person = factor(seq_len(n_p)), item = factor(seq_len(n_j)))
  dat$y <- rbinom(nrow(dat), 1, invlink(t$mu + nu[dat$person] + om[dat$item]))
  dat
}

# ---- one template per link ----------------------------------------------
message("Compiling templates (Study B)...")
templates <- list()
for (lk in c("probit", "logit")) {
  c0 <- which(design$link == lk)[1]
  templates[[lk]] <- brm(y ~ 1 + (1 | person) + (1 | item), data = gen_B(c0, 1L),
                         family = bernoulli(link = lk),
                         prior = c(prior(normal(0, 2.5), class = Intercept),
                                   prior(student_t(3, 0, 1), class = sd)),
                         chains = 1, iter = 200, refresh = 0, backend = "cmdstanr", silent = 2)
}

one_rep <- function(cell, rep) {
  t0 <- Sys.time()
  t <- targets[cell, ]; dat <- gen_B(cell, rep); seed <- 20260926L + 1000L * cell + rep
  f <- fit_with_refit(templates[[t$link]], dat, seed, chains = 4, iter = 1000, warmup = 500)
  set.seed(seed)
  res <- dgt::dgt_icc(f$fit, person_group = "person")
  d <- res$draws
  a <- post_summary(d$icc_Y_abs); r <- post_summary(d$icc_Y_rel)
  e <- post_summary(d$icc_eta);   i <- post_summary(d$icc_I)
  Y <- matrix(dat$y, nrow = n_p, ncol = n_j)          # expand.grid order: person varies fastest
  data.frame(cell = cell, rep = rep, seed = seed, link = t$link,
             abs_med = a[["med"]], abs_lo = a[["lo"]], abs_hi = a[["hi"]],
             rel_med = r[["med"]], rel_lo = r[["lo"]], rel_hi = r[["hi"]],
             eta_med = e[["med"]], eta_lo = e[["lo"]], eta_hi = e[["hi"]],
             I_med = i[["med"]],   I_lo = i[["lo"]],   I_hi = i[["hi"]],
             kappa_sample = fleiss_kappa_binary(Y), phi_pairwise = mean_pairwise_phi(Y),
             prev_sample = mean(Y),
             max_rhat = f$max_rhat, n_div = f$n_div, refit = f$refit, ok = f$ok,
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# ---- run ------------------------------------------------------------------
cl <- parallel::makeCluster(workers)
invisible(parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({library(brms); library(dgt)}); NULL }))
parallel::clusterExport(cl, c("templates", "targets", "design", "gen_B", "n_p", "n_j",
                              "fit_with_refit", "fit_diagnostics", "post_summary", "one_rep",
                              "fleiss_kappa_binary", "mean_pairwise_phi"))
for (cell in cells_sel) {
  done <- done_keys(out_csv)
  todo <- setdiff(seq_len(R_reps), as.integer(sub("^\\S+ ", "", done[grepl(paste0("^", cell, " "), done)])))
  if (length(todo) == 0) { message("cell ", cell, " complete"); next }
  message(sprintf("cell %d (%s, ICC_eta %.2f, P %.2f): %d replicates to run",
                  cell, targets$link[cell], targets$icc_eta[cell], targets$P[cell], length(todo)))
  for (chunk in split(todo, ceiling(seq_along(todo) / workers))) {
    rows <- parallel::parLapply(cl, chunk, function(r, cell) tryCatch(one_rep(cell, r),
                                error = function(e) data.frame(cell = cell, rep = r, error = conditionMessage(e))),
                                cell = cell)
    for (row in rows) {
      if ("error" %in% names(row)) { message("  rep ", row$rep, " ERROR: ", row$error); next }
      append_row(out_csv, row)
      message(sprintf("  cell %2d rep %3d  abs %.3f rel %.3f eta %.3f I %.3f  kappa %.3f  rhat %.3f div %d  %.0fs",
                      cell, row$rep, row$abs_med, row$rel_med, row$eta_med, row$I_med,
                      row$kappa_sample, row$max_rhat, row$n_div, row$secs))
    }
  }
}
parallel::stopCluster(cl)
message("Study B ", if (pilot) "pilot " else "", "finished: ", out_csv)
