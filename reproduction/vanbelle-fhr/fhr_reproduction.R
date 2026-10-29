## Reproduction of the fetal heart rate example of Vanbelle (2025), Section 6,
## with dgt 0.5.0.
##
## Data and model specification come from the author's own repository,
##   https://github.com/svanbelle/RELIABILITY.BINARY   (AGPL-3.0)
## cloned into ./RELIABILITY.BINARY (git-ignored; not vendored into this repo).
##
## Her published numbers (Vanbelle 2025, p. 14):
##   sigma_R^2  = .13 (.05, .33)
##   tau_abs    = .70 (.62, .77)     latent (probit link) scale
##   phi_abs    = .48 (.41, .55)     manifest (response) scale
##   kappa_abs  = .51 (.43, .58)     frequentist Fleiss / intraclass kappa
##
## Nothing in the dgt package or in any .tex file is modified by this script.

options(width = 120, digits = 6)

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dgt)
  library(posterior)
})

cat("=====================================================================\n")
cat("Vanbelle (2025) fetal heart rate reproduction\n")
cat("run at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("---------------------------------------------------------------------\n")
cat("R            :", R.version.string, "\n")
cat("dgt          :", as.character(packageVersion("dgt")), "\n")
cat("brms         :", as.character(packageVersion("brms")), "\n")
cat("cmdstanr     :", as.character(packageVersion("cmdstanr")), "\n")
cat("cmdstan      :", cmdstanr::cmdstan_version(), "\n")
cat("posterior    :", as.character(packageVersion("posterior")), "\n")
stopifnot(as.character(packageVersion("dgt")) == "0.5.0")
cat("CHECK: packageVersion(\"dgt\") == 0.5.0  -- OK\n")
cat("=====================================================================\n\n")

## ---------------------------------------------------------------------------
## 1. Data
## ---------------------------------------------------------------------------

repo_dir <- "RELIABILITY.BINARY"
if (!dir.exists(repo_dir)) {
  cat("Cloning svanbelle/RELIABILITY.BINARY ...\n")
  system2("git", c("clone", "https://github.com/svanbelle/RELIABILITY.BINARY",
                   repo_dir))
}
cat("source repository commit:\n")
system2("git", c("-C", repo_dir, "log", "-1",
                 shQuote("--format=%H | %ad | %s"), "--date=iso"))
cat("licence of that repository:",
    sub("^ *", "", readLines(file.path(repo_dir, "LICENSE"), n = 1)), "\n\n")

HR <- read.csv(file.path(repo_dir, "HR.csv"))
cat("HR.csv dimensions (recordings x midwives):", dim(HR), "\n")
cat("missing values in HR.csv:", sum(is.na(HR)), "\n")
cat("range of the raw FHR values (bpm):", range(HR), "\n")
stopifnot(nrow(HR) == 145L, ncol(HR) == 16L)
cat("CHECK: 145 recordings, 16 midwives -- OK\n\n")

## Her dichotomisation, EXAMPLES.R line 164:
##     HR_01 <- ifelse(HR_NA > 110 & HR_NA < 150, 0, 1)
## HR_NA is HR with unreadable entries set to NA; HR.csv as published has no
## missing entries, so HR_NA == HR here and the design is completely crossed.
HR_01 <- ifelse(HR > 110 & HR < 150, 0, 1)

cat("proportion classified abnormal, by midwife:\n")
print(round(colMeans(HR_01), 3))
cat("\noverall proportion abnormal P =", round(mean(HR_01), 4), "\n\n")

n    <- nrow(HR_01)
nrat <- ncol(HR_01)
dat  <- data.frame(
  y   = as.vector(HR_01),
  sub = factor(rep(seq_len(n), times = nrat)),
  met = factor(rep(seq_len(nrat), each = n))
)
cat("long data frame:", nrow(dat), "rows =", n, "x", nrat, "\n\n")

## ---------------------------------------------------------------------------
## 2. Sample coefficients from the raw data (no model)
## ---------------------------------------------------------------------------

## Fleiss kappa for m fixed raters, binary, complete design.
fleiss_kappa_binary <- function(M) {
  m  <- ncol(M)
  N  <- nrow(M)
  x1 <- rowSums(M)                 # number of "1" ratings per subject
  x0 <- m - x1
  P_i <- (x1 * (x1 - 1) + x0 * (x0 - 1)) / (m * (m - 1))
  Pbar <- mean(P_i)
  p1 <- sum(x1) / (N * m); p0 <- 1 - p1
  Pe <- p1^2 + p0^2
  (Pbar - Pe) / (1 - Pe)
}

## Mean of the pairwise phi (= Pearson) correlations over all rater pairs.
mean_pairwise_phi <- function(M) {
  m <- ncol(M)
  pr <- combn(m, 2)
  v <- apply(pr, 2, function(k) {
    a <- M[, k[1]]; b <- M[, k[2]]
    if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
    stats::cor(a, b)
  })
  list(mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
       range = range(v, na.rm = TRUE), n_pairs = ncol(pr))
}

k_fleiss <- fleiss_kappa_binary(HR_01)
phi_pair <- mean_pairwise_phi(HR_01)

cat("--- sample coefficients computed directly from the 145 x 16 table ---\n")
cat("Fleiss kappa (intraclass kappa, binary)      :", round(k_fleiss, 4), "\n")
cat("mean pairwise phi over", phi_pair$n_pairs, "rater pairs      :",
    round(phi_pair$mean, 4), "\n")
cat("   sd of the pairwise phi                    :", round(phi_pair$sd, 4), "\n")
cat("   range of the pairwise phi                 :",
    paste(round(phi_pair$range, 4), collapse = ", "), "\n")
cat("Vanbelle's published kappa_abs = .51 (.43, .58)\n")
cat("Vanbelle's published phi_abs   = .48 (.41, .55)\n\n")

## ---------------------------------------------------------------------------
## 3. The two-way probit Bernoulli model
## ---------------------------------------------------------------------------
## Her JAGS model probit6.txt (FUNCTIONS_EXAMPLES.R, lines 355-402):
##
##   y[i]   ~ dbern(p[i]);  p[i] <- phi(S[sub[i]] + R[met[i]])
##   S[j]   ~ dnorm(beta0, tau_S)      crossed subject intercepts
##   R[j]   ~ dnorm(0,     tau_R)      crossed midwife intercepts
##   beta0  ~ dnorm(0, 1)              precision 1, i.e. Normal(0, sd = 1)
##   sigma_S ~ dt(0, 1, 1) T(0,)       half-Cauchy(0, 1)
##   sigma_R ~ dt(0, 1, 1) T(0,)       half-Cauchy(0, 1)
##
## brms transcription: putting beta0 in the fixed part rather than in the mean
## of S[] is the same model.

prior_vb <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(cauchy(0, 1), class = "sd", group = "sub"),
  prior(cauchy(0, 1), class = "sd", group = "met")
)

cat("--- brms model and priors ---\n")
print(prior_vb)
cat("\n")

fit <- brm(
  y ~ 1 + (1 | sub) + (1 | met),
  data = dat,
  family = bernoulli(link = "probit"),
  prior = prior_vb,
  backend = "cmdstanr",
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4, seed = 20260911,
  control = list(adapt_delta = 0.95),
  refresh = 0,
  file = "fit_fhr_probit", file_refit = "on_change"
)

cat("\n--- model summary ---\n")
print(summary(fit))

cat("\n--- convergence ---\n")
cat("max Rhat :", round(max(rhat(fit), na.rm = TRUE), 4), "\n")
cat("min bulk ESS:", round(min(posterior::summarise_draws(
  posterior::as_draws_array(fit), "ess_bulk")$ess_bulk, na.rm = TRUE)), "\n")
cat("divergent transitions:",
    sum(nuts_params(fit)$Parameter == "divergent__" &
        nuts_params(fit)$Value > 0), "\n\n")

## ---------------------------------------------------------------------------
## 4. Variance components on her parameterisation
## ---------------------------------------------------------------------------

pd <- posterior::as_draws_df(fit)
sd_sub <- pd$sd_sub__Intercept
sd_met <- pd$sd_met__Intercept

q <- function(x) c(mean = mean(x), median = median(x),
                   lo = unname(quantile(x, 0.025)),
                   hi = unname(quantile(x, 0.975)))

cat("--- variance components (probit link, residual variance fixed at 1) ---\n")
cat("sigma_S  (sd, subjects) :", paste(round(q(sd_sub), 4), collapse = "  "), "\n")
cat("sigma_R  (sd, midwives) :", paste(round(q(sd_met), 4), collapse = "  "), "\n")
cat("sigma2_S (var, subjects):", paste(round(q(sd_sub^2), 4), collapse = "  "), "\n")
cat("sigma2_R (var, midwives):", paste(round(q(sd_met^2), 4), collapse = "  "), "\n")
cat("   [mean median lo hi]\n")
cat("Vanbelle's published sigma_R^2 = .13 (.05, .33)\n\n")

## Her Equation (11), the absolute latent-scale reliability.
tau_abs_manual <- sd_sub^2 / (sd_sub^2 + sd_met^2 + 1)
cat("her Eq. (11) recomputed by hand from our posterior:\n")
cat("tau_abs =", paste(round(q(tau_abs_manual), 4), collapse = "  "), "\n")
cat("Vanbelle's published tau_abs = .70 (.62, .77)\n\n")

## ---------------------------------------------------------------------------
## 5. dgt 0.5.0
## ---------------------------------------------------------------------------

cat("--- dgt::dgt_icc(fit, person_group = \"sub\") ---\n")
icc <- dgt::dgt_icc(fit, person_group = "sub", seed = 20260911)
print(icc)
cat("\nsummary table:\n")
print(icc$summary)

cat("\n--- per-draw ICC columns available ---\n")
print(names(icc$draws))

cat("\n--- posterior summaries of the dgt draws ---\n")
for (nm in c("icc_eta", "icc_Y_abs", "icc_Y_rel", "icc_I")) {
  cat(sprintf("%-10s : %s\n", nm, paste(round(q(icc$draws[[nm]]), 4),
                                        collapse = "  ")))
}
cat("   [mean median lo hi]\n\n")

cat("--- agreement between dgt icc_eta and her Eq. (11) ---\n")
cat("max |icc_eta - tau_abs_manual| =",
    max(abs(icc$draws$icc_eta - tau_abs_manual)), "\n\n")

## ---------------------------------------------------------------------------
## 6. Side-by-side table
## ---------------------------------------------------------------------------

fmt <- function(x) sprintf("%.3f (%.3f, %.3f)",
                           median(x), quantile(x, 0.025), quantile(x, 0.975))

comparison <- data.frame(
  quantity = c("sigma_R^2", "ICC_eta (latent)", "ICC_Y absolute (manifest)",
               "ICC_Y relative (manifest)", "Fleiss kappa (sample)",
               "mean pairwise phi (sample)"),
  ours = c(fmt(sd_met^2), fmt(icc$draws$icc_eta), fmt(icc$draws$icc_Y_abs),
           fmt(icc$draws$icc_Y_rel), sprintf("%.3f", k_fleiss),
           sprintf("%.3f", phi_pair$mean)),
  vanbelle = c(".13 (.05, .33)", "tau_abs = .70 (.62, .77)",
               "phi_abs = .48 (.41, .55)", "-- (not reported)",
               "kappa_abs = .51 (.43, .58)", "phi_abs = .48 (.41, .55)"),
  stringsAsFactors = FALSE
)
cat("--- side by side ---\n")
print(comparison, right = FALSE, row.names = FALSE)

saveRDS(list(icc = icc, comparison = comparison, fleiss = k_fleiss,
             phi_pair = phi_pair, q_sd_sub = q(sd_sub), q_sd_met = q(sd_met)),
        "results_fhr.rds")

cat("\ndone.\n")
cat("dgt version in this log:", as.character(packageVersion("dgt")), "\n")
