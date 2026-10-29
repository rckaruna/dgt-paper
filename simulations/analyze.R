# analyze.R -- summarise results/study_A.csv and results/study_B.csv into
# per-cell metrics (preregistration.md, "Metrics per cell"), apply the
# pre-registered success criteria, and write:
#   results/table_study_A.tex, results/table_study_B.tex, results/numbers.json
# Usage: Rscript analyze.R [--pilot]
args <- commandArgs(trailingOnly = TRUE); pilot <- "--pilot" %in% args
sfx <- if (pilot) "_pilot" else ""
suppressPackageStartupMessages(library(jsonlite))

mcse <- function(x) stats::sd(x) / sqrt(length(x))
cover <- function(lo, hi, truth) mean(lo <= truth & truth <= hi)
f3 <- function(x) sprintf("%.3f", round(x, 3) + 0)  # round first, "+ 0" turns -0 into 0, so a rounded zero prints as 0.000 rather than -0.000
f2 <- function(x) sprintf("%.2f", x)

numbers <- list()

# ---------------------------------------------------------------- Study A
if (file.exists(sprintf("results/study_A%s.csv", sfx))) {
  A  <- utils::read.csv(sprintf("results/study_A%s.csv", sfx))
  tA <- utils::read.csv("results/study_A_targets.csv")
  n_fail <- sum(!A$ok); A <- A[A$ok, ]
  sA <- do.call(rbind, lapply(split(A, A$cell), function(d) {
    t <- tA[tA$cell == d$cell[1], ]
    data.frame(cell = t$cell, icc_eta = t$icc_eta, s2eta = t$s2eta, icc_Y = t$icc_Y, O = t$O, R = nrow(d),
               bias_eta = mean(d$icc_eta_med) - t$icc_eta, mcse_eta = mcse(d$icc_eta_med),
               rmse_eta = sqrt(mean((d$icc_eta_med - t$icc_eta)^2)),
               cov_eta  = cover(d$icc_eta_lo, d$icc_eta_hi, t$icc_eta),
               bias_Y = mean(d$icc_Y_med) - t$icc_Y, mcse_Y = mcse(d$icc_Y_med),
               rmse_Y = sqrt(mean((d$icc_Y_med - t$icc_Y)^2)),
               cov_Y  = cover(d$icc_Y_lo, d$icc_Y_hi, t$icc_Y),
               width_Y = mean(d$icc_Y_hi - d$icc_Y_lo),
               O_mean = mean(d$O_med), rel_over = (t$icc_eta - t$icc_Y) / t$icc_Y,
               secs = mean(d$secs))
  }))
  sA$pass_bias <- abs(sA$bias_Y) <= 0.01; sA$pass_cov <- sA$cov_Y >= 0.93 & sA$cov_Y <= 0.97
  print(sA, digits = 3)
  cat(sprintf("Study A: %d fits excluded by the refit rule\n", n_fail))
  tex <- c("\\begin{table}[H]", "\\centering",
    "\\caption{Study A: Bayesian recovery of the lognormal one-facet ICCs ($n_p = 50$, $n_j = 30$, $R = 200$ per cell). Bias and RMSE of posterior medians against the analytic targets of Theorem~\\ref{thm:lognormal-icc}; coverage of 95\\% credible intervals; $\\mathcal{O}$ is the mean posterior-median overestimation ratio. Monte Carlo standard errors of the biases are at most " ,
    sprintf("%.4f.}", max(sA$mcse_Y, sA$mcse_eta)),
    "\\label{tab:sim-A}", "\\small", "\\begin{tabular}{cccccccccc}", "\\toprule",
    "$\\icc_\\eta$ & $\\sigma^2_\\eta$ & $\\icc_Y$ & Bias$_\\eta$ & RMSE$_\\eta$ & Cov$_\\eta$ & Bias$_Y$ & RMSE$_Y$ & Cov$_Y$ & $\\mathcal{O}$ \\\\", "\\midrule",
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", f2(sA$icc_eta), f2(sA$s2eta), f3(sA$icc_Y),
            f3(sA$bias_eta), f3(sA$rmse_eta), f2(sA$cov_eta), f3(sA$bias_Y), f3(sA$rmse_Y), f2(sA$cov_Y), f2(sA$O_mean)),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, sprintf("results/table_study_A%s.tex", sfx))
  numbers$study_A <- list(n_excluded = n_fail, max_abs_bias_Y = max(abs(sA$bias_Y)), min_cov_Y = min(sA$cov_Y),
                          max_cov_Y = max(sA$cov_Y), all_pass_bias = all(sA$pass_bias), all_pass_cov = all(sA$pass_cov),
                          O_range = range(sA$O_mean), rel_over_range = range(sA$rel_over), cells = sA)
}

# ---------------------------------------------------------------- Study B
if (file.exists(sprintf("results/study_B%s.csv", sfx))) {
  B  <- utils::read.csv(sprintf("results/study_B%s.csv", sfx), stringsAsFactors = FALSE)
  tB <- utils::read.csv("results/study_B_targets.csv", stringsAsFactors = FALSE)
  n_fail <- sum(!B$ok); B <- B[B$ok, ]
  sB <- do.call(rbind, lapply(split(B, B$cell), function(d) {
    t <- tB[tB$cell == d$cell[1], ]
    data.frame(cell = t$cell, link = t$link, icc_eta = t$icc_eta, P = t$P, R = nrow(d),
               t_abs = t$icc_Y_abs, t_rel = t$icc_Y_rel, t_I = t$icc_I, t_I_mcse = t$icc_I_mcse,
               bias_abs = mean(d$abs_med) - t$icc_Y_abs, mcse_abs = mcse(d$abs_med),
               rmse_abs = sqrt(mean((d$abs_med - t$icc_Y_abs)^2)), cov_abs = cover(d$abs_lo, d$abs_hi, t$icc_Y_abs),
               bias_rel = mean(d$rel_med) - t$icc_Y_rel, rmse_rel = sqrt(mean((d$rel_med - t$icc_Y_rel)^2)),
               cov_rel = cover(d$rel_lo, d$rel_hi, t$icc_Y_rel),
               bias_eta = mean(d$eta_med) - t$icc_eta, rmse_eta = sqrt(mean((d$eta_med - t$icc_eta)^2)),
               cov_eta = cover(d$eta_lo, d$eta_hi, t$icc_eta),
               bias_I = mean(d$I_med) - t$icc_I, rmse_I = sqrt(mean((d$I_med - t$icc_I)^2)),
               cov_I = cover(d$I_lo, d$I_hi, t$icc_I),
               kappa_mean = mean(d$kappa_sample), phi_mean = mean(d$phi_pairwise),
               kappa_gap = mean(d$kappa_sample) - t$icc_Y_abs, secs = mean(d$secs))
  }))
  sB$pass_bias <- abs(sB$bias_abs) <= 0.02; sB$pass_cov <- sB$cov_abs >= 0.92 & sB$cov_abs <= 0.98
  sB$pass_kappa <- abs(sB$kappa_gap) <= 0.02
  print(sB, digits = 3)
  cat(sprintf("Study B: %d fits excluded by the refit rule\n", n_fail))
  # link-invariance check: matched cells across links
  inv <- merge(sB[sB$link == "probit", c("icc_eta", "P", "t_abs", "bias_abs")],
               sB[sB$link == "logit",  c("icc_eta", "P", "t_abs", "bias_abs")], by = c("icc_eta", "P"), suffixes = c("_probit", "_logit"))
  tex <- c("\\begin{table}[H]", "\\centering",
    "\\caption{Study B: Bayesian recovery of the Bernoulli coefficients in a crossed persons $\\times$ items design ($n_p = 300$, $n_j = 20$, $\\sigma^2_o = 0.3$, $R = 100$ per cell). Targets from Proposition~\\ref{prop:binary-icc} by quadrature (probit absolute value checked against the exact bivariate-normal form). Bias and coverage of posterior medians and 95\\% intervals; $\\bar\\kappa$ is the mean sample Fleiss kappa, to be compared with the absolute $\\icc_Y$.}",
    "\\label{tab:sim-B}", "\\small", "\\begin{tabular}{llcccccccccc}", "\\toprule",
    "Link & $\\icc_\\eta$ & $P$ & $\\icc_Y$ & Bias$_Y$ & Cov$_Y$ & $\\icc_Y^{\\mathrm{rel}}$ & Bias$^{\\mathrm{rel}}$ & $\\icc_I$ & Bias$_I$ & Cov$_I$ & $\\bar\\kappa$ \\\\", "\\midrule",
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", sB$link, f2(sB$icc_eta), f2(sB$P), f3(sB$t_abs),
            f3(sB$bias_abs), f2(sB$cov_abs), f3(sB$t_rel), f3(sB$bias_rel), f3(sB$t_I), f3(sB$bias_I), f2(sB$cov_I), f3(sB$kappa_mean)),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, sprintf("results/table_study_B%s.tex", sfx))
  numbers$study_B <- list(n_excluded = n_fail, max_abs_bias_abs = max(abs(sB$bias_abs)), cov_abs_range = range(sB$cov_abs),
                          max_abs_kappa_gap = max(abs(sB$kappa_gap)), all_pass_bias = all(sB$pass_bias),
                          all_pass_cov = all(sB$pass_cov), all_pass_kappa = all(sB$pass_kappa),
                          link_invariance = inv, cells = sB)
}
write_json(numbers, sprintf("results/numbers%s.json", sfx), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("wrote results/numbers", sfx, ".json and the .tex tables\n", sep = "")
