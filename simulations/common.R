# common.R -- helpers shared by study_A_lognormal.R and study_B_bernoulli.R
# Pre-registered in preregistration.md (frozen 2026-09-25). Do not change
# designs here; record any deviation as an amendment in that file.

suppressPackageStartupMessages({
  library(brms); library(posterior); library(dgt)
})
stopifnot(packageVersion("dgt") >= "0.5.1")

# ---- resume-safe logging --------------------------------------------------
done_keys <- function(csv) {
  if (!file.exists(csv)) return(character(0))
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(character(0))
  paste(d$cell, d$rep)
}
append_row <- function(csv, row) {
  utils::write.table(row, csv, sep = ",", row.names = FALSE,
                     col.names = !file.exists(csv), append = file.exists(csv))
}

# ---- fit diagnostics and the pre-registered refit rule -------------------
fit_diagnostics <- function(fit) {
  rh <- suppressWarnings(brms::rhat(fit))
  nd <- tryCatch(sum(brms::nuts_params(fit, pars = "divergent__")$Value),
                 error = function(e) NA_real_)
  list(max_rhat = max(rh, na.rm = TRUE), n_div = nd)
}
fit_with_refit <- function(template, newdata, seed, adapt_delta = 0.9,
                           chains = 4, iter = 2000, warmup = 1000) {
  fit <- update(template, newdata = newdata, seed = seed, refresh = 0,
                chains = chains, iter = iter, warmup = warmup, cores = 1,
                control = list(adapt_delta = adapt_delta), silent = 2)
  dg <- fit_diagnostics(fit); refit <- 0L
  if (is.na(dg$n_div) || dg$n_div > 0 || dg$max_rhat > 1.05) {
    fit <- update(template, newdata = newdata, seed = seed + 1L, refresh = 0,
                  chains = chains, iter = iter, warmup = warmup, cores = 1,
                  control = list(adapt_delta = 0.99), silent = 2)
    dg <- fit_diagnostics(fit); refit <- 1L
  }
  ok <- !(is.na(dg$n_div) || dg$n_div > 0 || dg$max_rhat > 1.05)
  list(fit = fit, max_rhat = dg$max_rhat, n_div = dg$n_div, refit = refit, ok = ok)
}
post_summary <- function(x) {
  q <- stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  c(med = q[2], lo = q[1], hi = q[3])
}

# ---- sample binary agreement statistics ----------------------------------
# Y: persons x items 0/1 matrix, no missing values.
fleiss_kappa_binary <- function(Y) {
  n <- ncol(Y); x <- rowSums(Y); pbar <- mean(Y)
  1 - sum(x * (n - x)) / (nrow(Y) * n * (n - 1) * pbar * (1 - pbar))
}
mean_pairwise_phi <- function(Y) {
  m <- suppressWarnings(stats::cor(Y)); mean(m[upper.tri(m)], na.rm = TRUE)
}

# ---- Bernoulli population targets ----------------------------------------
gh <- function(n = 80) { x <- dgt:::.gh_nodes(n); x }   # dgt's own nodes/weights
solve_mu <- function(P, s2p, s2o, link) {
  if (link == "probit") return(stats::qnorm(P) * sqrt(1 + s2p + s2o))
  invlink <- dgt:::.inv_link_fun(link); g <- gh(80)
  marg <- function(mu) {
    eta <- outer(mu + sqrt(s2p) * g$x, sqrt(s2o) * g$x, "+")
    as.numeric(t(g$w) %*% invlink(eta) %*% g$w)
  }
  stats::uniroot(function(m) marg(m) - P, c(-15, 15), tol = 1e-10)$root
}
bernoulli_targets <- function(link, icc_eta, P, s2o = 0.3, seed = 1) {
  rg  <- dgt:::.link_resid_var(link)
  s2p <- icc_eta * (s2o + rg) / (1 - icc_eta)
  mu  <- solve_mu(P, s2p, s2o, link)
  invlink <- dgt:::.inv_link_fun(link)
  cc <- dgt:::.binomial_icc_core(mu, sqrt(s2p), sqrt(s2o), 1, n_nodes = 80, invlink = invlink)
  exact_abs <- NA_real_
  if (link == "probit") {
    V <- 1 + s2p + s2o; Pm <- stats::pnorm(mu / sqrt(V))
    P11 <- mvtnorm::pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), mean = c(mu, mu),
                            sigma = matrix(c(V, s2p, s2p, V), 2, 2))
    exact_abs <- as.numeric((P11 - Pm^2) / (Pm * (1 - Pm)))
    stopifnot(abs(exact_abs - cc[["icc_abs"]]) < 1e-4)
  }
  set.seed(seed)
  iccI <- replicate(5, dgt:::.binomial_info_facet(mu, sqrt(s2p), sqrt(s2o), 1,
                                                    K = 20000, M = 2000, invlink = invlink))
  list(link = link, icc_eta = icc_eta, P = P, s2p = s2p, s2o = s2o, mu = mu,
       icc_Y_abs = unname(cc[["icc_abs"]]), icc_Y_rel = unname(cc[["icc_rel"]]),
       icc_Y_abs_exact = exact_abs, icc_I = mean(iccI), icc_I_mcse = stats::sd(iccI),
       mean_pi = unname(cc[["mean_pi"]]))
}
