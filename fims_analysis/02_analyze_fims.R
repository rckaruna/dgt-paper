# 02_analyze_fims.R -- Illustration 3: coefficients, D-study, SEM (dgt >= 0.5.1)
suppressPackageStartupMessages({ library(brms); library(dgt); library(posterior) })
stopifnot(packageVersion("dgt") >= "0.5.1")
fits <- list(aus_logit = readRDS("fims_aus_logit_fit.rds"),
             aus_probit = readRDS("fims_aus_probit_fit.rds"),
             jpn_logit = readRDS("fims_jpn_logit_fit.rds"))
long <- readRDS("fims_long.rds")
q <- function(x, d = 3) { v <- stats::quantile(x, c(.5, .025, .975), names = FALSE); list(median = v[1], lo = v[2], hi = v[3]) }

res <- list()
for (nm in names(fits)) {
  fit <- fits[[nm]]; link <- family(fit)$link
  message("== ", nm, " (", link, ")")
  set.seed(1)
  ip <- dgt_icc(fit, person_group = "person_id")     # persons as object, items as facet
  set.seed(1)
  ii <- dgt_icc(fit, person_group = "item_id")       # items as object, persons as facet
  dp <- ip$draws; di <- ii$draws
  # D-study on the response scale for a proportion correct over n items (Theorem 3),
  # and on the link scale, both from the person-as-object draws; required n for 0.70/0.80
  n_grid <- 1:60
  sb <- function(icc, n) outer(icc, n, function(r, k) k * r / (1 + (k - 1) * r))
  curve_Y   <- sb(dp$icc_Y_abs, n_grid); curve_eta <- sb(dp$icc_eta, n_grid)
  need <- function(curve, c) { m <- apply(curve, 2, median); i <- which(m >= c); if (length(i)) n_grid[i[1]] else NA }
  # individual-score SEM on the response scale for the observed proportion over 14 items:
  # Var(Ybar_14) = P(1-P)/14 + (13/14) Cov, Cov = ICC_Y * P(1-P); SEM = sqrt(Var(Ybar)(1 - Erho2(14)))
  P <- dp$mean_pi
  var_ybar <- P * (1 - P) / 14 + (13 / 14) * dp$icc_Y_abs * P * (1 - P)
  erho14 <- 14 * dp$icc_Y_abs / (1 + 13 * dp$icc_Y_abs)
  sem14 <- sqrt(var_ybar * (1 - erho14))
  # small-effect bound of Proposition (d) evaluated at the posterior means (probit only)
  bound <- NA
  if (link == "probit") {
    post <- as_draws_df(fit); so <- post[["sd_item_id__Intercept"]]; mu <- post[["b_Intercept"]]
    z0 <- mu / sqrt(1 + so^2); bound <- q(dnorm(z0)^2 / (pnorm(z0) * (1 - pnorm(z0))))
  }
  res[[nm]] <- list(
    link = link, n_persons = nlevels(droplevels(fit$data$person_id)), n_items = nlevels(droplevels(fit$data$item_id)),
    mean_correct = mean(fit$data$correct),
    persons = list(icc_eta = q(dp$icc_eta), icc_Y = q(dp$icc_Y_abs), icc_Y_rel = q(dp$icc_Y_rel), icc_I = q(dp$icc_I),
                   gap_eta_I = q(dp$icc_eta - dp$icc_I), ratio_Y_eta = q(dp$icc_Y_abs / dp$icc_eta),
                   n80_Y = need(curve_Y, 0.80), n80_eta = need(curve_eta, 0.80),
                   n70_Y = need(curve_Y, 0.70), n70_eta = need(curve_eta, 0.70),
                   erho2_14_Y = q(erho14), erho2_14_eta = q(14 * dp$icc_eta / (1 + 13 * dp$icc_eta)),
                   sem14 = q(sem14), prevalence = q(P), bound_d = bound),
    items = list(icc_eta = q(di$icc_eta), icc_Y = q(di$icc_Y_abs), icc_Y_rel = q(di$icc_Y_rel), icc_I = q(di$icc_I),
                 gap_eta_I = q(di$icc_eta - di$icc_I)),
    varcomp = list(sd_person = q(as_draws_df(fit)[["sd_person_id__Intercept"]]),
                   sd_item = q(as_draws_df(fit)[["sd_item_id__Intercept"]])),
    curves = list(n = n_grid, Y_med = apply(curve_Y, 2, median), Y_lo = apply(curve_Y, 2, quantile, .025),
                  Y_hi = apply(curve_Y, 2, quantile, .975), eta_med = apply(curve_eta, 2, median),
                  eta_lo = apply(curve_eta, 2, quantile, .025), eta_hi = apply(curve_eta, 2, quantile, .975)))
  cat(sprintf("  persons: ICC_eta %.3f  ICC_Y %.3f (rel %.3f)  ICC_I %.3f | n*(0.80) Y %s, eta %s | SEM14 %.3f\n",
              res[[nm]]$persons$icc_eta$median, res[[nm]]$persons$icc_Y$median, res[[nm]]$persons$icc_Y_rel$median,
              res[[nm]]$persons$icc_I$median, res[[nm]]$persons$n80_Y, res[[nm]]$persons$n80_eta, res[[nm]]$persons$sem14$median))
  cat(sprintf("  items:   ICC_eta %.3f  ICC_Y %.3f (rel %.3f)  ICC_I %.3f\n", res[[nm]]$items$icc_eta$median,
              res[[nm]]$items$icc_Y$median, res[[nm]]$items$icc_Y_rel$median, res[[nm]]$items$icc_I$median))
}
# link-invariance on the same data: Australia logit vs probit, persons as object
res$link_check <- list(
  icc_Y_logit = res$aus_logit$persons$icc_Y$median, icc_Y_probit = res$aus_probit$persons$icc_Y$median,
  icc_eta_logit = res$aus_logit$persons$icc_eta$median, icc_eta_probit = res$aus_probit$persons$icc_eta$median)
cat(sprintf("\nLink invariance (Australia, persons): ICC_Y logit %.3f vs probit %.3f; ICC_eta logit %.3f vs probit %.3f\n",
            res$link_check$icc_Y_logit, res$link_check$icc_Y_probit, res$link_check$icc_eta_logit, res$link_check$icc_eta_probit))
res$dgt_version <- as.character(packageVersion("dgt"))
saveRDS(res, "fims_results.rds")
cat("Saved fims_results.rds. Next: Rscript 03_make_outputs_fims.R\n")
