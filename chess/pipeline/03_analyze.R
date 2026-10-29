## 03_analyze.R
## ---------------------------------------------------------------------
## Load cached brms fits and compute DGT quantities for Section 8.2.
##
## For response times (lognormal): uses dgt::dgt_info_icc (analytic
## branch via transformation invariance).
##
## For accuracy (Bernoulli): dgt does not yet support this family for
## information ICC, so we compute I(nu_p; Y) by Monte Carlo in-script.
##
## Usage:  Rscript 03_analyze.R
## ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(brms)
  library(dgt)
  library(dplyr)
  library(tibble)
})

set.seed(20260418)

summarise_posterior <- function(x, digits = 3) {
  q <- quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  list(
    median = unname(q[2]),
    lower  = unname(q[1]),
    upper  = unname(q[3]),
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    pretty = sprintf(sprintf("%%.%df [%%.%df, %%.%df]", digits, digits, digits),
                     q[2], q[1], q[3])
  )
}

## =====================================================================
## 1. Response-time analysis (lognormal model)
## =====================================================================

chess_rt_fit <- readRDS("chess_rt_fit.rds")
post_rt      <- as_draws_df(chess_rt_fit)

sp2   <- post_rt$sd_person_id__Intercept^2
sw2   <- post_rt$sd_item_id__Intercept^2
se2   <- post_rt$sigma^2
seta2 <- sp2 + sw2 + se2

## --- ICCs with PERSON as object --------------------------------------
## dgt_info_icc with person_group = "person_id" uses the person SD as
## s2p and treats (item + residual) variance as s2e.

info_rt_person <- dgt::dgt_info_icc(
  fit          = chess_rt_fit,
  person_group = "person_id"
)

icc_eta_p <- info_rt_person$draws$icc_eta
icc_I_p   <- info_rt_person$draws$icc_I

icc_Y_p <- (exp(sp2) - 1) / (exp(seta2) - 1)
O_p     <- icc_eta_p / icc_Y_p

## --- ICCs with ITEM as object ----------------------------------------
## Pass person_group = "item_id".  dgt_info_icc treats it generically
## as "the grouping variable whose variance we want as signal."

info_rt_item <- dgt::dgt_info_icc(
  fit          = chess_rt_fit,
  person_group = "item_id"
)

icc_eta_i <- info_rt_item$draws$icc_eta
icc_I_i   <- info_rt_item$draws$icc_I
icc_Y_i   <- (exp(sw2) - 1) / (exp(seta2) - 1)
O_i       <- icc_eta_i / icc_Y_i

## --- D-study: required occasions for target E-rho^2 ------------------

d_study_n <- function(icc_draws, target) {
  n <- ceiling(target * (1 - icc_draws) / (icc_draws * (1 - target)))
  list(
    target        = target,
    n_star_median = median(n, na.rm = TRUE),
    n_star_lower  = quantile(n, 0.025, na.rm = TRUE),
    n_star_upper  = quantile(n, 0.975, na.rm = TRUE)
  )
}

dstudy_p_80 <- d_study_n(icc_Y_p, 0.80)
dstudy_p_70 <- d_study_n(icc_Y_p, 0.70)
dstudy_i_80 <- d_study_n(icc_Y_i, 0.80)
dstudy_i_70 <- d_study_n(icc_Y_i, 0.70)

## --- D-study curves --------------------------------------------------

dstudy_curve <- function(icc_draws, n_grid) {
  out <- matrix(NA_real_, nrow = length(n_grid), ncol = 3)
  for (k in seq_along(n_grid)) {
    n  <- n_grid[k]
    er <- n * icc_draws / (1 + (n - 1) * icc_draws)
    out[k, ] <- quantile(er, c(0.025, 0.5, 0.975), na.rm = TRUE)
  }
  tibble(n = n_grid, lower = out[, 1], median = out[, 2], upper = out[, 3])
}

n_grid      <- 1:80
curve_p_eta <- dstudy_curve(icc_eta_p, n_grid)
curve_p_y   <- dstudy_curve(icc_Y_p,   n_grid)
curve_i_eta <- dstudy_curve(icc_eta_i, n_grid)
curve_i_y   <- dstudy_curve(icc_Y_i,   n_grid)

## =====================================================================
## 2. Accuracy analysis (Bernoulli model, Theorem 5.5)
## =====================================================================
##
## For each posterior draw s of (alpha, tau_obj, tau_facet):
##   1. Draw K_obj object-level effects u from N(0, tau_obj^2).
##   2. For each u_k, inner Monte Carlo over v ~ N(0, tau_facet^2)
##      gives p(Y = 1 | u_k) = E_v[plogis(alpha + u_k + v)].
##   3. p(Y = 1) marginal = mean over k of p(Y = 1 | u_k).
##   4. H(Y) from p(Y = 1) marginal; H(Y | u) from mean H over u_k.
##   5. I = H(Y) - H(Y | u); icc_I = 1 - exp(-2 I).

chess_acc_fit <- readRDS("chess_acc_fit.rds")
post_acc      <- as_draws_df(chess_acc_fit)

alpha  <- post_acc$b_Intercept
tau_p  <- post_acc$sd_person_id__Intercept
tau_w  <- post_acc$sd_item_id__Intercept
tau_p2 <- tau_p^2
tau_w2 <- tau_w^2

H_bern <- function(p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  -p * log(p) - (1 - p) * log(1 - p)
}

info_icc_bernoulli <- function(alpha, sd_object, sd_facet,
                               K_obj = 500, K_facet = 500, seed = 20260418) {
  set.seed(seed)
  S      <- length(alpha)
  icc_I  <- numeric(S)
  I_vals <- numeric(S)

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_object[s]
    sf <- sd_facet[s]

    u <- rnorm(K_obj, 0, so)

    p_cond <- numeric(K_obj)
    for (k in seq_len(K_obj)) {
      v          <- rnorm(K_facet, 0, sf)
      p_cond[k]  <- mean(plogis(a + u[k] + v))
    }

    p_marg <- mean(p_cond)

    H_Y         <- H_bern(p_marg)
    H_Y_given_u <- mean(H_bern(p_cond))

    I         <- max(H_Y - H_Y_given_u, 0)
    I_vals[s] <- I
    icc_I[s]  <- 1 - exp(-2 * I)
  }

  list(I = I_vals, icc_I = icc_I)
}

logit_resid_var <- pi^2 / 3
tau_tot_var     <- tau_p2 + tau_w2 + logit_resid_var

icc_eta_acc_p <- tau_p2 / tau_tot_var
icc_eta_acc_i <- tau_w2 / tau_tot_var

message("Computing Bernoulli information ICC (person as object)...")
info_acc_p <- info_icc_bernoulli(alpha, tau_p, tau_w, K_obj = 500, K_facet = 500)

message("Computing Bernoulli information ICC (item as object)...")
info_acc_i <- info_icc_bernoulli(alpha, tau_w, tau_p, K_obj = 500, K_facet = 500)

info_loss_p <- icc_eta_acc_p - info_acc_p$icc_I
info_loss_i <- icc_eta_acc_i - info_acc_i$icc_I

## ---------------------------------------------------------------------
## 2b. Response-scale ICC_Y for accuracy (revision, dgt >= 0.5.0)
## ---------------------------------------------------------------------
## Definition 1 on the 0/1 scale: Cov(Y_ij, Y_ij') / Var(Y_ij) with the
## facet integrated out (absolute coefficient). For a single Bernoulli
## trial this is the population intraclass kappa (Fleiss) and the phi
## coefficient under compound symmetry. dgt_icc() also returns the
## relative coefficient (facet main effect removed from the error term)
## and reads the link from the fit; the chess fit uses the logit link.

stopifnot(packageVersion("dgt") >= "0.5.0")
stopifnot(identical(family(chess_acc_fit)$link, "logit"))

message("Computing Bernoulli response-scale ICC_Y (person as object)...")
set.seed(20260911)
icc_acc_p <- dgt::dgt_icc(chess_acc_fit, person_group = "person_id")

message("Computing Bernoulli response-scale ICC_Y (item as object)...")
set.seed(20260911)
icc_acc_i <- dgt::dgt_icc(chess_acc_fit, person_group = "item_id")

icc_Y_acc_p     <- icc_acc_p$draws$icc_Y_abs
icc_Y_rel_acc_p <- icc_acc_p$draws$icc_Y_rel
icc_Y_acc_i     <- icc_acc_i$draws$icc_Y_abs
icc_Y_rel_acc_i <- icc_acc_i$draws$icc_Y_rel

## Consistency guard: dgt's own ICC_eta must agree with the hand
## computation above (same posterior draws, same pi^2/3 convention).
stopifnot(max(abs(icc_acc_p$draws$icc_eta - icc_eta_acc_p)) < 1e-8)
stopifnot(max(abs(icc_acc_i$draws$icc_eta - icc_eta_acc_i)) < 1e-8)

## ---------------------------------------------------------------------
## 2c. Accuracy D-study and individual-score SEM (revision, AE-7)
## ---------------------------------------------------------------------
## Theorem 3 on the 0/1 scale: the reliability of the proportion correct
## over n items is the Spearman-Brown transform of the absolute ICC_Y
## (= population intraclass kappa); the link-scale curve uses ICC_eta.
## SEM of an individual's proportion correct over the 40 ACT items:
##   Var(Ybar_40) = P(1-P)/40 + (39/40) ICC_Y P(1-P),  P = marginal P(Y=1)
##   SEM_Y(40)    = sqrt(Var(Ybar_40) (1 - Erho2_Y(40)))
## (same construction as fims_analysis/02_analyze_fims.R for 14 items).
n_items_acc <- 40L
sb <- function(icc, n) n * icc / (1 + (n - 1) * icc)
P_acc <- icc_acc_p$draws$mean_pi
stopifnot(length(P_acc) == length(icc_Y_acc_p))
erho2_40_Y   <- sb(icc_Y_acc_p, n_items_acc)
erho2_40_eta <- sb(icc_eta_acc_p, n_items_acc)
var_ybar_40  <- P_acc * (1 - P_acc) / n_items_acc +
                ((n_items_acc - 1) / n_items_acc) * icc_Y_acc_p * P_acc * (1 - P_acc)
sem_40_Y     <- sqrt(var_ybar_40 * (1 - erho2_40_Y))

dstudy_acc_p_80_Y   <- d_study_n(icc_Y_acc_p,   0.80)
dstudy_acc_p_70_Y   <- d_study_n(icc_Y_acc_p,   0.70)
dstudy_acc_p_80_eta <- d_study_n(icc_eta_acc_p, 0.80)
dstudy_acc_p_70_eta <- d_study_n(icc_eta_acc_p, 0.70)
dstudy_acc_i_80_Y   <- d_study_n(icc_Y_acc_i,   0.80)
dstudy_acc_i_80_eta <- d_study_n(icc_eta_acc_i, 0.80)

## =====================================================================
## 3. Load data summaries
## =====================================================================

rt_data  <- readRDS("chess_rt_data.rds")
acc_data <- readRDS("chess_acc_data.rds")

## =====================================================================
## 4. Assemble results
## =====================================================================

results <- list(
  rt = list(
    variance_components = list(
      sigma_p2   = summarise_posterior(sp2),
      sigma_w2   = summarise_posterior(sw2),
      sigma_eps2 = summarise_posterior(se2),
      sigma_eta2 = summarise_posterior(seta2)
    ),
    icc_person = list(
      icc_eta = summarise_posterior(icc_eta_p),
      icc_Y   = summarise_posterior(icc_Y_p),
      icc_I   = summarise_posterior(icc_I_p),
      O       = summarise_posterior(O_p, digits = 2)
    ),
    icc_item = list(
      icc_eta = summarise_posterior(icc_eta_i),
      icc_Y   = summarise_posterior(icc_Y_i),
      icc_I   = summarise_posterior(icc_I_i),
      O       = summarise_posterior(O_i, digits = 2)
    ),
    dstudy = list(
      person_80 = dstudy_p_80,
      person_70 = dstudy_p_70,
      item_80   = dstudy_i_80,
      item_70   = dstudy_i_70
    ),
    dstudy_curves = list(
      person_classical = curve_p_eta,
      person_dgt       = curve_p_y,
      item_classical   = curve_i_eta,
      item_dgt         = curve_i_y
    )
  ),
  acc = list(
    variance_components = list(
      tau_p2 = summarise_posterior(tau_p2),
      tau_w2 = summarise_posterior(tau_w2)
    ),
    icc_person = list(
      icc_eta   = summarise_posterior(icc_eta_acc_p),
      icc_Y     = summarise_posterior(icc_Y_acc_p),
      icc_Y_rel = summarise_posterior(icc_Y_rel_acc_p),
      icc_I     = summarise_posterior(info_acc_p$icc_I),
      I         = summarise_posterior(info_acc_p$I, digits = 3),
      info_loss = summarise_posterior(info_loss_p)
    ),
    icc_item = list(
      icc_eta   = summarise_posterior(icc_eta_acc_i),
      icc_Y     = summarise_posterior(icc_Y_acc_i),
      icc_Y_rel = summarise_posterior(icc_Y_rel_acc_i),
      icc_I     = summarise_posterior(info_acc_i$icc_I),
      I         = summarise_posterior(info_acc_i$I, digits = 3),
      info_loss = summarise_posterior(info_loss_i)
    ),
    dstudy = list(
      person_80_Y   = dstudy_acc_p_80_Y,
      person_70_Y   = dstudy_acc_p_70_Y,
      person_80_eta = dstudy_acc_p_80_eta,
      person_70_eta = dstudy_acc_p_70_eta,
      item_80_Y     = dstudy_acc_i_80_Y,
      item_80_eta   = dstudy_acc_i_80_eta,
      n_items       = n_items_acc,
      erho2_40_Y    = summarise_posterior(erho2_40_Y),
      erho2_40_eta  = summarise_posterior(erho2_40_eta),
      sem_40_Y      = summarise_posterior(sem_40_Y),
      prevalence    = summarise_posterior(P_acc)
    ),
    dgt_version = as.character(packageVersion("dgt")),
    link        = family(chess_acc_fit)$link
  ),
  data_summary = list(
    n_persons_rt  = nlevels(rt_data$person_id),
    n_items_rt    = nlevels(rt_data$item_id),
    n_obs_rt      = nrow(rt_data),
    n_persons_acc = nlevels(acc_data$person_id),
    n_items_acc   = nlevels(acc_data$item_id),
    n_obs_acc     = nrow(acc_data)
  )
)

saveRDS(results, "chess_results.rds")

## =====================================================================
## 5. Print headline numbers
## =====================================================================

cat("\n============================================================\n")
cat("Section 8.2 --- Headline results\n")
cat("============================================================\n\n")

cat("RESPONSE TIMES (Lognormal; Theorems 3.1, 5.3, 6.1(a))\n")
cat("-----------------------------------------------------\n")
cat("  Persons as object:\n")
cat(sprintf("    ICC_eta = %s\n", results$rt$icc_person$icc_eta$pretty))
cat(sprintf("    ICC_Y   = %s\n", results$rt$icc_person$icc_Y$pretty))
cat(sprintf("    ICC_I   = %s  (should equal ICC_eta by Cor 5.4)\n",
            results$rt$icc_person$icc_I$pretty))
cat(sprintf("    O       = %s\n", results$rt$icc_person$O$pretty))
cat(sprintf("    n* (E-rho^2 >= 0.80) = %d\n", dstudy_p_80$n_star_median))
cat(sprintf("    n* (E-rho^2 >= 0.70) = %d\n", dstudy_p_70$n_star_median))

cat("  Items as object:\n")
cat(sprintf("    ICC_eta = %s\n", results$rt$icc_item$icc_eta$pretty))
cat(sprintf("    ICC_Y   = %s\n", results$rt$icc_item$icc_Y$pretty))
cat(sprintf("    ICC_I   = %s  (should equal ICC_eta by Cor 5.4)\n",
            results$rt$icc_item$icc_I$pretty))
cat(sprintf("    O       = %s\n", results$rt$icc_item$O$pretty))
cat(sprintf("    n* (E-rho^2 >= 0.80) = %d\n", dstudy_i_80$n_star_median))
cat(sprintf("    n* (E-rho^2 >= 0.70) = %d\n", dstudy_i_70$n_star_median))

cat("\nACCURACY (Bernoulli; Theorems 5.5, 6.1(b))\n")
cat("------------------------------------------\n")
cat("  Persons as object:\n")
cat(sprintf("    ICC_eta       = %s\n", results$acc$icc_person$icc_eta$pretty))
cat(sprintf("    ICC_I         = %s\n", results$acc$icc_person$icc_I$pretty))
cat(sprintf("    Information I = %s nats\n", results$acc$icc_person$I$pretty))
cat(sprintf("    Info loss     = %s  (positive by Thm 5.5)\n",
            results$acc$icc_person$info_loss$pretty))

cat("  Items as object:\n")
cat(sprintf("    ICC_eta       = %s\n", results$acc$icc_item$icc_eta$pretty))
cat(sprintf("    ICC_I         = %s\n", results$acc$icc_item$icc_I$pretty))
cat(sprintf("    Information I = %s nats\n", results$acc$icc_item$I$pretty))
cat(sprintf("    Info loss     = %s  (positive by Thm 5.5)\n",
            results$acc$icc_item$info_loss$pretty))
cat("  Accuracy D-study (persons as object, proportion correct):\n")
cat(sprintf("    ICC_Y (abs)   = %s;  ICC_Y (rel) = %s\n",
            results$acc$icc_person$icc_Y$pretty, results$acc$icc_person$icc_Y_rel$pretty))
cat(sprintf("    n*(0.80): Y %g, eta %g;  n*(0.70): Y %g, eta %g\n",
            results$acc$dstudy$person_80_Y$n_star_median, results$acc$dstudy$person_80_eta$n_star_median,
            results$acc$dstudy$person_70_Y$n_star_median, results$acc$dstudy$person_70_eta$n_star_median))
cat(sprintf("    E-rho^2(40): Y %s, eta %s;  SEM_Y(40) = %s;  P(correct) = %s\n",
            results$acc$dstudy$erho2_40_Y$pretty, results$acc$dstudy$erho2_40_eta$pretty,
            results$acc$dstudy$sem_40_Y$pretty, results$acc$dstudy$prevalence$pretty))
cat(sprintf("  Items as object: n*(0.80): Y %g, eta %g\n",
            results$acc$dstudy$item_80_Y$n_star_median, results$acc$dstudy$item_80_eta$n_star_median))

message("\nDone. Next: Rscript 04_make_outputs.R")
