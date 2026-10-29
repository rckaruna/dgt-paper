## 02_fit_acc.R
## ---------------------------------------------------------------------
## Fit Bernoulli crossed-effects model to Amsterdam Chess Test (ACT)
## item correctness.  Overnight batch script; caches to chess_acc_fit.rds.
##
## Usage:  Rscript 02_fit_acc.R
## ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(LNIRT)
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
})

set.seed(20260418)

## ---- Load and preprocess ---------------------------------------------

data(AmsterdamChess)
act <- AmsterdamChess

y_cols  <- grep("^Y([0-9]+)$",  names(act), value = TRUE)
stopifnot(length(y_cols) == 40L)

y_wide <- act[, y_cols, drop = FALSE]
y_wide[y_wide == 9] <- NA                     # LNIRT missing code for Y

all_missing <- apply(y_wide, 1, function(r) all(is.na(r)))
message(sprintf("Dropping %d fully-missing persons (rows %s).",
                sum(all_missing),
                paste(which(all_missing), collapse = ", ")))
y_keep <- !all_missing

y_long <- y_wide[y_keep, , drop = FALSE] |>
  as.data.frame() |>
  mutate(person_id = factor(seq_len(sum(y_keep)))) |>
  pivot_longer(
    cols      = all_of(y_cols),
    names_to  = "item_id",
    values_to = "correct"
  ) |>
  filter(!is.na(correct)) |>
  mutate(
    item_id = factor(item_id, levels = y_cols),
    correct = as.integer(correct)
  )

stopifnot(all(y_long$correct %in% c(0L, 1L)))

message(sprintf(
  "Accuracy dataset: N=%d observations; %d persons; %d items; overall p(correct) = %.3f",
  nrow(y_long), nlevels(y_long$person_id), nlevels(y_long$item_id),
  mean(y_long$correct)))

saveRDS(y_long, "chess_acc_data.rds")

## ---- Priors ----------------------------------------------------------

acc_priors <- c(
  prior(normal(0, 2),       class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sd")
)

## ---- Fit -------------------------------------------------------------

chess_acc_fit <- brm(
  formula = correct ~ 1 + (1 | person_id) + (1 | item_id),
  data    = y_long,
  family  = bernoulli(link = "logit"),
  prior   = acc_priors,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 1000,
  backend = "cmdstanr",
  seed    = 20260418,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file    = "chess_acc_fit"
)

## ---- Diagnostics -----------------------------------------------------

cat("\n---- brms summary ----\n")
print(summary(chess_acc_fit), digits = 3)

cat("\n---- Posterior variance components (logit scale) ----\n")
post <- as_draws_df(chess_acc_fit)
tau_p2 <- post$sd_person_id__Intercept^2
tau_w2 <- post$sd_item_id__Intercept^2

qp <- function(x) sprintf("%.3f [%.3f, %.3f]",
                          median(x), quantile(x, 0.025), quantile(x, 0.975))
cat(sprintf("tau_p^2 (person, logit): %s\n", qp(tau_p2)))
cat(sprintf("tau_w^2 (item, logit)  : %s\n", qp(tau_w2)))

message("\nDone. Next: Rscript 03_analyze.R")
