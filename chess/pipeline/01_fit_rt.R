## 01_fit_rt.R
## ---------------------------------------------------------------------
## Fit lognormal crossed-effects model to Amsterdam Chess Test (ACT)
## response times.  Overnight batch script; caches to chess_rt_fit.rds.
##
## Usage:  Rscript 01_fit_rt.R
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

rt_cols <- grep("^RT([0-9]+)$", names(act), value = TRUE)
stopifnot(length(rt_cols) == 40L)

rt_wide <- act[, rt_cols, drop = FALSE]
rt_wide[rt_wide == 10000.000] <- NA           # LNIRT missing code

all_missing <- apply(rt_wide, 1, function(r) all(is.na(r)))
message(sprintf("Dropping %d fully-missing persons (rows %s).",
                sum(all_missing),
                paste(which(all_missing), collapse = ", ")))
rt_keep <- !all_missing

## Build a long-format tibble with person and item factors.
rt_long <- rt_wide[rt_keep, , drop = FALSE] |>
  as.data.frame() |>
  mutate(person_id = factor(seq_len(sum(rt_keep)))) |>
  pivot_longer(
    cols      = all_of(rt_cols),
    names_to  = "item_id",
    values_to = "rt_sec"
  ) |>
  filter(!is.na(rt_sec), rt_sec > 0) |>
  mutate(
    item_id = factor(item_id, levels = rt_cols),
    log_rt  = log(rt_sec)
  )

message(sprintf(
  "RT dataset: N=%d observations; %d persons; %d items; mean RT = %.1f s",
  nrow(rt_long), nlevels(rt_long$person_id), nlevels(rt_long$item_id),
  mean(rt_long$rt_sec)))

saveRDS(rt_long, "chess_rt_data.rds")

## ---- Priors ----------------------------------------------------------

rt_priors <- c(
  prior(normal(4, 2),       class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(student_t(3, 0, 1), class = "sigma")
)

## ---- Fit -------------------------------------------------------------

chess_rt_fit <- brm(
  formula = log_rt ~ 1 + (1 | person_id) + (1 | item_id),
  data    = rt_long,
  family  = gaussian(),
  prior   = rt_priors,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 1000,
  backend = "cmdstanr",
  seed    = 20260418,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file    = "chess_rt_fit"
)

## ---- Diagnostics -----------------------------------------------------

cat("\n---- brms summary ----\n")
print(summary(chess_rt_fit), digits = 3)

cat("\n---- Posterior variance components ----\n")
post <- as_draws_df(chess_rt_fit)
sp2   <- post$sd_person_id__Intercept^2
sw2   <- post$sd_item_id__Intercept^2
se2   <- post$sigma^2
seta2 <- sp2 + sw2 + se2

qp <- function(x) sprintf("%.3f [%.3f, %.3f]",
                          median(x), quantile(x, 0.025), quantile(x, 0.975))
cat(sprintf("sigma_p^2     : %s\n", qp(sp2)))
cat(sprintf("sigma_w^2     : %s\n", qp(sw2)))
cat(sprintf("sigma_eps^2   : %s\n", qp(se2)))
cat(sprintf("sigma_eta^2   : %s\n", qp(seta2)))

message("\nDone. Next: Rscript 02_fit_acc.R")
