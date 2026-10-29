# 01_fit_fims.R -- Illustration 3: FIMS mathematics items (TAM::data.fims.Aus.Jpn.scored)
# Fits three Bernoulli crossed-random-effects models with brms/cmdstanr and
# caches them with brms's file= argument (rerunning this script loads the
# cache; delete the .rds to refit).
#   AUS-logit : Australian students, logit link   (main analysis)
#   AUS-probit: Australian students, probit link  (same data: link-invariance test)
#   JPN-logit : Japanese students, logit link     (second population)
# Priors: Intercept ~ normal(0, 2.5); sd ~ student_t(3, 0, 1) (both facets).
# Sampler: 4 chains x 2000 iterations (1000 warmup), adapt_delta 0.9.

suppressPackageStartupMessages({ library(TAM); library(brms); library(dplyr); library(tidyr) })
set.seed(20260926)
data("data.fims.Aus.Jpn.scored", package = "TAM")
raw <- data.fims.Aus.Jpn.scored
items <- grep("^M1PTI", names(raw), value = TRUE)
stopifnot(length(items) == 14, nrow(raw) == 6371, !anyNA(raw[, items]))

long <- raw |>
  mutate(person_id = factor(seq_len(n())),
         country   = factor(country, levels = c(1, 2), labels = c("Australia", "Japan"))) |>
  select(person_id, country, all_of(items)) |>
  pivot_longer(all_of(items), names_to = "item_id", values_to = "correct") |>
  mutate(item_id = factor(item_id), correct = as.integer(correct))
saveRDS(long, "fims_long.rds")
cat(sprintf("Long data: %d rows; Australia %d persons, Japan %d persons; 14 items; mean correct %.3f\n",
            nrow(long), n_distinct(long$person_id[long$country == "Australia"]),
            n_distinct(long$person_id[long$country == "Japan"]), mean(long$correct)))

pri <- c(prior(normal(0, 2.5), class = Intercept), prior(student_t(3, 0, 1), class = sd))
fit_one <- function(dat, link, file) {
  brm(correct ~ 1 + (1 | person_id) + (1 | item_id), data = dat,
      family = bernoulli(link = link), prior = pri,
      chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 20260926,
      control = list(adapt_delta = 0.9), backend = "cmdstanr", refresh = 200,
      file = file, file_refit = "never")
}
aus <- long |> filter(country == "Australia") |> droplevels()
jpn <- long |> filter(country == "Japan")     |> droplevels()

fit_aus_logit  <- fit_one(aus, "logit",  "fims_aus_logit_fit")
fit_aus_probit <- fit_one(aus, "probit", "fims_aus_probit_fit")
fit_jpn_logit  <- fit_one(jpn, "logit",  "fims_jpn_logit_fit")

for (nm in c("fit_aus_logit", "fit_aus_probit", "fit_jpn_logit")) {
  f <- get(nm); rh <- brms::rhat(f); nd <- sum(nuts_params(f, pars = "divergent__")$Value)
  cat(sprintf("%-15s max Rhat %.3f  divergences %d  min ESS bulk %.0f\n", nm, max(rh, na.rm = TRUE), nd,
              min(posterior::summarise_draws(f, "ess_bulk")$ess_bulk, na.rm = TRUE)))
}
writeLines(capture.output(sessionInfo()), "sessionInfo_fit.txt")
cat("Done. Next: Rscript 02_analyze_fims.R\n")
