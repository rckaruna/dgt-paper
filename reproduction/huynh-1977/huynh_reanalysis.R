## Re-analysis of the Huynh (1977) 12 x 6 dichotomous data set printed as
## Figure 1 on p. 25 of Choi (2013), with dgt 0.5.0.
##
## Target of comparison -- Choi (2013), Table 4, printed p. 26, computed at
## I' = 6 items (her Note c):
##
##            Method A   Method B   Method C   Method D
##   Person     0.047      0.052      0.052      0.978
##   Item       0.012      0.016      0.057      0.246
##   PxItem     0.184      0.179      0.248      0.000
##   Error                            0.187      3.287
##   E(rho^2)   0.605      0.635      0.418      0.641
##   Phi        0.560      0.615      0.388      0.624
##
##   Method C = GIRM (Briggs & Wilson 2007), response scale.
##   Method D = logistic mixed model, latent (logit) scale.
##
## Nothing in the dgt package and no .tex file is modified by this script.

options(width = 120, digits = 6)

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dgt)
  library(posterior)
})

cat("=====================================================================\n")
cat("Huynh (1977) re-analysis, data from Choi (2013) Figure 1, p. 25\n")
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

wide <- read.csv("huynh1977.csv")
cat("--- transcribed matrix (see TRANSCRIPTION.md) ---\n")
print(wide, row.names = FALSE)
cat("\npersons:", nrow(wide), "  items:", ncol(wide) - 1L, "\n")
cat("missing cells:", sum(is.na(wide[, -1])),
    "  observed responses:", sum(!is.na(wide[, -1])), "\n")
stopifnot(nrow(wide) == 12L, ncol(wide) == 7L, sum(is.na(wide[, -1])) == 13L)
cat("CHECK: 12 persons x 6 items, 13 missing, 59 observed -- OK\n\n")

long <- reshape(wide, direction = "long", varying = paste0("i", 1:6),
                v.names = "y", timevar = "item", idvar = "person",
                times = 1:6)
long <- long[!is.na(long$y), ]
long$person <- factor(long$person)
long$item   <- factor(long$item)
rownames(long) <- NULL

cat("long data frame:", nrow(long), "rows\n")
cat("responses per person:\n"); print(table(long$person))
cat("responses per item:\n");   print(table(long$item))
cat("overall proportion correct:", round(mean(long$y), 4), "\n\n")

## ---------------------------------------------------------------------------
## 2. Model
## ---------------------------------------------------------------------------
## y ~ 1 + (1 | person) + (1 | item), bernoulli(logit), weakly informative
## priors. With one binary response per cell the person-by-item interaction is
## confounded with binomial error, which is why Choi's Method D reports
## person x item = 0.000 and error 3.287 (= pi^2/3, the logistic residual
## variance). Our model is the same specification.

prior_h <- c(
  prior(normal(0, 2),          class = "Intercept"),
  prior(student_t(3, 0, 2.5),  class = "sd", group = "person"),
  prior(student_t(3, 0, 2.5),  class = "sd", group = "item")
)

cat("--- brms model and priors ---\n")
print(prior_h)
cat("\n")

fit <- brm(
  y ~ 1 + (1 | person) + (1 | item),
  data = long,
  family = bernoulli(link = "logit"),
  prior = prior_h,
  backend = "cmdstanr",
  chains = 4, iter = 3000, warmup = 1000,
  cores = 4, seed = 19770101,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  refresh = 0,
  file = "fit_huynh_logit", file_refit = "on_change"
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

pd <- posterior::as_draws_df(fit)
sd_p <- pd$sd_person__Intercept
sd_i <- pd$sd_item__Intercept

q <- function(x) c(median = median(x), lo = unname(quantile(x, 0.025)),
                   hi = unname(quantile(x, 0.975)))
fmt <- function(x) sprintf("%.3f (%.3f, %.3f)", median(x),
                           quantile(x, 0.025), quantile(x, 0.975))

cat("--- variance components, logit scale (residual fixed at pi^2/3 =",
    round(pi^2/3, 4), ") ---\n")
cat("sigma^2 person :", fmt(sd_p^2), "   Choi Method D: 0.978\n")
cat("sigma^2 item   :", fmt(sd_i^2), "   Choi Method D: 0.246\n")
cat("sigma^2 error  : pi^2/3 = 3.290 (fixed)  Choi Method D: 3.287\n\n")

## ---------------------------------------------------------------------------
## 3. dgt_icc: single-item coefficients
## ---------------------------------------------------------------------------

cat("--- dgt::dgt_icc(fit, person_group = \"person\") ---\n")
icc <- dgt::dgt_icc(fit, person_group = "person", seed = 19770101)
print(icc)
cat("\nsummary table:\n")
print(icc$summary)

cat("\nposterior summaries of the per-draw ICCs [median, 2.5%, 97.5%]:\n")
for (nm in c("icc_eta", "icc_Y_abs", "icc_Y_rel", "icc_I")) {
  cat(sprintf("  %-10s : %s\n", nm, fmt(icc$draws[[nm]])))
}

cat("\ncheck that ICC_eta uses the pi^2/3 convention:\n")
icc_eta_manual <- sd_p^2 / (sd_p^2 + sd_i^2 + pi^2/3)
cat("  max |dgt icc_eta - sd_p^2/(sd_p^2 + sd_i^2 + pi^2/3)| =",
    max(abs(icc$draws$icc_eta - icc_eta_manual)), "\n\n")

## ---------------------------------------------------------------------------
## 4. dgt_dstudy: coefficients at 6 items
## ---------------------------------------------------------------------------

cat("--- dgt::dgt_dstudy(fit, n_grid = c(1:12, 25, 50), person_group = \"person\") ---\n")
ds <- dgt::dgt_dstudy(fit, n_grid = c(1:12, 25, 50), person_group = "person",
                      seed = 19770101)
cat("curve columns:", paste(names(ds$curves), collapse = ", "), "\n\n")
print(ds$curves[ds$curves$n %in% c(1, 6, 12, 25, 50), ], row.names = FALSE)

## Choi (pp. 27-28) reports that the vertical gap between Method C and the other
## methods "does not disappear as the D study item sample size increases to
## many, such as 25 items". The DGT scale gap must close as items accumulate,
## because both coefficients tend to 1. Tabulate it.
cat("\n--- the scale gap as items accumulate (dgt_dstudy medians) ---\n")
gap_tab <- do.call(rbind, lapply(c(1, 6, 12, 25, 50), function(nn) {
  lk <- ds$curves$median[ds$curves$n == nn & ds$curves$type == "link-scale"]
  rp <- ds$curves$median[ds$curves$n == nn & ds$curves$type == "response-scale"]
  data.frame(n = nn, link = round(lk, 3), response = round(rp, 3),
             gap = round(lk - rp, 3), ratio = round(rp / lk, 3))
}))
print(gap_tab, row.names = FALSE)

cat("\n--- the two coefficients at n = 6 items ---\n")
d6 <- ds$curves[ds$curves$n == 6, ]
print(d6, row.names = FALSE)

## dgt_dstudy's curves put the item main-effect variance in the error term
## (link scale: sigma_p^2 / (sigma_p^2 + (sigma_i^2 + pi^2/3)/n)), i.e. they are
## absolute / dependability coefficients, Choi's Phi. The relative /
## generalizability counterpart, Choi's E(rho^2), drops sigma_i^2 from the
## denominator. It is recomputed here by hand from the same posterior draws;
## dgt_dstudy() does not return it.

n_d <- 6
link_abs_6 <- sd_p^2 / (sd_p^2 + (sd_i^2 + pi^2/3) / n_d)
link_rel_6 <- sd_p^2 / (sd_p^2 + (pi^2/3) / n_d)

cat("\nhand check of the dgt link-scale curve at n = 6:\n")
cat("  dgt link-scale estimate      :",
    sprintf("%.4f", d6$median[d6$type == "link-scale"]), "\n")
cat("  sd_p^2/(sd_p^2+(sd_i^2+pi^2/3)/6) posterior median:",
    sprintf("%.4f", median(link_abs_6)), "\n")

## Response-scale relative coefficient at n = 6, by the same Gauss-Hermite
## machinery dgt uses at n = 1, with the item main effect removed from the
## error term.
gh <- dgt:::.gh_nodes(40)      # dgt's own nodes, read-only use of the package
resp_rel_6 <- resp_abs_6 <- numeric(length(sd_p))
alpha <- pd$b_Intercept
for (s in seq_along(sd_p)) {
  z <- gh$x; w <- gh$w
  eta <- outer(alpha[s] + sd_p[s] * z, sd_i[s] * z, "+")
  pg  <- plogis(eta)
  mu_p <- as.vector(pg %*% w)
  mu_i <- as.vector(w %*% pg)
  m    <- sum(w * mu_p)
  v_p  <- sum(w * (mu_p - m)^2)
  v_i  <- sum(w * (mu_i - m)^2)
  v_t  <- sum(outer(w, w) * (pg - m)^2)
  v_int <- max(v_t - v_p - v_i, 0)
  e_b  <- sum(outer(w, w) * pg * (1 - pg))
  resp_abs_6[s] <- v_p / (v_p + (v_i + v_int + e_b) / n_d)
  resp_rel_6[s] <- v_p / (v_p + (v_int + e_b) / n_d)
}

cat("\n=====================================================================\n")
cat("COEFFICIENTS AT 6 ITEMS\n")
cat("=====================================================================\n")
cat("absolute (Phi), link  scale, dgt_dstudy :",
    sprintf("%.3f  [%.3f, %.3f]", d6$median[d6$type == "link-scale"],
            d6$lower[d6$type == "link-scale"],
            d6$upper[d6$type == "link-scale"]),
    "   Choi Method D Phi = 0.624\n")
cat("absolute (Phi), resp. scale, dgt_dstudy :",
    sprintf("%.3f  [%.3f, %.3f]", d6$median[d6$type == "response-scale"],
            d6$lower[d6$type == "response-scale"],
            d6$upper[d6$type == "response-scale"]),
    "   Choi Method C Phi = 0.388\n")
cat("relative (E rho^2), link  scale, by hand:", fmt(link_rel_6),
    "   Choi Method D E(rho^2) = 0.641\n")
cat("relative (E rho^2), resp. scale, by hand:", fmt(resp_rel_6),
    "   Choi Method C E(rho^2) = 0.418\n")
cat("(hand response-scale absolute at n=6, as a check on dgt_dstudy:",
    fmt(resp_abs_6), ")\n")
cat("=====================================================================\n\n")

## ---------------------------------------------------------------------------
## 5. Assembled comparison
## ---------------------------------------------------------------------------

comparison <- data.frame(
  coefficient = c("E(rho^2), latent/logit scale", "Phi, latent/logit scale",
                  "E(rho^2), response scale",     "Phi, response scale"),
  ours_at_6_items = c(fmt(link_rel_6),
                      sprintf("%.3f (%.3f, %.3f)",
                              d6$median[d6$type == "link-scale"],
                              d6$lower[d6$type == "link-scale"],
                              d6$upper[d6$type == "link-scale"]),
                      fmt(resp_rel_6),
                      sprintf("%.3f (%.3f, %.3f)",
                              d6$median[d6$type == "response-scale"],
                              d6$lower[d6$type == "response-scale"],
                              d6$upper[d6$type == "response-scale"])),
  choi_method_D_latent = c("0.641", "0.624", "--", "--"),
  choi_method_C_response = c("--", "--", "0.418", "0.388"),
  stringsAsFactors = FALSE
)
cat("--- comparison table ---\n")
print(comparison, right = FALSE, row.names = FALSE)

cat("\n--- does the response-scale coefficient land near Method C? ---\n")
inCI <- function(target, draws)
  target >= quantile(draws, 0.025) && target <= quantile(draws, 0.975)
cat("Choi's Method C E(rho^2) = 0.418 inside our response-scale relative 95% CI:",
    inCI(0.418, resp_rel_6), "\n")
cat("Choi's Method C Phi      = 0.388 inside our response-scale absolute 95% CI:",
    inCI(0.388, resp_abs_6), "\n")
cat("Choi's Method D E(rho^2) = 0.641 inside our link-scale relative 95% CI:",
    inCI(0.641, link_rel_6), "\n")
cat("Choi's Method D Phi      = 0.624 inside our link-scale absolute 95% CI:",
    inCI(0.624, link_abs_6), "\n")
cat("posterior P(response-scale relative < link-scale relative) =",
    mean(resp_rel_6 < link_rel_6), "\n")

saveRDS(list(icc = icc, dstudy = ds, comparison = comparison,
             link_rel_6 = link_rel_6, link_abs_6 = link_abs_6,
             resp_rel_6 = resp_rel_6, resp_abs_6 = resp_abs_6),
        "results_huynh.rds")

cat("\ndone.\n")
cat("dgt version in this log:", as.character(packageVersion("dgt")), "\n")
