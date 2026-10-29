## 01_lexdec.R -- Illustration 1 (Section 9.1): lexical decision reaction times.
## Reproduces the numbers of Section 9.1 of PSY-2026-0092 from the languageR
## lexdec data with brms and dgt. Run from ~/dev/dgt-revision/lexdec/:
##   Rscript 01_lexdec.R 2>&1 | tee log_01.txt
## Every number quoted in Section 9.1 is printed at the end under "HEADLINE".
## The fit is cached in lexdec_fit.rds (ignored by git); delete it to refit.

suppressMessages({ library(brms); library(dgt); library(posterior) })
if (!requireNamespace("languageR", quietly = TRUE))
  install.packages("languageR", repos = "https://cloud.r-project.org")
cat("dgt version:", as.character(packageVersion("dgt")), "\n")

d <- languageR::lexdec
d <- d[d$Correct == "correct", c("Subject", "Word", "RT")]
d$rt_ms <- exp(d$RT)             # lexdec$RT is the natural log of the RT in ms
cat(sprintf("%d correct-response observations; %d subjects; %d words; mean RT %.0f ms, median %.0f ms\n",
            nrow(d), nlevels(droplevels(d$Subject)), nlevels(droplevels(d$Word)),
            mean(d$rt_ms), median(d$rt_ms)))

fit <- brm(rt_ms ~ 1 + (1 | Subject) + (1 | Word), data = d, family = lognormal(),
           chains = 4, iter = 4000, warmup = 1000, seed = 20260429,
           backend = "cmdstanr", refresh = 0, file = "lexdec_fit")
print(summary(fit))

q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .025), quantile(x, .975))
q3 <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x), quantile(x, .025), quantile(x, .975))
dr <- as_draws_df(fit)
sd_p <- dr$sd_Subject__Intercept; sd_w <- dr$sd_Word__Intercept; sd_e <- dr$sigma

set.seed(20260429)
icc <- dgt_icc(fit, person_group = "Subject")
print(icc)
ds <- dgt_dstudy(fit, person_group = "Subject", n_grid = 1:40)
rq <- ds$required_n
print(rq)

cat("\nHEADLINE (Section 9.1)\n")
cat("sigma_p      ", q(sd_p), "\n")
cat("sigma_w      ", q(sd_w), "\n")
cat("sigma_eps    ", q(sd_e), "\n")
cat("sigma^2_eta  ", sprintf("%.3f", median(sd_p^2 + sd_w^2 + sd_e^2)), "\n")
cat("ICC_eta      ", q3(icc$draws$icc_eta), "\n")
cat("ICC_Y        ", q3(icc$draws$icc_Y), "\n")
cat("O (median of overestimation draws)", sprintf("%.3f", median(icc$draws$overestimation)), "\n")
cat("n* for 0.80: link", rq$link_scale[rq$threshold == 0.80],
    " response", rq$arith_mean[rq$threshold == 0.80], "\n")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
