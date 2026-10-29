## check_dgt_060_handball.R
## Refits the nine-meter shooting efficiency models of the 7 Oct 2026
## handball reanalysis and checks dgt 0.6.x (working tree ~/dev/dgt)
## against the numbers that exposed the three 0.5.1 defects:
##   relative coefficient with the team:game term  0.150 (0.5.1 gave 0.175)
##   matches for 0.80, with the term                about 23
##   matches for 0.80, without the term             about 17 (reanalysis table;
##                                                  Spearman-Brown on the relative
##                                                  coefficient; absolute gives 20)
## and, for 0.6.1, the fourth defect: the link-scale curve of 0.6.0 treated
## each match as one shot (0.70 only at 162 matches); 0.6.1 divides the latent
## residual by the trials, so the link-scale required n is now of the same
## order as the response-scale one.
## Console output is also written to verification/dgt_handball_<version>.txt.
## Data: outputs/cache/clean.rds of the handball G-study project.
## Run from ~/dev/dgt-revision:
##   Rscript verification/check_dgt_060_handball.R
## Saves the two brms fits beside the data so the package test can be
## pointed at them with DGT_HANDBALL_FIT.

suppressMessages(devtools::load_all("~/dev/dgt", quiet = TRUE))
stopifnot(packageVersion("dgt") >= "0.6.0")
library(brms)

hb_dir <- path.expand("~/Documents/research/handball_tier1_gstudy")
dat <- readRDS(file.path(hb_dir, "outputs/cache/clean.rds"))
d <- data.frame(team = factor(dat$team), game = factor(dat$game),
                num = dat$g9, den = dat$s9)
d <- d[d$den > 0 & !is.na(d$num) & !is.na(d$den), ]
message(sprintf("%d team-games, %d teams, %d games; shots per match %g to %g (harmonic mean %.2f)",
                nrow(d), nlevels(d$team), nlevels(d$game), min(d$den), max(d$den), 1 / mean(1 / d$den)))

fit_one <- function(formula, file) {
  if (file.exists(file)) return(readRDS(file))
  f <- brm(formula, data = d, family = binomial(), chains = 4, iter = 1500,
           seed = 2026, refresh = 0, silent = 2)
  saveRDS(f, file); f
}
out_dir <- file.path(hb_dir, "outputs/cache")
fit_with    <- fit_one(num | trials(den) ~ 1 + (1 | team) + (1 | game) + (1 | team:game),
                       file.path(out_dir, "dgt060_nine_meter_with_teamgame.rds"))
fit_without <- fit_one(num | trials(den) ~ 1 + (1 | team) + (1 | game),
                       file.path(out_dir, "dgt060_nine_meter_without_teamgame.rds"))

report <- function(fit, label) {
  set.seed(2026)
  icc <- dgt_icc(fit, person_group = "team", K = 200)
  ds  <- dgt_dstudy(fit, n_grid = 1:200, person_group = "team", K = 200)
  rq  <- ds$required_n
  cat(sprintf("\n== %s ==\n", label))
  print(icc)
  cat(sprintf("required matches (median curve): 0.70 -> %s, 0.80 -> %s (response scale); link scale 0.70 -> %s, 0.80 -> %s\n",
              rq$response_scale[rq$threshold == 0.70], rq$response_scale[rq$threshold == 0.80],
              rq$link_scale[rq$threshold == 0.70], rq$link_scale[rq$threshold == 0.80]))
  invisible(list(icc = icc, ds = ds))
}
out_txt <- sprintf("verification/dgt_handball_%s.txt", packageVersion("dgt"))
sink(out_txt, split = TRUE)
cat(sprintf("dgt %s (%s), run %s\n", packageVersion("dgt"),
            if (nzchar(Sys.getenv("DGT_COMMIT"))) Sys.getenv("DGT_COMMIT") else "~/dev/dgt working tree",
            format(Sys.time(), "%Y-%m-%d %H:%M")))
r_with    <- report(fit_with,    "nine-meter efficiency, with team:game")
r_without <- report(fit_without, "nine-meter efficiency, without team:game")
sink()
message("wrote ", out_txt)

cat("\nPackage test with the saved fit:\n")
cat(sprintf("DGT_HANDBALL_FIT=%s Rscript -e 'devtools::test(\"~/dev/dgt\", filter = \"interaction\")'\n",
            file.path(out_dir, "dgt060_nine_meter_with_teamgame.rds")))
