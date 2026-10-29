## check_dgt_060_regression.R
## Regression check for dgt 0.6.0 against 0.5.1 on the paper's stored
## Bernoulli fits (chess accuracy, FIMS, Vanbelle fetal heart rate).
## Every fit is single-trial with crossed object and facet main effects
## only, so the three 0.6.0 corrections must leave every draw unchanged.
##
## Run twice from ~/dev/dgt-revision, once per version, then compare:
##   Rscript verification/check_dgt_060_regression.R installed   # dgt 0.5.1 from the library
##   Rscript verification/check_dgt_060_regression.R devel       # ~/dev/dgt working tree (0.6.0)
##   Rscript verification/check_dgt_060_regression.R compare
## Outputs: verification/dgt_regression_<mode>.rds and, from compare, the
## printed table, also written to verification/dgt_regression_compare_<a>_vs_<b>.txt.
## verification/compare_regression_rds.py reproduces the compare step without R.

mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) stop("mode: installed | devel | compare")

fits <- list(
  chess_acc_person = list(path = "chess/pipeline/chess_acc_fit.rds",              obj = "person_id"),
  chess_acc_item   = list(path = "chess/pipeline/chess_acc_fit.rds",              obj = "item_id"),
  fims_aus_logit   = list(path = "fims_analysis/fims_aus_logit_fit.rds",          obj = "person_id"),
  fims_aus_probit  = list(path = "fims_analysis/fims_aus_probit_fit.rds",         obj = "person_id"),
  fims_jpn_logit   = list(path = "fims_analysis/fims_jpn_logit_fit.rds",          obj = "person_id"),
  fhr_probit       = list(path = "reproduction/vanbelle-fhr/fit_fhr_probit.rds",  obj = NULL)
)

if (mode %in% c("installed", "devel")) {
  if (mode == "devel") {
    suppressMessages(devtools::load_all("~/dev/dgt", quiet = TRUE))
  } else {
    library(dgt)
  }
  ver <- as.character(packageVersion("dgt"))
  message("dgt version in use: ", ver)
  out <- list(version = ver)
  for (nm in names(fits)) {
    f <- readRDS(fits[[nm]]$path)
    obj <- fits[[nm]]$obj
    if (is.null(obj)) obj <- names(brms::ranef(f))[1]
    set.seed(20261010)
    r <- suppressMessages(dgt::dgt_icc(f, person_group = obj))
    d <- r$draws
    out[[nm]] <- d[, intersect(c("icc_Y_abs", "icc_Y_rel", "icc_eta", "icc_I",
                                 "var_obj", "var_facet", "var_int", "mean_pi"), names(d))]
    message(sprintf("%-18s median ICC_Y_abs %.4f  rel %.4f  eta %.4f  I %.4f",
                    nm, median(d$icc_Y_abs), median(d$icc_Y_rel),
                    median(d$icc_eta), median(d$icc_I)))
  }
  saveRDS(out, sprintf("verification/dgt_regression_%s.rds", mode))
  message("wrote verification/dgt_regression_", mode, ".rds")
}

if (mode == "compare") {
  a <- readRDS("verification/dgt_regression_installed.rds")
  b <- readRDS("verification/dgt_regression_devel.rds")
  out_txt <- sprintf("verification/dgt_regression_compare_%s_vs_%s.txt", a$version, b$version)
  sink(out_txt, split = TRUE)
  cat("installed ", a$version, " vs devel ", b$version, "\n", sep = "")
  ok <- TRUE
  for (nm in names(fits)) {
    for (col in names(a[[nm]])) {
      d <- max(abs(a[[nm]][[col]] - b[[nm]][[col]]))
      flag <- if (d < 1e-10) "identical" else sprintf("DIFFERS (max |diff| = %.3g)", d)
      if (d >= 1e-10) ok <- FALSE
      cat(sprintf("%-18s %-10s %s\n", nm, col, flag))
    }
  }
  cat(if (ok) "\nAll draws identical between the two versions.\n" else "\nDIFFERENCES FOUND.\n")
  sink()
  message("wrote ", out_txt)
}
