# hurdle_table.R -- write table_hurdle.tex from hurdle_demo_results.rds (Section 9.4)
r <- readRDS("hurdle_demo_results.rds")
s <- r$icc$summary; v <- r$var
get <- function(pat) { i <- grep(pat, s$measure, fixed = TRUE)[1]; s[i, ] }
f3 <- function(x) sprintf("%.3f", x)
ci <- function(row) sprintf("%s [%s, %s]", f3(row$estimate), f3(row$lower), f3(row$upper))
comp <- v$components   # columns: component, estimate, mean, sd, lower, upper
comp$share <- comp$estimate / sum(comp$estimate)
p <- r$params
lines <- c("\\begin{table}[t]", "\\centering",
  sprintf("\\caption{Hurdle variance decomposition for the simulated practice-minutes data ($%d$ persons $\\times$ $%d$ occasions; generating values $\\alpha = %.2f$, $\\tau = %.2f$, $\\mu = %.2f$, $\\sigma_p = %.2f$, $\\sigma_e = %.2f$, seed %d). Posterior medians with 95\\%% credible intervals.}",
          p[["N_P"]], p[["N_J"]], p[["ALPHA"]], p[["TAU"]], p[["MU"]], p[["SD_P"]], p[["SD_E"]], p[["SEED"]]),
  "\\label{tab:hurdle}", "\\begin{tabular}{lrr}", "\\toprule", "Component & Value [95\\% CI] & \\% of total \\\\", "\\midrule")
for (i in seq_len(nrow(comp))) lines <- c(lines, sprintf("%s & %.2f [%.2f, %.2f] & %.1f\\%% \\\\", comp$component[i], comp$estimate[i], comp$lower[i], comp$upper[i], 100 * comp$share[i]))
lines <- c(lines, "\\midrule", "\\multicolumn{3}{l}{Reliability} \\\\",
  sprintf("$\\icc_Z$ (engagement) & \\multicolumn{2}{l}{%s} \\\\", ci(get("ICC_Z"))),
  sprintf("$\\icc_{Y^*}$ (intensity) & \\multicolumn{2}{l}{%s} \\\\", ci(get("ICC_Y*"))),
  sprintf("$\\icc^{\\mathrm{comp}}$ (composite) & \\multicolumn{2}{l}{%s} \\\\", ci(get("ICC_comp"))),
  sprintf("$\\icc_I$ (information) & \\multicolumn{2}{l}{%s} \\\\", ci(r$inf$summary[grep("ICC_I", r$inf$summary$measure)[1], ])),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, "table_hurdle.tex"); cat("Wrote table_hurdle.tex\n"); cat(lines, sep = "\n")
