# 03_make_outputs_fims.R -- tables and figure for Section 8.3
suppressPackageStartupMessages({ library(ggplot2) })
r <- readRDS("fims_results.rds")
cell <- function(x, d = 3) sprintf("\\begin{tabular}[t]{@{}c@{}}%s\\\\ {\\footnotesize[%s, %s]}\\end{tabular}",
                                   formatC(x$median, d, format = "f"), formatC(x$lo, d, format = "f"), formatC(x$hi, d, format = "f"))
# One row = 8 cells, matching the 8-column header (Sample, Object, five coefficients, n*).
row <- function(label, object, obj, extra = "") sprintf("%s & %s & %s & %s & %s & %s & %s%s \\\\", label, object,
  cell(obj$icc_eta), cell(obj$icc_Y), cell(obj$icc_Y_rel), cell(obj$icc_I), cell(obj$gap_eta_I), extra)

lines <- c(
  "\\begin{table}[t]", "\\centering",
  "\\caption{Reliability of FIMS mathematics items (Australia, $4{,}320$ students $\\times$ 14 items; Japan, $2{,}051$ students) under a crossed Bernoulli model. Each cell shows the posterior median (top) and 95\\% credible interval (bottom). $\\icc_Y$ is the absolute response-scale ICC (equal to the intraclass kappa) and $\\icc_Y^{\\mathrm{rel}}$ the relative coefficient; $\\icc_I$ estimated by nested Monte Carlo. The two Australian rows fit the same data under the logit and probit links.}",
  "\\label{tab:fims}",
  "\\resizebox{\\textwidth}{!}{\\begin{tabular}{llcccccc}", "\\toprule",
  "Sample (link) & Object & $\\icc_\\eta$ & $\\icc_Y$ & $\\icc_Y^{\\mathrm{rel}}$ & $\\icc_I$ & $\\icc_\\eta - \\icc_I$ & $n^*_Y$ ($c = .80$) \\\\", "\\midrule",
  row("Australia (logit)", "Persons", r$aus_logit$persons,  sprintf(" & %s", r$aus_logit$persons$n80_Y)),
  row("", "Items", r$aus_logit$items, " & "),
  row("Australia (probit)", "Persons", r$aus_probit$persons, sprintf(" & %s", r$aus_probit$persons$n80_Y)),
  row("", "Items", r$aus_probit$items, " & "),
  row("Japan (logit)", "Persons", r$jpn_logit$persons, sprintf(" & %s", r$jpn_logit$persons$n80_Y)),
  row("", "Items", r$jpn_logit$items, " & "),
  "\\bottomrule", "\\end{tabular}}",
  "\\end{table}")
writeLines(lines, "table_fims.tex"); cat("Wrote table_fims.tex\n")

# D-study figure: Australia (logit), persons as object, both scales
cv <- r$aus_logit$curves
df <- rbind(data.frame(n = cv$n, med = cv$eta_med, lo = cv$eta_lo, hi = cv$eta_hi, scale = "Classical (link scale)"),
            data.frame(n = cv$n, med = cv$Y_med,  lo = cv$Y_lo,  hi = cv$Y_hi,  scale = "DGT (response scale: proportion correct)"))
g <- ggplot(df, aes(n, med, colour = scale, fill = scale)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) + geom_line(linewidth = 1) +
  geom_hline(yintercept = c(.7, .8), linetype = 3) + geom_vline(xintercept = 14, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c("#b2182b", "#2166ac")) + scale_fill_manual(values = c("#b2182b", "#2166ac")) +
  labs(x = "Number of items (n)", y = expression(E*rho^2), colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figure_fims_dstudy.pdf", g, width = 6.5, height = 4)
cat("Wrote figure_fims_dstudy.pdf\n")

# headline numbers for the text (read these into Section 8.3; do not retype)
p <- r$aus_logit$persons
cat(sprintf("\nAustralia logit, persons: ICC_eta %.3f [%.3f, %.3f]; ICC_Y %.3f [%.3f, %.3f]; ICC_I %.3f; n*(.80) Y %s vs eta %s; Erho2(14) Y %.3f vs eta %.3f; SEM14 %.3f; ratio Y/eta %.3f\n",
    p$icc_eta$median, p$icc_eta$lo, p$icc_eta$hi, p$icc_Y$median, p$icc_Y$lo, p$icc_Y$hi, p$icc_I$median,
    p$n80_Y, p$n80_eta, p$erho2_14_Y$median, p$erho2_14_eta$median, p$sem14$median, p$ratio_Y_eta$median))
cat(sprintf("Australia probit, persons: ICC_eta %.3f; ICC_Y %.3f; small-effect bound phi(z0)^2/(P(1-P)) %.3f (ratio observed %.3f)\n",
    r$aus_probit$persons$icc_eta$median, r$aus_probit$persons$icc_Y$median, r$aus_probit$persons$bound_d$median, r$aus_probit$persons$ratio_Y_eta$median))
cat(sprintf("Link check: ICC_Y %.3f (logit) vs %.3f (probit); ICC_eta %.3f vs %.3f\n", r$link_check$icc_Y_logit, r$link_check$icc_Y_probit, r$link_check$icc_eta_logit, r$link_check$icc_eta_probit))
writeLines(capture.output(sessionInfo()), "sessionInfo_outputs.txt")
