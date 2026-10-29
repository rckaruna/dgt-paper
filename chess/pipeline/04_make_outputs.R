## 04_make_outputs.R
## ---------------------------------------------------------------------
## Generate LaTeX tables and figures for Section 8.2 from cached
## chess_results.rds.
##
## Usage:  Rscript 04_make_outputs.R
## ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

results <- readRDS("chess_results.rds")

## =====================================================================
## Table 1: Response-time ICCs
## =====================================================================

## Two-line cells: posterior median above, 95% interval below in footnotesize
## (same layout as fims_analysis/03_make_outputs_fims.R), so that the
## tables fit the text width without shrinking.
fmt <- function(x) {
  sprintf("\\begin{tabular}[t]{@{}c@{}}%.3f\\\\ {\\footnotesize[%.3f, %.3f]}\\end{tabular}",
          x$median, x$lower, x$upper)
}

fmt_O <- function(x) {
  sprintf("\\begin{tabular}[t]{@{}c@{}}%.2f\\\\ {\\footnotesize[%.2f, %.2f]}\\end{tabular}",
          x$median, x$lower, x$upper)
}

r <- results$rt

## Row 1: Person as object
row_person <- sprintf(
  "Person & %s & %s & %s & %s & %s & %s & %s & %d \\\\",
  fmt(r$variance_components$sigma_p2),        # sigma_p^2
  fmt(r$variance_components$sigma_w2),        # sigma_w^2 (secondary facet)
  fmt(r$variance_components$sigma_eps2),      # sigma_eps^2
  fmt(r$variance_components$sigma_eta2),      # sigma_eta^2
  fmt(r$icc_person$icc_eta),
  fmt(r$icc_person$icc_Y),
  fmt_O(r$icc_person$O),
  r$dstudy$person_80$n_star_median
)

## Row 2: Item as object (rearranges which variance is "primary")
row_item <- sprintf(
  "Item & %s & %s & %s & %s & %s & %s & %s & %d \\\\",
  fmt(r$variance_components$sigma_w2),        # sigma_w^2 (primary, item var)
  fmt(r$variance_components$sigma_p2),        # sigma_p^2 (secondary)
  fmt(r$variance_components$sigma_eps2),
  fmt(r$variance_components$sigma_eta2),
  fmt(r$icc_item$icc_eta),
  fmt(r$icc_item$icc_Y),
  fmt_O(r$icc_item$O),
  r$dstudy$item_80$n_star_median
)

rt_table <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Response-Time Reliability for the Amsterdam Chess Test, Both Objects of Measurement. Entries are posterior medians with 95\\% credible intervals.}",
  "\\label{tab:chess-rt}",
  "\\small",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lcccccccr}",
  "\\toprule",
  "Object & $\\sigma^2_{\\text{obj}}$ & $\\sigma^2_{\\text{facet}}$ & $\\sigma^2_\\epsilon$ & $\\sigma^2_\\eta$ & $\\icc_\\eta$ & $\\icc_Y^{\\mathrm{LN}}$ & $\\mathcal{O}$ & $n^*_Y$ ($c=.80$) \\\\",
  "\\midrule",
  row_person,
  row_item,
  "\\bottomrule",
  "\\end{tabular}}",
  "\\begin{flushleft}",
  "\\small\\textit{Note.} $\\sigma^2_{\\text{obj}}$ is the variance of the object of measurement (persons in the first row, items in the second) and $\\sigma^2_{\\text{facet}}$ that of the other facet; $n^*_Y$ is the number of conditions required for $E\\rho^2 \\ge .80$ on the response scale (posterior median of the draw-wise requirement).",
  "\\end{flushleft}",
  "\\end{table}"
)

writeLines(rt_table, "table_chess_rt.tex")
message("Wrote table_chess_rt.tex")

## =====================================================================
## Table 2: Accuracy ICCs
## =====================================================================

a <- results$acc

row_acc_person <- sprintf(
  "Person & %s & %s & %s & %s & %s & %s & %g (%g) \\\\",
  fmt(a$variance_components$tau_p2),
  fmt(a$variance_components$tau_w2),
  fmt(a$icc_person$icc_eta),
  fmt(a$icc_person$icc_Y),
  fmt(a$icc_person$icc_I),
  fmt(a$icc_person$info_loss),
  a$dstudy$person_80_Y$n_star_median,
  a$dstudy$person_80_eta$n_star_median
)

row_acc_item <- sprintf(
  "Item & %s & %s & %s & %s & %s & %s & %g (%g) \\\\",
  fmt(a$variance_components$tau_w2),
  fmt(a$variance_components$tau_p2),
  fmt(a$icc_item$icc_eta),
  fmt(a$icc_item$icc_Y),
  fmt(a$icc_item$icc_I),
  fmt(a$icc_item$info_loss),
  a$dstudy$item_80_Y$n_star_median,
  a$dstudy$item_80_eta$n_star_median
)

acc_table <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Accuracy Reliability for the Amsterdam Chess Test, Illustrating Theorem~\\ref{thm:discrete-loss} and Proposition~\\ref{prop:binary-icc}. Logit-scale variance components; $\\icc_Y$ is the absolute response-scale ICC on the 0/1 scale (equal to the population intraclass kappa), computed by Gauss--Hermite quadrature from posterior draws; $\\icc_I$ estimated by Monte Carlo from posterior draws; $n^*_Y$ is the number of conditions required for $E\\rho^2 \\ge .80$ for the observed proportion, with the link-scale requirement in parentheses.}",
  "\\label{tab:chess-acc}",
  "\\small",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lccccccc}",
  "\\toprule",
  "Object & $\\tau^2_{\\text{obj}}$ & $\\tau^2_{\\text{facet}}$ & $\\icc_\\eta$ & $\\icc_Y$ & $\\icc_I$ & $\\icc_\\eta - \\icc_I$ & $n^*_Y$ ($n^*_\\eta$) \\\\",
  "\\midrule",
  row_acc_person,
  row_acc_item,
  "\\bottomrule",
  "\\end{tabular}}",
  "\\end{table}"
)

writeLines(acc_table, "table_chess_acc.tex")
message("Wrote table_chess_acc.tex")

## =====================================================================
## Figure: D-study curves
## =====================================================================

curves <- results$rt$dstudy_curves

plot_df <- bind_rows(
  curves$person_classical |> mutate(object = "Person", method = "Classical (link-scale)"),
  curves$person_dgt       |> mutate(object = "Person", method = "DGT (response-scale)"),
  curves$item_classical   |> mutate(object = "Item",   method = "Classical (link-scale)"),
  curves$item_dgt         |> mutate(object = "Item",   method = "DGT (response-scale)")
)

plot_df <- plot_df |> mutate(object = factor(object, levels = c("Person", "Item")))

p <- ggplot(plot_df, aes(x = n, y = median, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0.70, linetype = "dotted", color = "grey40") +
  facet_wrap(~ object) +
  scale_color_manual(values = c("Classical (link-scale)" = "#b22222",
                                "DGT (response-scale)"   = "#1f3a5f")) +
  scale_fill_manual(values = c("Classical (link-scale)" = "#b22222",
                               "DGT (response-scale)"   = "#1f3a5f")) +
  labs(
    x = "Number of items (n)",
    y = expression(E*rho^2),
    color = NULL, fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

ggsave("figure_chess_dstudy.pdf", p, width = 6.5, height = 3.2)
message("Wrote figure_chess_dstudy.pdf")

## =====================================================================
## Numbers quoted in the manuscript text (Section 8.2), written by script
## =====================================================================

a <- results$acc; r <- results$rt
num <- list(
  dgt_version = a$dgt_version,
  rt = list(person = r$icc_person, item = r$icc_item,
            n80_person = r$dstudy$person_80$n_star_median, n70_person = r$dstudy$person_70$n_star_median,
            n80_item = r$dstudy$item_80$n_star_median, n70_item = r$dstudy$item_70$n_star_median),
  acc = list(person = a$icc_person, item = a$icc_item, dstudy = a$dstudy),
  data = results$data_summary
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(num, "chess_numbers.json", auto_unbox = TRUE, pretty = TRUE, digits = 6)
  message("Wrote chess_numbers.json")
} else {
  dput(num, file = "chess_numbers.R")
  message("jsonlite not installed; wrote chess_numbers.R instead")
}

## =====================================================================
## sessionInfo for reproducibility
## =====================================================================

sink("sessionInfo.txt")
print(sessionInfo())
sink()
message("Wrote sessionInfo.txt")

message("\nAll outputs generated. Files for paper:")
message("  - table_chess_rt.tex")
message("  - table_chess_acc.tex")
message("  - figure_chess_dstudy.pdf")
