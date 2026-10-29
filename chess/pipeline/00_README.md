# Section 8.2 Analysis — Amsterdam Chess Test

Reproducible pipeline for the DGT paper's Illustration 2: response-time
and accuracy reliability in a psychometric chess expertise test.

## Data

`AmsterdamChess` from the LNIRT R package (van der Maas & Wagenmakers,
2005). 259 chess players × 40 chess-tactics items; each cell gives a
response time (seconds) and a correctness score (0/1). Missing RTs are
coded 10000 and missing responses are coded 9; three persons (rows 147,
201, 209) are entirely missing and are dropped.

After preprocessing the design is approximately 256 × 40 = 10,240 cells
with item-level missingness around 3–5% of cells.

## Model

### Response times (Section 8.2.1)

Lognormal crossed random-effects:

    log(RT_{pi}) = mu + nu_p + omega_i + eps_{pi}
    nu_p   ~ N(0, sigma_p^2)
    omega_i ~ N(0, sigma_w^2)
    eps_pi ~ N(0, sigma_eps^2)

Fitted in `brms` with `cmdstanr`, weakly informative penalized-complexity
priors on SD parameters.

### Accuracy (Section 8.2.2)

Binomial (Bernoulli per cell) crossed random-effects:

    Y_{pi}           ~ Bernoulli(pi_{pi})
    logit(pi_{pi})   = alpha + u_p + v_i
    u_p ~ N(0, tau_p^2), v_i ~ N(0, tau_w^2)

The discrete sampling step makes this a Theorem 5.5 example: ICC_I <
ICC_eta strictly.

## Pipeline

    Rscript 01_fit_rt.R         # overnight (~15 min)
    Rscript 02_fit_acc.R        # overnight (~20 min)
    Rscript 03_analyze.R        # lightweight; produces chess_results.rds
    Rscript 04_make_outputs.R   # generates LaTeX tables and figures

## Outputs

    chess_rt_fit.rds          brms fit, response-time model
    chess_acc_fit.rds         brms fit, accuracy model
    chess_results.rds         list of posterior summaries, ICC computations,
                              D-study curves
    table_chess_rt.tex        LaTeX table for Section 8.2.1
    table_chess_acc.tex       LaTeX table for Section 8.2.2
    figure_chess_dstudy.pdf   D-study curves (optional)

## Reproducibility

R 4.4+, brms 2.21+, cmdstanr with cmdstan 2.35+. All seeds set to
20260418. `sessionInfo()` is captured at the end of `04_make_outputs.R`.
