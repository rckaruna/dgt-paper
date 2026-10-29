# Supplementary repository for "Distributional Generalizability Theory: Reliability for Non-Gaussian Measurements"

Scripts, logs, and script-made outputs for every numerical result in the
paper (Psychometrika manuscript PSY-2026-0092, revised version) and its
Supplementary Material. The companion R package is `dgt`
(https://github.com/rckaruna/dgt); the paper cites version 0.6.1 (tag
`v0.6.1`). The analyses were run with `dgt` 0.5.0 or 0.5.1 as each log
records; the verification scripts show that 0.6.0 and 0.6.1 reproduce the
stored posterior draws of those runs to machine precision (see
`verification/`).

Each folder has its own README or header comment with the exact commands.
Fitted models (`*.rds`) are not tracked; every script refits from the data
it names, and `brms` caches the fit beside the script on first run.

| Folder | Paper | What it holds |
|---|---|---|
| `lexdec/` | Section 9.1 (Illustration 1, lexical decision RTs) | `01_lexdec.R`: lognormal fit to `languageR::lexdec` and the ICCs, overestimation ratio, and D-study |
| `chess/pipeline/` | Section 9.2 (Illustration 2, Amsterdam Chess Test); Tables 3 and 4, Figure 1 | `01_fit_rt.R`, `02_fit_acc.R`, `03_analyze.R`, `04_make_outputs.R`, logs, `sessionInfo.txt`, `chess_numbers.json`, the two tables and the figure |
| `fims_analysis/` | Section 9.3 (Illustration 3, FIMS mathematics items); Table 5, Figure 2. Section 10.4 (hurdle demonstration); Table 6, Figure 3 | `01_fit_fims.R`, `02_analyze_fims.R`, `03_make_outputs_fims.R`; `hurdle_demo.R`, `hurdle_table.R`; logs and `sessionInfo` files; tables and figures |
| `simulations/` | Section 10 (Studies A and B); Appendix C Tables C.1 and C.2; Supplement Section S3 | `preregistration.md` with its sha256 sidecar (frozen before any run; amendments dated at the end), `common.R`, `study_A_lognormal.R`, `study_B_bernoulli.R`, `analyze.R`, `RUN.md`, `logs/`, `results/` (per-condition CSVs, `numbers.json`, the LaTeX tables) |
| `reproduction/vanbelle-fhr/` | Section 8 (reproduction paragraph) and Supplement S5.1 (reproduction of Vanbelle 2026 on her fetal heart rate data) | `fhr_reproduction.R`, `console-log.txt`, `REPORT.md`. Her data and model file are fetched from https://github.com/svanbelle/RELIABILITY.BINARY (AGPL-3.0) and are not copied here |
| `reproduction/huynh-1977/` | Section 2.3 (Choi 2013, Table 4, on Huynh 1977) | `huynh1977.csv` (transcribed from Choi 2013, Figure 1, and checked against Brennan 2001, Table 7.7; see `TRANSCRIPTION.md`), `huynh_reanalysis.R`, `console-log.txt`, `REPORT.md` |
| `reviewer2/` | Supplement S5.2 (Reviewer 2's simulated design) | `log-vanbelle.txt`, `vanbelle_design_results.csv`, `log-chess-guard.txt` |
| `verification/` | Proposition 2 and Remark 3; Appendix A item A11; Table 2; the supplement's tables; package regression checks | `verify_binary_proposition.py` (47 numerical checks on 504 parameter configurations), `verify_appendix_A11.py`, `table1_lognormal_grid.py`, `supplement_tables.py`; `check_dgt_060_regression.R` with `compare_regression_rds.py` and the saved comparisons (`dgt_regression_compare_*.txt`); `check_dgt_060_handball.R` with `dgt_handball_0.6.1.txt` (this one refits a model on data from a separate handball project that is not public; its output is kept for the record) |

Software: R (4.5 or later), `brms`, `cmdstanr`, `posterior`, `dgt`; `TAM`,
`LNIRT`, and `languageR` for the data of Sections 9.3, 9.2, and 9.1; Python 3
with `numpy` and `scipy` for the scripts in `verification/`. Section numbers
refer to the revised manuscript; a few older READMEs inside the folders
carry the submitted version's numbering (Section 8.x for the illustrations).

Fixed states of this repository are tagged by submission round
(`psychometrika-r1` for the revised version). The paper's numbers are read
from the files here and never retyped.

License: MIT (see `LICENSE`) for the scripts and outputs in this
repository. Data remain under their own licenses: `languageR`, `LNIRT`, and
`TAM` as on CRAN; Vanbelle's fetal heart rate data under AGPL-3.0 at her
repository; the Huynh (1977) matrix as printed in Choi (2013) and Brennan
(2001).
