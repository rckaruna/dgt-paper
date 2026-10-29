# Illustration 3 (FIMS) and the hurdle demonstration -- PSY-2026-0092 revision

Run from ~/dev/dgt-revision/fims_analysis/ with dgt >= 0.5.1, brms + cmdstanr,
TAM, dplyr, tidyr, posterior, ggplot2.

    Rscript -e 'for (p in c("TAM","dplyr","tidyr","ggplot2")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")'

## FIMS (Section 8.3)

    nohup sh -c 'Rscript 01_fit_fims.R > log_01.txt 2>&1; Rscript 02_analyze_fims.R > log_02.txt 2>&1; Rscript 03_make_outputs_fims.R > log_03.txt 2>&1' > log_run.txt 2>&1 &
    tail -f log_01.txt        # three fits: Australia logit, Australia probit, Japan logit

Data: TAM::data.fims.Aus.Jpn.scored (6,371 students, 14 dichotomous mathematics
items, no missing cells; 4,320 Australia, 2,051 Japan). Fits are cached as
fims_*_fit.rds (brms file=). Outputs: fims_results.rds, table_fims.tex,
figure_fims_dstudy.pdf, and a headline print at the end of log_03.txt whose
numbers go into Section 8.3 (read from the file, never retyped).

Expected time: 20-40 minutes per fit on an M-series Mac with 4 chains in
parallel; the Australian fits have 4,320 person intercepts.

## Hurdle demonstration (Section 9.4)

    Rscript hurdle_demo.R 2>&1 | tee log_hurdle.txt

Every generating parameter is stated at the top of the script (alpha 0.70,
tau 0.80, mu 3.50, sd_p 0.40, sd_e 0.60, 80 x 60, seed 20260926). Outputs:
hurdle_demo_results.rds, figure_hurdle_dstudy.pdf (more occasions versus
halved intensity noise), and the dgt printouts. The script cross-checks its
closed-form five components against dgt's decomposition at the posterior
means; the two should agree to the third decimal.

## Commit

    cd ~/dev/dgt-revision && printf 'fims_analysis/*.rds\n' >> .gitignore && git add -A && git commit -q -m "D3: FIMS illustration pipeline and hurdle demonstration with stated parameters" && git log --oneline -1
