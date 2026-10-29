# Running the pre-registered simulations (Studies A and B)

Everything below runs in ~/dev/dgt-revision/simulations. Requires dgt >= 0.5.1
installed (remotes::install_github("rckaruna/dgt", ref = "v0.5.1")), brms with
cmdstanr, posterior, mvtnorm, jsonlite.

## 0. Freeze the pre-registration (once, before any run)

    cd ~/dev/dgt-revision/simulations
    shasum -a 256 preregistration.md > preregistration.md.sha256
    cat preregistration.md.sha256
    cd ~/dev/dgt-revision && git add -A && git commit -q -m "simulations: pre-registration frozen (sha256 recorded) and scripts" && git log --oneline -1

The sha256 is re-checked before analysis: `shasum -a 256 -c preregistration.md.sha256`.

## 1. Pilot (3 replications per cell; timing and sanity)

    cd ~/dev/dgt-revision/simulations
    Rscript study_A_lognormal.R --pilot --workers 4 2>&1 | tee logs/A_pilot.log
    Rscript study_B_bernoulli.R --pilot --workers 4 2>&1 | tee logs/B_pilot.log
    Rscript analyze.R --pilot

Check in the logs: seconds per fit (expect 10-20 s for A, 60-90 s for B
including dgt_icc), max Rhat near 1.00, zero divergences, and that the
pilot estimates are near the targets printed at the top of each log.
Report the mean `secs` from results/study_A_pilot.csv and study_B_pilot.csv
before launching the full runs; the timing decides whether to run both
studies in one night or two.

## 2. Full runs (unattended; resume-safe)

    cd ~/dev/dgt-revision/simulations
    nohup sh -c 'Rscript study_A_lognormal.R --workers 4 > logs/A.log 2>&1; Rscript study_B_bernoulli.R --workers 4 > logs/B.log 2>&1' > logs/run.log 2>&1 &
    # A (~1.5 h) then B (~5 h) run back to back; check with: tail -3 logs/A.log logs/B.log

Progress: `tail -f logs/A.log`. Completed (cell, rep) rows are appended to
results/study_A.csv as they finish; if the Mac sleeps or the process dies,
rerun the same command and it continues from where it stopped. Keep the
laptop plugged in and set it not to sleep (System Settings > Battery >
Options > Prevent automatic sleeping on power adapter).

## 3. Analysis

    shasum -a 256 -c preregistration.md.sha256
    Rscript analyze.R
    ls results/

Outputs: results/table_study_A.tex, results/table_study_B.tex,
results/numbers.json. The pre-registered success criteria are printed as
pass/fail columns; cells that fail are reported, not dropped.

## 4. Commit

    cd ~/dev/dgt-revision && git add -A && git commit -q -m "simulations: Studies A and B complete; tables and numbers.json" && git log --oneline -1
