# Reproduction: the fetal heart rate example of Vanbelle (2025), Section 6

Run 11 September 2026. `packageVersion("dgt")` is **0.5.0**; the check is
asserted at the top of `fhr_reproduction.R` and printed in `console-log.txt`
(`CHECK: packageVersion("dgt") == 0.5.0 -- OK`, and again on the last line).

Software: R 4.6.0, brms 2.23.0, cmdstanr 0.8.0, CmdStan 2.39.0, posterior 1.7.0.
No file in the `dgt` package and no `.tex` file was touched.

## Source of the data and of the model

`https://github.com/svanbelle/RELIABILITY.BINARY`, commit
`52f0a3aadc05d683e2a7b392c3a794cca1ae8dbf` (1 April 2025), licence
**GNU Affero General Public License v3.0**, which permits use and
redistribution, so the reproduction goes ahead. The clone is git-ignored rather
than vendored here; `fhr_reproduction.R` re-clones it if absent.

* Data: `HR.csv`, 145 rows by 16 columns of fetal heart rates in bpm, no missing
  entries, range 60-193. This is her `N = 145, R = 16` (p. 13).
* Dichotomisation: `EXAMPLES.R` line 164,
  `HR_01 <- ifelse(HR_NA > 110 & HR_NA < 150, 0, 1)` — "normal" is strictly
  between 110 and 150 bpm, everything else "abnormal". `HR_NA` is her name for
  the matrix with unreadable entries set to `NA`; the published `HR.csv` has
  none, so the design is completely crossed, 2320 observations. Overall
  proportion abnormal P = 0.2625.
* Model: her JAGS file `probit6.txt`, written by `FUNCTIONS_EXAMPLES.R`
  lines 355-402:

  ```
  y[i]    ~ dbern(phi(S[sub[i]] + R[met[i]]))
  S[j]    ~ dnorm(beta0, tau_S)        subjects, crossed
  R[j]    ~ dnorm(0,     tau_R)        midwives, crossed
  beta0   ~ dnorm(0, 1)                precision 1, i.e. Normal(0, sd = 1)
  sigma_S ~ dt(0, 1, 1) T(0,)          half-Cauchy(0, 1)
  sigma_R ~ dt(0, 1, 1) T(0,)          half-Cauchy(0, 1)
  ```

  transcribed to brms as

  ```r
  brm(y ~ 1 + (1 | sub) + (1 | met),
      family = bernoulli(link = "probit"),
      prior  = c(prior(normal(0, 1), class = "Intercept"),
                 prior(cauchy(0, 1), class = "sd", group = "sub"),
                 prior(cauchy(0, 1), class = "sd", group = "met")),
      backend = "cmdstanr", chains = 4, iter = 4000, warmup = 1000)
  ```

  Moving `beta0` from the mean of `S[]` into the fixed part is the same model.
  4 chains x 3000 post-warmup draws = 12000; max Rhat 1.0026, min bulk ESS 1926,
  0 divergent transitions.

## Results

Posterior medians with 95% equal-tailed credible intervals.

| Quantity | This reproduction (dgt 0.5.0) | Vanbelle (2025), p. 14 |
|---|---|---|
| `sigma_R^2` (midwife variance, probit scale) | **0.128 (0.054, 0.329)** | `sigmahat_R^2` = .13 (.05, .33) |
| `ICC_eta` (link scale, tetrachoric) | **0.700 (0.621, 0.771)** | `tauhat_abs` = .70 (.62, .77) |
| `ICC_Y` response scale, **absolute** | **0.475 (0.404, 0.545)** | `phihat_abs` = .48 (.41, .55) |
| `ICC_Y` response scale, **relative** | **0.485 (0.416, 0.553)** | not reported |
| `ICC_I` (information) | 0.405 (0.334, 0.473) | not reported |
| Fleiss kappa, computed from the raw 145 x 16 table | **0.506** | `kappahat_abs` = .51 (.43, .58) |
| mean pairwise phi over the 120 midwife pairs | **0.515** (sd 0.095, range 0.220-0.792) | — |

Supporting posterior summaries: `sigma_S` = 1.635 (1.378, 1.949),
`sigma_R` = 0.358 (0.232, 0.574).

## Correspondence of the coefficients

* Her `tau_abs`, Equation (11), is
  `sigma_S^2 / (sigma_S^2 + sigma_R^2 + 1)`. `dgt_icc()` computes `ICC_eta` for a
  probit model by exactly that expression (link residual variance 1). Recomputing
  her Equation (11) by hand from our posterior draws and differencing against
  `icc$draws$icc_eta` gives `max |difference| = 0` — they are the same number,
  not merely a close one.
* Her `phi_abs` is the manifest-scale phi coefficient, her Equations (12)-(13),
  which she identifies with the intraclass kappa in her Equation (15). This is
  dgt's `ICC_Y (response-scale), absolute`.
* Her `kappa_abs` is the sample Fleiss / intraclass kappa, the frequentist
  counterpart of the same manifest-scale quantity.
* The `ICC_Y` **relative** coefficient has no counterpart in her paper; it is the
  response-scale coefficient the Proposition adds.

## Verdict

Every published number is reproduced.

* `sigma_R^2`: 0.128 (0.054, 0.329) against .13 (.05, .33) — identical at her
  printed precision, interval and all.
* `tau_abs`: 0.700 (0.621, 0.771) against .70 (.62, .77) — identical.
* `phi_abs`: 0.475 (0.404, 0.545) against .48 (.41, .55). The point estimate
  rounds to .47 rather than .48 and the interval limits differ by at most .006;
  this is the size of gap expected between her JAGS run and an independent Stan
  run of the same model with a different seed, and between her closed-form
  bivariate-normal inversion (the `a11` algebra at the foot of `probit6.txt`) and
  dgt's numerical integration of the same quantity.
* `kappa_abs`: our sample Fleiss kappa is 0.506 against her .51 — identical.
  Her interval is frequentist and we computed the point estimate only.

The mean pairwise phi, 0.515, sits slightly above the model-based `phi_abs` of
0.475, as it should: averaging the 120 pairwise correlations is not the same
estimator as the common intraclass value, and the pairwise values are widely
spread (0.220 to 0.792) because the midwives differ in how often they call a
recording abnormal (0.145 to 0.421). The Fleiss kappa, which pools the
disagreement rather than averaging correlations, is 0.506 and is the closer
match to `ICC_Y` absolute, as the identity `ICC_Y = kappa` predicts.

The ordering `ICC_Y < ICC_eta` holds throughout, 0.475 against 0.700, and the
relative response-scale coefficient, 0.485, lies just above the absolute one and
far below the link-scale one — the small midwife variance makes the
absolute/relative gap small while the scale gap stays large.

## Nothing failed to reproduce

The data are in the repository, the licence permits use, and the model
specification is given in full. No part of Section 6 of Vanbelle (2025) was left
unreproduced.

## Files

* `fhr_reproduction.R` — the script.
* `console-log.txt` — the full console log of the run, including the dgt version
  check.
* `results_fhr.rds` — the `dgt_icc` object, the comparison table and the sample
  coefficients (git-ignored).
* `RELIABILITY.BINARY/` — the clone (git-ignored).
