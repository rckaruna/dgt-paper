# Re-analysis of Huynh (1977) with dgt 0.5.0

The check called for in Section 6f of the revision plan: re-fit the 12-person by
6-item dichotomous data set that Choi (2013) prints as Figure 1 on p. 25, and see
whether the response-scale D-study coefficients land near her Method C.

Run 11 September 2026. `packageVersion("dgt")` is **0.5.0**; the check is
asserted at the top of `huynh_reanalysis.R` and printed twice in
`console-log.txt` (`CHECK: packageVersion("dgt") == 0.5.0 -- OK`, and the last
line). R 4.6.0, brms 2.23.0, cmdstanr 0.8.0, CmdStan 2.39.0, posterior 1.7.0.
No file in the `dgt` package and no `.tex` file was modified.

## Data

Transcribed from a 300 dpi rendering of printed p. 25 of `reading/Choi2013.pdf`,
with the figure block re-rendered at 600 dpi for cell-by-cell reading. The
procedure, the transcribed matrix and an independent cross-check against the
PDF's text layer are in `TRANSCRIPTION.md`; the matrix is recorded in
`huynh1977.csv`.

**No cell was illegible.** All 72 cells read unambiguously at both resolutions,
and the two independent readings (glyph rendering, embedded text layer) agree
cell for cell. 13 cells are `NA`, leaving **59 observed responses**; responses
per person range from 4 to 6, per item from 8 to 11. Overall proportion correct
0.627.

## Model

```r
brm(y ~ 1 + (1 | person) + (1 | item),
    family = bernoulli(link = "logit"),
    prior  = c(prior(normal(0, 2),         class = "Intercept"),
               prior(student_t(3, 0, 2.5), class = "sd", group = "person"),
               prior(student_t(3, 0, 2.5), class = "sd", group = "item")),
    backend = "cmdstanr", chains = 4, iter = 3000, warmup = 1000,
    control = list(adapt_delta = 0.99, max_treedepth = 12))
```

8000 post-warmup draws; max Rhat 1.0022, min bulk ESS 1461, 0 divergent
transitions. With one binary response per cell the person-by-item interaction is
not separable from binomial error, which is why Choi's Method D reports
person x item = 0.000 and error 3.287 (= pi^2/3); our model is the same
specification.

Variance components on the logit scale, posterior median (95% CrI):

| | ours | Choi Method D |
|---|---|---|
| sigma^2 person | 1.823 (0.056, 9.232) | 0.978 |
| sigma^2 item | 0.709 (0.004, 6.582) | 0.246 |
| sigma^2 error | pi^2/3 = 3.290 (fixed) | 3.287 |

Her ML point estimates sit inside our intervals, which span two orders of
magnitude. Twelve persons and six items do not identify a variance component.

## Single-item coefficients, `dgt_icc()`

| measure | posterior median (95% CrI) |
|---|---|
| `ICC_eta` (link scale, pi^2/3 convention) | **0.295 (0.013, 0.674)** |
| `ICC_Y` response scale, **absolute** | **0.200 (0.009, 0.480)** |
| `ICC_Y` response scale, **relative** | **0.225 (0.011, 0.503)** |
| `ICC_I` (information) | 0.193 (0.010, 0.435) |

`ICC_eta` uses the pi^2/3 residual convention: differencing the dgt draws
against `sd_p^2 / (sd_p^2 + sd_i^2 + pi^2/3)` gives `max |difference| = 0`.

## Which of our coefficients correspond to which of hers

Choi's Table 4 note c says the coefficients were computed at I' = 6 items. Her
formulas, recovered from her own Method D column and confirmed arithmetically,
are the standard p x i D-study pair:

* `E(rho^2)` (generalizability, **relative**) `= s_p^2 / (s_p^2 + (s_pi^2 + s_e^2)/n)`.
  Method D: `0.978/(0.978 + 3.287/6) = 0.641`. Method C:
  `0.052/(0.052 + (0.248+0.187)/6) = 0.418`.
* `Phi` (dependability, **absolute**) `= s_p^2 / (s_p^2 + (s_i^2 + s_pi^2 + s_e^2)/n)`.
  Method D: `0.978/(0.978 + (0.246+3.287)/6) = 0.624`. Method C:
  `0.052/(0.052 + (0.057+0.248+0.187)/6) = 0.388`.

`dgt_dstudy()` returns one curve per scale, and both are the **absolute** form:
its link-scale curve is `s_p^2 / (s_p^2 + (s_i^2 + pi^2/3)/n)`, which is exactly
Choi's `Phi`. So:

| Choi | our coefficient | source |
|---|---|---|
| `Phi`, Method D (latent) = 0.624 | `dgt_dstudy()` **link-scale** at n = 6 | dgt 0.5.0 |
| `Phi`, Method C (response) = 0.388 | `dgt_dstudy()` **response-scale** at n = 6 | dgt 0.5.0 |
| `E(rho^2)`, Method D (latent) = 0.641 | link-scale **relative** at n = 6 | computed by hand from the same posterior; `dgt_dstudy()` does not return a relative curve |
| `E(rho^2)`, Method C (response) = 0.418 | response-scale **relative** at n = 6 | computed by hand, same Gauss-Hermite nodes dgt uses, item main effect dropped from the error term |

Two checks that the hand computations are the same machinery: the hand link-scale
absolute at n = 6 reproduces the `dgt_dstudy()` median to four decimals (0.7150
against 0.7150), and the hand response-scale absolute reproduces it to within
Monte Carlo noise (0.600 against 0.601).

## Results at 6 items

Posterior median (95% credible interval).

| coefficient | ours at 6 items | Choi Method D (latent) | Choi Method C (response) |
|---|---|---|---|
| `E(rho^2)`, latent / logit scale | **0.769 (0.093, 0.944)** | 0.641 | — |
| `Phi`, latent / logit scale | **0.715 (0.076, 0.925)** | 0.624 | — |
| `E(rho^2)`, response scale | **0.635 (0.064, 0.859)** | — | 0.418 |
| `Phi`, response scale | **0.601 (0.054, 0.847)** | — | 0.388 |

## Do the response-scale coefficients land near Method C?

**No.** Our response-scale coefficients at 6 items are 0.635 and 0.601. Choi's
Method C values are 0.418 and 0.388. Our response-scale numbers sit roughly 0.21
above hers on both coefficients, and land instead within 0.01-0.02 of her
*Method D* values (0.641 and 0.624) — the opposite of what the conjecture in
plan 6f would need.

That is the whole finding, and nothing beyond it should be claimed. The data set
is 59 observed binary responses, and the credible intervals are correspondingly
useless for discrimination: all four of Choi's numbers, Method C and Method D
alike, fall inside our corresponding 95% intervals. The re-analysis therefore
does not *exclude* Method C either. It simply gives no support to the reading
that Method C's low coefficients are the response-scale quantity.

A second observation, independent of the interval width, points the same way.
Our scale gap closes as items accumulate, because both coefficients tend to 1:

| n items | link-scale | response-scale | gap | ratio |
|---|---|---|---|---|
| 1 | 0.295 | 0.201 | 0.094 | 0.681 |
| 6 | 0.715 | 0.601 | 0.114 | 0.841 |
| 12 | 0.834 | 0.751 | 0.083 | 0.901 |
| 25 | 0.913 | 0.863 | 0.050 | 0.945 |
| 50 | 0.954 | 0.926 | 0.028 | 0.971 |

Choi reports (pp. 27-28) that the gap between Method C and the other methods
"does not disappear as the D study item sample size increases to many, such as 25
items". A pure scale difference does disappear: by 25 items ours is down to 0.050
and the ratio is 0.945. A persistent gap is not what the scale difference
predicts.

## Conclusion for the revision

The conditional in plan 6f — "if the response-scale coefficient lands near
Method C, this is a published, unexplained discrepancy the Proposition explains;
if not, Choi's artefact explanation stands and nothing is said" — resolves to
**"if not"**. Choi's unbalanced-design explanation for Method C stands, and the
manuscript should say nothing about her Table 4. The strict ordering
`ICC_Y < ICC_eta` does hold here, with posterior probability 1 at every item
count examined, but that ordering is a property of the Proposition, not evidence
about Method C.

## Files

* `huynh1977.csv` — the transcribed data, wide form, `NA` for missing.
* `TRANSCRIPTION.md` — rendering procedure, legibility note, cross-check.
* `huynh_reanalysis.R` — the script.
* `console-log.txt` — the full console log, including the dgt version check.
* `results_huynh.rds` — the `dgt_icc` and `dgt_dstudy` objects and the per-draw
  coefficients (git-ignored).
