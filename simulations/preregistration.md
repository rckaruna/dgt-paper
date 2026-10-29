# Pre-registration: simulation studies for the PSY-2026-0092 revision

Frozen before any run. A sha256 of this file is stored beside it
(preregistration.md.sha256). Any change after the first run is recorded
as a dated amendment at the end of this file, never by editing the frozen
text.

Date frozen: 2026-09-25. Software: R (>= 4.5), brms, cmdstanr, posterior,
dgt 0.5.1 (github.com/rckaruna/dgt, tag v0.5.1). Hardware: one M-series
Mac, four parallel workers, one chain per core.

## Purpose

Reviewer 2 (comment 6) and the AE ask that the simulation be conducted in
the Bayesian framework the paper recommends, and that the binary case be
studied. Study A replaces the REML lognormal study of the submitted
Section 9 with a Bayesian one; the REML results move to an appendix as a
frequentist robustness check. Study B is new and tests recovery of the
Bernoulli coefficients of the revised paper's binary Proposition,
including the identity ICC_Y (absolute) = intraclass kappa.

Both studies test estimation, not the theorems: the population values are
known in closed form (Study A) or by quadrature and exact bivariate-normal
probabilities (Study B), and the question is whether the Bayesian
posterior summaries recover them with negligible bias and nominal
coverage at realistic sample sizes.

## Study A: lognormal one-facet recovery (Bayesian)

Generating model. Y_ij = exp(mu + nu_i + eps_ij), nu_i ~ N(0, s2p),
eps_ij ~ N(0, s2eps), i = 1..n_p persons, j = 1..n_j conditions, mu = 6.5.
No facet main effect (one-facet design, as in the submitted study, so
Appendix B remains comparable).

Cells. ICC_eta in {0.10, 0.25, 0.50} x s2eta in {0.25, 1.00, 2.00}, with
s2p = ICC_eta * s2eta and s2eps = s2eta - s2p. Nine cells. The grid brackets
the two empirical illustrations (lexdec s2eta about 0.06 with ICC_eta 0.40;
chess s2eta 0.56 with ICC_eta 0.17 for persons and 0.43 for items) and
extends to the s2eta = 2 regime where Solomon and Taylor's first-order
invariance fails.

Sample size. n_p = 50, n_j = 30 (1,500 observations), as in the submitted
study.

Replications. R = 200 per cell (1,800 fits). Monte Carlo standard error of
an estimated bias is SD/sqrt(200) <= 0.05/14 = 0.0035 for the largest
plausible SD, so a bias of 0.01 is detectable.

Fit. brms, family lognormal, formula y ~ 1 + (1 | person), priors
Intercept ~ normal(6.5, 5), sd ~ student_t(3, 0, 1), sigma ~
student_t(3, 0, 1); 2 chains x 1,000 iterations, 500 warmup, adapt_delta
0.9; backend cmdstanr; the model is compiled once and refitted with
update(). A fit with any Rhat > 1.05 or any divergent transition is
refitted once with adapt_delta 0.99; if it still fails it is recorded and
excluded, and the count is reported.

Estimators. From each fit, dgt_icc() posterior draws; the point estimate
is the posterior median of ICC_eta and of ICC_Y; the interval is the
central 95% credible interval.

Targets. ICC_eta = s2p / s2eta; ICC_Y = (exp(s2p) - 1) / (exp(s2eta) - 1)
(Theorem 1); O = ICC_eta / ICC_Y.

Metrics per cell. Bias (mean of posterior medians minus target), RMSE,
95% interval coverage, mean posterior interval width, and the mean
overestimation ratio, for ICC_eta and ICC_Y; the relative overestimation
(ICC_eta - ICC_Y) / ICC_Y as in the submitted Table 4. Reported with Monte
Carlo standard errors.

Success criterion, stated in advance. |bias| <= 0.01 for ICC_Y in every
cell and coverage within 0.93-0.97 are reported as "recovered with
negligible bias and nominal coverage"; a cell outside either bound is
reported as such, with the cell named, and is not dropped.

Seeds. Data seed for cell c (1..9), replicate r (1..200): 20260925 +
1000 * c + r. Sampler seed: the same integer. Recorded per row.

## Study B: Bernoulli crossed-design recovery (Bayesian)

Generating model. Y_ij | nu_i, omega_j ~ Bernoulli(h(mu + nu_i + omega_j)),
nu_i ~ N(0, s2p), omega_j ~ N(0, s2o), i = 1..n_p persons, j = 1..n_j
items; h the inverse probit or inverse logit link. Facet variance fixed at
s2o = 0.3 in every cell.

Cells. Link in {probit, logit} x ICC_eta in {0.10, 0.25, 0.50} x marginal
prevalence P in {0.50, 0.20}. Twelve cells. For each cell s2p is set from
ICC_eta = s2p / (s2p + s2o + r_g), r_g = 1 (probit) or pi^2/3 (logit),
i.e. s2p = ICC_eta (s2o + r_g) / (1 - ICC_eta), and mu is set so that the
marginal prevalence E[h(mu + nu + omega)] equals P: mu = Phi^{-1}(P)
sqrt(1 + s2p + s2o) for the probit, and by root-finding on the
Gauss-Hermite quadrature of the logit mean otherwise.

Sample size. n_p = 300, n_j = 20 (6,000 observations).

Replications. R = 100 per cell (1,200 fits).

Fit. brms, family bernoulli(link), formula y ~ 1 + (1 | person) +
(1 | item), priors Intercept ~ normal(0, 2.5), sd ~ student_t(3, 0, 1);
2 chains x 1,000 iterations, 500 warmup, adapt_delta 0.9; cmdstanr; same
refit-once rule as Study A.

Estimators. dgt_icc() with the fitted link: posterior medians and 95%
intervals of ICC_eta (link-scale, absolute), ICC_Y absolute, ICC_Y
relative, ICC_I. Also, from the raw data of each replicate, the sample
Fleiss kappa (intraclass kappa for many raters, computed from the
per-person counts) and the mean pairwise phi across item pairs.

Targets. ICC_eta by the formula above. ICC_Y absolute and relative by
dgt:::.binomial_icc_core(mu, sqrt(s2p), sqrt(s2o), 1, n_nodes = 80,
invlink) at the true parameters; for the probit, ICC_Y absolute is also
checked against the exact bivariate-normal orthant form (P11 - P^2)/(P(1-P))
and the two must agree to 1e-4. ICC_I reference value by
dgt:::.binomial_info_facet at the true parameters with K = 20,000 and
M = 2,000, repeated five times with different seeds; the mean is the
reference and its SD is reported as reference-value Monte Carlo error.

Metrics per cell. Bias, RMSE, coverage, and mean interval width for each
of the four coefficients; mean sample kappa and mean pairwise phi against
the population absolute and relative ICC_Y respectively; and, across the
two links at matched ICC_eta and P, the difference in posterior-median
ICC_Y (link-invariance check; expected near zero) against the difference
in ICC_eta (expected non-zero by construction of r_g).

Success criterion, stated in advance. |bias| <= 0.02 for ICC_Y absolute
and coverage within 0.92-0.98 (wider than Study A because the information
estimator is Monte Carlo); the identity ICC_Y (absolute) = intraclass
kappa is judged supported if the mean sample kappa lies within 0.02 of
the population ICC_Y absolute in every cell. Cells outside are reported
as such.

Seeds. Data seed for cell c (1..12), replicate r (1..100): 20260926 +
1000 * c + r; sampler seed the same. Recorded per row.

## What will be reported

Section 9 of the revised manuscript: a table for Study A (nine cells:
targets, bias, RMSE, coverage for ICC_eta and ICC_Y, mean O) and a table
for Study B (twelve cells: targets and bias, RMSE, coverage for ICC_Y
absolute, ICC_Y relative, ICC_eta, ICC_I, plus mean kappa and mean
pairwise phi), both generated by script into .tex files, with the
success criteria applied as stated. Appendix B: the submitted 15-cell
REML study unchanged, as a robustness check. Every number in the text is
read from the script output files, none typed.

## Exclusions and amendments

No cell, replicate, or coefficient is dropped after the fact. Failed
fits are counted and reported. Any deviation from this document is
recorded below with its date and reason.

Amendment 1 (2026-09-25, after the pilots, before any full run).
Sampler settings changed from 2 chains x 1,000 iterations (500 warmup) to
4 chains x 2,000 iterations (1,000 warmup) for Study A and 4 chains x
1,000 iterations (500 warmup) for Study B. Reason: in the Study A pilot
(27 fits) one fit reached Rhat 1.053 on a person intercept with 2 x 500
post-warmup draws and was excluded by the refit rule; at 1,800 fits this
would exclude on the order of sixty fits for sampling noise rather than
model failure. Pilot timing (about 3 s per Study A fit, about 30 s per
Study B fit) shows the larger runs are affordable. Nothing else changes:
cells, replications, priors, seeds, targets, metrics and success criteria
are as frozen. The refit rule is unchanged. sha256 of this file after the
amendment is recorded as a second line of preregistration.md.sha256.
