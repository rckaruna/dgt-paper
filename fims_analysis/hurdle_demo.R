# hurdle_demo.R -- Section 9 hurdle demonstration (Theorem 4) with known truth
# Scenario: daily minutes spent on a practice task, recorded by ecological
# momentary assessment; a day with no practice is a structural zero.
# Generating model (all parameters stated; Reviewer 1, C3):
#   n_p = 80 persons, n_j = 60 days
#   engagement: Z_ij ~ Bernoulli(pi_i), logit(pi_i) = alpha + u_i, u_i ~ N(0, tau^2)
#   intensity:  log Y*_ij = mu + nu_i + eps_ij, nu_i ~ N(0, s2p), eps_ij ~ N(0, s2e)
#   u_i and nu_i independent; Y_ij = Z_ij * Y*_ij
ALPHA <- 0.70; TAU <- 0.80; MU <- 3.50; SD_P <- 0.40; SD_E <- 0.60
N_P <- 80L; N_J <- 60L; SEED <- 20260926L
suppressPackageStartupMessages({ library(brms); library(dgt); library(posterior); library(ggplot2) })
stopifnot(packageVersion("dgt") >= "0.5.1")
set.seed(SEED)
u  <- rnorm(N_P, 0, TAU); nu <- rnorm(N_P, 0, SD_P)
dat <- expand.grid(person = factor(seq_len(N_P)), day = seq_len(N_J))
pi_i <- plogis(ALPHA + u); z <- rbinom(nrow(dat), 1, pi_i[dat$person])
ystar <- exp(MU + nu[dat$person] + rnorm(nrow(dat), 0, SD_E))
dat$y <- z * ystar
cat(sprintf("Simulated %d x %d; zero proportion %.3f; median positive minutes %.1f\n",
            N_P, N_J, mean(dat$y == 0), median(dat$y[dat$y > 0])))
saveRDS(dat, "hurdle_demo_data.rds")

fit <- brm(bf(y ~ 1 + (1 | person), hu ~ 1 + (1 | person)), data = dat,
           family = hurdle_lognormal(),
           prior = c(prior(normal(3.5, 5), class = Intercept), prior(normal(0, 2.5), class = Intercept, dpar = "hu"),
                     prior(student_t(3, 0, 1), class = sd), prior(student_t(3, 0, 1), class = sd, dpar = "hu"),
                     prior(student_t(3, 0, 1), class = sigma)),
           chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = SEED, backend = "cmdstanr",
           refresh = 0, file = "hurdle_demo_fit", file_refit = "never")
cat(sprintf("max Rhat %.3f, divergences %d\n", max(brms::rhat(fit), na.rm = TRUE), sum(nuts_params(fit, pars = "divergent__")$Value)))

set.seed(1); icc <- dgt_icc(fit, person_group = "person");   print(icc)
set.seed(1); var <- dgt_variance(fit, person_group = "person"); print(var)
set.seed(1); inf <- dgt_info_icc(fit, person_group = "person"); print(inf)

# ---- closed-form five components at given parameters (for the counterfactual curves) ----
# Var(pi) and E[pi(1-pi)] by Monte Carlo over u; lognormal moments in closed form.
comp5 <- function(alpha, tau, mu, sp, se, M = 200000, seed = 1) {
  set.seed(seed); p <- plogis(alpha + rnorm(M, 0, tau))
  Epi <- mean(p); Vpi <- var(p); Epq <- mean(p * (1 - p))
  Emu  <- exp(mu + sp^2 / 2 + se^2 / 2)                       # E[mu*_i]
  Vmu  <- exp(2 * mu + se^2) * exp(sp^2) * (exp(sp^2) - 1)    # Var(mu*_i)
  Emu2 <- Vmu + Emu^2
  EVarY <- exp(2 * mu + 2 * sp^2 + se^2) * (exp(se^2) - 1)     # E[Var(Y*|nu)] = E[mu*^2 (e^{se^2}-1)]
  V <- c(V1 = Epq * Emu2, V2 = Epi * EVarY, V3 = Vpi * Emu^2, V4 = Epi^2 * Vmu, V5 = Vpi * Vmu)
  list(V = V, icc_comp = sum(V[3:5]) / sum(V), icc_Z = Vpi / (Vpi + Epq),
       icc_Ystar = (exp(sp^2) - 1) / (exp(sp^2 + se^2) - 1))
}
truth <- comp5(ALPHA, TAU, MU, SD_P, SD_E)
cat("Closed-form components at the true parameters:\n"); print(round(truth$V / sum(truth$V), 3))
cat(sprintf("ICC_comp %.3f, ICC_Z %.3f, ICC_Y* %.3f (truth)\n", truth$icc_comp, truth$icc_Z, truth$icc_Ystar))
# cross-check the closed form against dgt at posterior means
post <- as_draws_df(fit)
pm <- c(alpha = mean(post$b_hu_Intercept), tau = mean(post$sd_person__hu_Intercept), mu = mean(post$b_Intercept),
        sp = mean(post$sd_person__Intercept), se = mean(post$sigma))
# brms hu models P(Y = 0): engagement pi = 1 - plogis(hu), i.e. alpha_engage = -b_hu
chk <- comp5(-pm["alpha"], pm["tau"], pm["mu"], pm["sp"], pm["se"])
cat(sprintf("Closed form at posterior means: ICC_comp %.3f (dgt reports %s)\n", chk$icc_comp,
            format(round(icc$summary$estimate[grepl("comp", icc$summary$measure)][1], 3))))

# ---- D-study figure: more occasions vs halved intensity noise (Reviewer 1, m6) ----
sb <- function(r, n) n * r / (1 + (n - 1) * r)
n_grid <- 1:60
base <- comp5(-pm["alpha"], pm["tau"], pm["mu"], pm["sp"], pm["se"])$icc_comp
half <- comp5(-pm["alpha"], pm["tau"], pm["mu"], pm["sp"], pm["se"] / 2)$icc_comp
df <- rbind(data.frame(n = n_grid, r = sb(base, n_grid), scenario = "Current intensity noise"),
            data.frame(n = n_grid, r = sb(half, n_grid), scenario = "Intensity noise halved (sigma_e / 2)"))
g <- ggplot(df, aes(n, r, colour = scenario)) + geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = 3) + scale_colour_manual(values = c("#2166ac", "#b2182b")) +
  labs(x = "Number of occasions (n)", y = expression(E*rho^2*" (composite)"), colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figure_hurdle_dstudy.pdf", g, width = 6.5, height = 4)
n80 <- function(r) { v <- which(sb(r, n_grid) >= 0.8); if (length(v)) n_grid[v[1]] else NA }
cat(sprintf("Composite ICC %.3f -> n*(.80) = %s occasions; with sigma_e halved, ICC %.3f -> n*(.80) = %s\n",
            base, n80(base), half, n80(half)))
saveRDS(list(params = c(ALPHA = ALPHA, TAU = TAU, MU = MU, SD_P = SD_P, SD_E = SD_E, N_P = N_P, N_J = N_J, SEED = SEED),
             zero_prop = mean(dat$y == 0), icc = icc, var = var, inf = inf, truth = truth, post_means = pm,
             base_icc = base, half_icc = half, n80_base = n80(base), n80_half = n80(half)), "hurdle_demo_results.rds")
cat("Wrote hurdle_demo_results.rds and figure_hurdle_dstudy.pdf\n")
