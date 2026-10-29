#!/usr/bin/env python3
"""verify_binary_proposition.py -- numerical verification of the Bernoulli
section of the PSY-2026-0092 revision (Proposition: binary response-scale
ICC; Remark: small effects and the information ICC).

Checks, each printed PASS/FAIL:
  A1  absolute ICC_Y (quadrature) == (P11 - P^2)/(P(1-P)) == Fleiss kappa == phi
  A2  relative ICC_Y == Var m / E_omega[P_omega(1-P_omega)], rel >= abs, equal iff sd_o = 0
  C1  probit: P = Phi(z), P11 = Phi_2(z,z;rho), rho = sd_p^2/V
  C2  probit: Mehler/tetrachoric series matches exact P11 - P^2 for rho < 1
  C3  probit: sum_k a_k -> P(1-P) as rho -> 1 (Abel), and ICC_Y < rho
  D1  small-effect limit of ICC_Y/ICC_eta matches the closed form (probit, logit)
  D2  probit limit <= 2/pi; logit (no facet) limit <= pi^2/12
  R1  ICC_I / ICC_Y -> 1 as sd_p -> 0
  R2  ICC_I <= ICC_Y on a 504-cell grid over both links (empirical regularity)
  R3  ICC_Y < ICC_eta on the same grid, both links (ordering (b); proved for probit only)

Requires numpy, scipy. Deterministic (Gauss-Hermite quadrature); no seeds.
"""
import itertools, sys
import numpy as np
from scipy.stats import norm, multivariate_normal
from scipy.special import expit

X, W = np.polynomial.hermite_e.hermegauss(120); W = W / np.sqrt(2 * np.pi)
results = []

def check(name, cond, detail=""):
    results.append(cond)
    print(f"{'PASS' if cond else 'FAIL'} {name}  {detail}")

def moments(mu, sp, so, inv):
    eta = np.add.outer(mu + sp * X, so * X)
    pig = inv(eta)
    m_nu = pig @ W; P_om = W @ pig; P = W @ m_nu
    var_obj = W @ (m_nu - P) ** 2
    H = lambda p: -(p * np.log(p) + (1 - p) * np.log(1 - p))
    I = H(P) - W @ H(np.clip(m_nu, 1e-12, 1 - 1e-12))
    return dict(P=P, var_obj=var_obj,
                icc_abs=var_obj / (P * (1 - P)),
                icc_rel=var_obj / (W @ (P_om * (1 - P_om))),
                icc_I=1 - np.exp(-2 * I))

def probit_exact(mu, sp, so):
    V = 1 + sp ** 2 + so ** 2; z = mu / np.sqrt(V); rho = sp ** 2 / V
    P = norm.cdf(z)
    P11 = multivariate_normal(mean=[-mu, -mu], cov=[[V, sp ** 2], [sp ** 2, V]]).cdf([0, 0])
    return z, rho, P, P11

def a_terms(z, K):
    psi = np.empty(K); psi[0] = 1.0; psi[1] = z
    for n in range(1, K - 1):
        psi[n + 1] = (z * psi[n] - np.sqrt(n) * psi[n - 1]) / np.sqrt(n + 1)
    return norm.pdf(z) ** 2 * psi ** 2 / np.arange(1, K + 1)

# ---- A1, A2, C1 ----
for mu, sp, so in [(-1.19, 1.66, 0.37), (0.5, 1.0, 1.0), (-2.0, 3.0, 0.2), (0.0, 0.7, 0.0)]:
    z, rho, P, P11 = probit_exact(mu, sp, so)
    q = moments(mu, sp, so, norm.cdf)
    p_o = 2 * P11 + 1 - 2 * P; p_e = P ** 2 + (1 - P) ** 2
    kappa = (p_o - p_e) / (1 - p_e); phi = (P11 - P ** 2) / (P * (1 - P))
    check(f"A1 kappa=phi=ICC_Y abs (mu={mu},sp={sp},so={so})",
          abs(kappa - phi) < 1e-12 and abs(phi - q["icc_abs"]) < 2e-4, f"{kappa:.5f} {phi:.5f} {q['icc_abs']:.5f}")
    check(f"A2 rel>=abs, equality iff so=0 (so={so})",
          (q["icc_rel"] > q["icc_abs"] + 1e-9) if so > 0 else abs(q["icc_rel"] - q["icc_abs"]) < 1e-12,
          f"rel={q['icc_rel']:.5f} abs={q['icc_abs']:.5f}")
    check(f"C1 probit P=Phi(z) (mu={mu},sp={sp},so={so})", abs(q["P"] - P) < 1e-6, f"{q['P']:.6f} vs {P:.6f}")

# ---- C2, C3 ----
for z in [-0.603, 0.289, 0.0, 1.2]:
    a = a_terms(z, 400000); P = norm.cdf(z)
    ok = True
    for rho in [0.3, 0.7079, 0.99]:
        S = np.sum(rho ** np.arange(1, len(a) + 1) * a)
        exact = multivariate_normal(mean=[-z, -z], cov=[[1, rho], [rho, 1]]).cdf([0, 0]) - P ** 2
        ok &= abs(S - exact) < 2e-6
    check(f"C2 tetrachoric series == exact P11-P^2 (z={z})", ok)
    # Abel limit: partial sums approach P(1-P) from below, slowly (k^-1/2 tail); accept 0.2% at K=4e5
    check(f"C3 sum a_k -> P(1-P) (z={z})", 0 < P * (1 - P) - a.sum() < 0.002 * P * (1 - P),
          f"sum={a.sum():.6f} target={P*(1-P):.6f}; a_2>0 or a_3>0: {a[1] > 0 or a[2] > 0}")
    for rho in [0.2, 0.7079]:
        icc = np.sum(rho ** np.arange(1, len(a) + 1) * a) / (P * (1 - P))
        check(f"C3 ICC_Y < rho (z={z}, rho={rho})", icc < rho, f"{icc:.4f} < {rho}")

# ---- D1, D2 ----
for inv, rvar, name in [(norm.cdf, 1.0, "probit"), (expit, np.pi ** 2 / 3, "logit")]:
    for mu in [0.0, -1.0, 1.5]:
        so = 0.6
        sp = 0.01
        q = moments(mu, sp, so, inv); ratio = q["icc_abs"] / (sp ** 2 / (sp ** 2 + so ** 2 + rvar))
        if name == "probit":
            z0 = mu / np.sqrt(1 + so ** 2); pred = norm.pdf(z0) ** 2 / (norm.cdf(z0) * (1 - norm.cdf(z0)))
        else:
            d = W @ (expit(mu + so * X) * (1 - expit(mu + so * X))); P = W @ expit(mu + so * X)
            pred = d ** 2 * (so ** 2 + rvar) / (P * (1 - P))
        check(f"D1 small-effect limit ({name}, mu={mu})", abs(ratio - pred) < 2e-3, f"ratio={ratio:.4f} pred={pred:.4f}")
        if name == "probit":
            check(f"D2 probit limit <= 2/pi (mu={mu})", pred <= 2 / np.pi + 1e-12, f"{pred:.4f}")
# logit bound with no facet
for mu in [0.0, -1.0, 1.5, 3.0]:
    P = expit(mu); lim = P * (1 - P) * np.pi ** 2 / 3
    check(f"D2 logit no-facet limit <= pi^2/12 (mu={mu})", lim <= np.pi ** 2 / 12 + 1e-12, f"{lim:.4f}")

# ---- R1, R2 ----
for inv, name in [(norm.cdf, "probit"), (expit, "logit")]:
    for mu in [0.0, -1.0]:
        q = moments(mu, 0.03, 0.6, inv)
        check(f"R1 ICC_I/ICC_Y -> 1 ({name}, mu={mu}, sp=0.03)", abs(q["icc_I"] / q["icc_abs"] - 1) < 2e-3,
              f"{q['icc_I']/q['icc_abs']:.5f}")
viol = 0; tot = 0
for inv, name in [(norm.cdf, "probit"), (expit, "logit")]:
    for mu, sp, so in itertools.product([-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3], [0.05, 0.2, 0.5, 1, 2, 3, 5], [0, 0.3, 1, 2]):
        q = moments(mu, sp, so, inv); tot += 1
        if q["icc_I"] > q["icc_abs"] + 1e-9: viol += 1
check(f"R2 ICC_I <= ICC_Y on {tot}-cell grid (empirical regularity)", viol == 0, f"violations={viol}")
# ---- R3: ICC_Y (absolute) < ICC_eta on the same grid, both links (ordering (b) for the logit is not proved) ----
viol3 = 0; tot3 = 0; min_ratio = 1.0
for inv, rvar, name in [(norm.cdf, 1.0, "probit"), (expit, np.pi ** 2 / 3, "logit")]:
    for mu, sp, so in itertools.product([-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3], [0.05, 0.2, 0.5, 1, 2, 3, 5], [0, 0.3, 1, 2]):
        q = moments(mu, sp, so, inv); tot3 += 1
        icc_eta = sp ** 2 / (sp ** 2 + so ** 2 + rvar)
        if q["icc_abs"] >= icc_eta: viol3 += 1
        min_ratio = min(min_ratio, icc_eta / q["icc_abs"] if q["icc_abs"] > 0 else 1.0)
check(f"R3 ICC_Y < ICC_eta on {tot3}-cell grid, both links (ordering (b))", viol3 == 0, f"violations={viol3}, min ICC_eta/ICC_Y={min_ratio:.4f}")

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
