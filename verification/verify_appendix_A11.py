#!/usr/bin/env python3
"""verify_appendix_A11.py -- the counterexample in Appendix A11 of the
PSY-2026-0092 revision: for an absolutely continuous bimodal pair,
rho^2 > 1 - exp(-2 I(X;Y)), so no general inequality between the squared
correlation and the information ICC exists (SELF-4 correction to the
discussion under Theorem 9(c) of the submitted version).

X ~ 0.5 N(-1, d^2) + 0.5 N(1, d^2), Y = X + N(0, s^2). I(X;Y) = h(Y) - h(Y|X)
with h(Y|X) = 0.5 ln(2 pi e s^2); h(Y) by numerical integration. Deterministic.
"""
import sys
import numpy as np
from scipy import integrate
from scipy.stats import norm

def case(d, s):
    fY = lambda y: 0.5*norm.pdf(y, -1, np.hypot(d, s)) + 0.5*norm.pdf(y, 1, np.hypot(d, s))
    varX = 1 + d*d; rho2 = varX**2 / (varX * (varX + s*s))
    hY = -integrate.quad(lambda y: fY(y)*np.log(fY(y)), -12, 12, limit=400)[0]
    I = hY - 0.5*np.log(2*np.pi*np.e*s*s)
    return rho2, 1 - np.exp(-2*I)

ok = True
for d, s in [(0.05, 1.0), (0.05, 1.7), (0.3, 1.0), (0.6, 1.0)]:
    rho2, iccI = case(d, s)
    viol = rho2 > iccI + 1e-6
    ok &= viol
    print(f"{'PASS' if viol else 'FAIL'} d={d} s={s}: rho^2={rho2:.4f} > 1-exp(-2I)={iccI:.4f}")
print("The manuscript quotes the d=0.05, s=1.0 case: rho^2 = 0.501, 1 - e^{-2I} = 0.491")
sys.exit(0 if ok else 1)
