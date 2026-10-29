#!/usr/bin/env python3
"""Table 1 of the PSY-2026-0092 revision: lognormal response-scale ICC on a
factorial grid (Reviewer 1 C1, Reviewer 2 m2, Reviewer 3 m6).

Closed forms only (Theorem 1, Definition of O, Theorem 3 inverted):
  ICC_Y  = (exp(s2p) - 1) / (exp(s2eta) - 1),  s2p = ICC_eta * s2eta
  O      = ICC_eta / ICC_Y
  n*(c)  = ceil( c (1 - ICC) / (ICC (1 - c)) ),  c = 0.80

Writes manuscript/generated/table_lognormal_grid.tex. Nothing in the table is
typed by hand. Run from anywhere:  python3 verification/table1_lognormal_grid.py
"""
import math
import os

ICC_ETA = [0.10, 0.25, 0.50, 0.75]
S2ETA = [0.1, 0.5, 1.0, 2.0, 4.0]
C = 0.80


def icc_y(icc_eta, s2eta):
    s2p = icc_eta * s2eta
    return math.expm1(s2p) / math.expm1(s2eta)


def n_star(icc, c=C):
    return math.ceil(c * (1.0 - icc) / (icc * (1.0 - c)) - 1e-12)


def fmt3(x):
    s = f"{x:.3f}"
    return s[1:] if s.startswith("0.") else s


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "..", "manuscript", "generated", "table_lognormal_grid.tex")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    rows = []
    rows.append(r"\begin{table}[t]")
    rows.append(r"\centering")
    rows.append(r"\caption{Response-Scale ICC for Lognormal Measurements on a Factorial Grid of Link-Scale ICC and Total Log-Scale Variance}")
    rows.append(r"\label{tab:lognormal-illustration}")
    rows.append(r"\begin{tabular}{l" + "c" * len(S2ETA) + "}")
    rows.append(r"\toprule")
    rows.append(r" & \multicolumn{" + str(len(S2ETA)) + r"}{c}{$\sigma^2_\eta$} \\")
    rows.append(r"\cmidrule(lr){2-" + str(len(S2ETA) + 1) + "}")
    rows.append(r"$\icc_\eta$ & " + " & ".join(f"{s:g}" for s in S2ETA) + r" \\")
    rows.append(r"\midrule")
    rows.append(r"\multicolumn{" + str(len(S2ETA) + 1) + r"}{l}{\textit{Response-scale ICC $\icc_Y^{\mathrm{LN}}$ (overestimation ratio $\mathcal{O}$)}} \\")
    for r in ICC_ETA:
        cells = []
        for s in S2ETA:
            y = icc_y(r, s)
            cells.append(f"{fmt3(y)} ({r / y:.2f})")
        rows.append(f"{fmt3(r)} & " + " & ".join(cells) + r" \\")
    rows.append(r"\addlinespace")
    rows.append(r"\multicolumn{" + str(len(S2ETA) + 1) + r"}{l}{\textit{Occasions required for $E\rho^2 \ge .80$: response scale $n^*_Y$ (link scale $n^*_\eta$)}} \\")
    for r in ICC_ETA:
        cells = []
        ne = n_star(r)
        for s in S2ETA:
            y = icc_y(r, s)
            cells.append(f"{n_star(y)} ({ne})")
        rows.append(f"{fmt3(r)} & " + " & ".join(cells) + r" \\")
    rows.append(r"\bottomrule")
    rows.append(r"\end{tabular}")
    rows.append(r"\begin{flushleft}")
    rows.append(r"\small\textit{Note.} Each cell is computed from Theorem~\ref{thm:lognormal-icc} with $\sigma^2_p = \icc_\eta\,\sigma^2_\eta$ and $\sigma^2_e = (1 - \icc_\eta)\,\sigma^2_\eta$. The overestimation ratio is $\mathcal{O} = \icc_\eta/\icc_Y^{\mathrm{LN}}$ (Definition~\ref{def:ratios}). $n^*$ is the smallest number of occasions for which the Spearman--Brown reliability of the mean reaches $.80$ (Theorem~\ref{thm:dstudy-arith}); the classical D-study based on $\icc_\eta$ reports the value in parentheses. The two empirical illustrations have $\sigma^2_\eta \approx 0.06$ (Section~\ref{sec:rt}) and $\sigma^2_\eta = 0.56$ (Section~\ref{sec:chess}).")
    rows.append(r"\end{flushleft}")
    rows.append(r"\end{table}")
    with open(out, "w") as f:
        f.write("\n".join(rows) + "\n")
    print("wrote", os.path.normpath(out))
    # echo the grid for the log
    for r in ICC_ETA:
        print(fmt3(r), [f"{fmt3(icc_y(r, s))}/{r / icc_y(r, s):.2f}/{n_star(icc_y(r, s))}" for s in S2ETA], "n*_eta", n_star(r))


if __name__ == "__main__":
    main()
