#!/usr/bin/env python3
"""Tables for the Supplementary Material of PSY-2026-0092 (revision).

Reads only files written by the analysis scripts and writes LaTeX tables
into manuscript/generated/. Nothing here is typed by hand.

  simulations/results/numbers.json            -> table_study_A_full.tex, table_study_B_full.tex
  reviewer2/vanbelle_design_results.csv       -> table_r2_design.tex
  reproduction/vanbelle-fhr/console-log.txt   -> table_vanbelle_fhr.tex

Run from anywhere:  python3 verification/supplement_tables.py
"""
import csv
import json
import os
import re
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "manuscript", "generated")
os.makedirs(OUT, exist_ok=True)


def f3(x):
    return f"{round(x, 3) + 0:.3f}"


def f2(x):
    return f"{round(x, 2) + 0:.2f}"


def f4(x):
    return f"{round(x, 4) + 0:.4f}"


def write(name, lines):
    path = os.path.join(OUT, name)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", os.path.relpath(path, ROOT))


# ---------------------------------------------------------------- Study A
num = json.load(open(os.path.join(ROOT, "simulations", "results", "numbers.json")))
A = num["study_A"]
lines = [
    r"\begin{table}[tp]",
    r"\centering",
    r"\caption{Study A in full: every per-cell metric of the pre-registration ($n_p = 50$, $n_j = 30$, $R = 200$ per cell; one fit excluded by the refit rule). Targets are the analytic values of Theorem~1 of the main text. Bias is the mean of posterior medians minus the target, with its Monte Carlo standard error (MCSE); Cov is the coverage of the central 95\% credible interval and Width its mean width; $\bar{\mathcal{O}}$ is the mean posterior-median overestimation ratio and Rel.\ the relative overestimation $(\icc_\eta - \icc_Y)/\icc_Y$ at the targets.}",
    r"\label{tab:S-studyA}",
    r"\small",
    r"\resizebox{\linewidth}{!}{%",
    r"\begin{tabular}{rcccrrrccrrrcccr}",
    r"\toprule",
    r"Cell & $\icc_\eta$ & $\sigma^2_\eta$ & $\icc_Y$ & Bias$_\eta$ & MCSE & RMSE$_\eta$ & Cov$_\eta$ & Bias$_Y$ & MCSE & RMSE$_Y$ & Cov$_Y$ & Width$_Y$ & $\bar{\mathcal{O}}$ & Rel. & $R$ \\",
    r"\midrule",
]
for c in A["cells"]:
    lines.append(" & ".join([
        str(c["cell"]), f2(c["icc_eta"]), f2(c["s2eta"]), f3(c["icc_Y"]),
        f4(c["bias_eta"]), f4(c["mcse_eta"]), f3(c["rmse_eta"]), f3(c["cov_eta"]),
        f4(c["bias_Y"]), f4(c["mcse_Y"]), f3(c["rmse_Y"]), f3(c["cov_Y"]), f3(c["width_Y"]),
        f2(c["O_mean"]), f"{100 * c['rel_over']:.0f}\\%", str(c["R"]),
    ]) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]
write("table_study_A_full.tex", lines)

# ---------------------------------------------------------------- Study B
B = num["study_B"]
lines = [
    r"\begin{table}[tp]",
    r"\centering",
    r"\caption{Study B in full ($n_p = 300$, $n_j = 20$, $\sigma^2_o = 0.3$, $R = 100$ per cell; no fit excluded). Targets: absolute and relative $\icc_Y$ by quadrature at the true parameters, $\icc_I$ by five long nested Monte Carlo runs (its own Monte Carlo error in the column MCSE$_I$). Bias, RMSE and coverage (Cov) of posterior medians and 95\% intervals for the four coefficients; $\bar\kappa$ and $\bar\varphi$ are the mean sample Fleiss kappa and mean pairwise phi of the raw data, to be compared with the absolute and relative targets.}",
    r"\label{tab:S-studyB}",
    r"\footnotesize",
    r"\setlength{\tabcolsep}{3.5pt}",
    r"\resizebox{\linewidth}{!}{%",
    r"\begin{tabular}{rlccccccrccrccrccrccc}",
    r"\toprule",
    r"Cell & Link & $\icc_\eta$ & $P$ & $\icc_Y$ & $\icc_Y^{\mathrm{rel}}$ & $\icc_I$ & MCSE$_I$ & Bias$_Y$ & RMSE$_Y$ & Cov$_Y$ & Bias$^{\mathrm{rel}}$ & RMSE$^{\mathrm{rel}}$ & Cov$^{\mathrm{rel}}$ & Bias$_\eta$ & RMSE$_\eta$ & Cov$_\eta$ & Bias$_I$ & RMSE$_I$ & Cov$_I$ & $\bar\kappa$ / $\bar\varphi$ \\",
    r"\midrule",
]
for c in B["cells"]:
    lines.append(" & ".join([
        str(c["cell"]), c["link"], f2(c["icc_eta"]), f2(c["P"]), f3(c["t_abs"]), f3(c["t_rel"]), f3(c["t_I"]), f4(c["t_I_mcse"]),
        f4(c["bias_abs"]), f3(c["rmse_abs"]), f2(c["cov_abs"]),
        f4(c["bias_rel"]), f3(c["rmse_rel"]), f2(c["cov_rel"]),
        f4(c["bias_eta"]), f3(c["rmse_eta"]), f2(c["cov_eta"]),
        f4(c["bias_I"]), f3(c["rmse_I"]), f2(c["cov_I"]),
        f"{f3(c['kappa_mean'])} / {f3(c['phi_mean'])}",
    ]) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]
write("table_study_B_full.tex", lines)

# link-invariance comparison table (pre-registered, reported as uninformative in the main text)
lines = [
    r"\begin{table}[tp]",
    r"\centering",
    r"\caption{The pre-registered comparison of the absolute $\icc_Y$ across links at matched $\icc_\eta$ and prevalence: population targets and biases under each link. The targets differ by construction because the two residual conventions define different populations; the comparison does not test link-invariance, which is a statement about one data set fitted under both links (Sections~8 and~9.3 of the main text).}",
    r"\label{tab:S-linkcompare}",
    r"\small",
    r"\begin{tabular}{cccccc}",
    r"\toprule",
    r"$\icc_\eta$ & $P$ & Target (probit) & Target (logit) & Bias (probit) & Bias (logit) \\",
    r"\midrule",
]
for r in B["link_invariance"]:
    lines.append(" & ".join([f2(r["icc_eta"]), f2(r["P"]), f3(r["t_abs_probit"]), f3(r["t_abs_logit"]),
                             f4(r["bias_abs_probit"]), f4(r["bias_abs_logit"])]) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
write("table_study_B_links.tex", lines)

# ---------------------------------------------------------------- Reviewer 2's design
rows = list(csv.DictReader(open(os.path.join(ROOT, "reviewer2", "vanbelle_design_results.csv"))))
cols = ["kappa_sample", "phi_pairwise", "probit_icc_Y_abs", "probit_icc_Y_rel", "probit_icc_eta", "probit_icc_I",
        "probit_pre_0_5_0", "logit_icc_Y_abs", "logit_icc_Y_rel", "logit_icc_eta", "logit_icc_I", "logit_pre_0_5_0"]
lines = [
    r"\begin{table}[tp]",
    r"\centering",
    r"\caption{Reproduction of the simulated design in Reviewer 2's R file (five replicates of her generating model, fitted under the probit and the logit link with \texttt{dgt} 0.5.0). Sample kappa and mean pairwise phi are computed from the raw data; the remaining columns are posterior medians from \texttt{dgt\_icc()}; ``pre-0.5.0'' is the value the package returned before the correction described in Section~8 of the main text. The last row is the mean over replicates.}",
    r"\label{tab:S-r2design}",
    r"\footnotesize",
    r"\setlength{\tabcolsep}{3.5pt}",
    r"\begin{tabular}{rcc|ccccc|ccccc}",
    r"\toprule",
    r" & & & \multicolumn{5}{c|}{Probit fit} & \multicolumn{5}{c}{Logit fit} \\",
    r"Rep. & $\hat\kappa$ & $\bar\varphi$ & $\icc_Y$ & $\icc_Y^{\mathrm{rel}}$ & $\icc_\eta$ & $\icc_I$ & pre-0.5.0 & $\icc_Y$ & $\icc_Y^{\mathrm{rel}}$ & $\icc_\eta$ & $\icc_I$ & pre-0.5.0 \\",
    r"\midrule",
]
for r in rows:
    lines.append(" & ".join([r["rep"]] + [f3(float(r[c])) for c in cols]) + r" \\")
lines.append(r"\midrule")
lines.append(" & ".join(["Mean"] + [f3(statistics.mean(float(r[c]) for r in rows)) for c in cols]) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
write("table_r2_design.tex", lines)

# ---------------------------------------------------------------- Vanbelle (2025) FHR reproduction
log = open(os.path.join(ROOT, "reproduction", "vanbelle-fhr", "console-log.txt")).read()
side = log[log.index("--- side by side ---"):]
side = side[:side.index("done.")]
entries = []
for line in side.splitlines()[2:]:
    if not line.strip():
        continue
    # columns are fixed-width: quantity | ours | vanbelle
    m = re.match(r"^\s(.{26})\s(.{20})\s(.*)$", line)
    if not m:
        continue
    q, ours, hers = (m.group(1).strip(), m.group(2).strip(), m.group(3).strip())
    entries.append((q, ours, hers))
kappa = re.search(r"Fleiss kappa \(intraclass kappa, binary\)\s*:\s*([0-9.]+)", log).group(1)
rhat = re.search(r"max Rhat\s*:\s*([0-9.]+)", log).group(1)
ess = re.search(r"min bulk ESS:\s*([0-9]+)", log).group(1)
div = re.search(r"divergent transitions:\s*([0-9]+)", log).group(1)
dgtv = re.search(r"dgt\s*:\s*([0-9.]+)", log).group(1)
lines = [
    r"\begin{table}[tp]",
    r"\centering",
    r"\caption{Reproduction of the fetal heart rate example of Vanbelle (2025, Section~6): 145 recordings classified as normal or abnormal by 16 midwives, probit crossed random-effects model with her priors, refitted in \texttt{brms} and passed to \texttt{dgt\_icc()} (\texttt{dgt} " + dgtv + r"; 12,000 post-warmup draws, maximum $\hat R$ " + rhat + r", minimum bulk ESS " + ess + r", " + div + r" divergent transitions). Posterior medians with 95\% intervals; her values as printed on her p.~14. The sample Fleiss kappa of the raw table is " + kappa + r".}",
    r"\label{tab:S-fhr}",
    r"\small",
    r"\begin{tabular}{lll}",
    r"\toprule",
    r"Quantity & This reproduction & Vanbelle (2025) \\",
    r"\midrule",
]
rename = {
    "sigma_R^2": r"$\sigma^2_R$ (midwife variance, probit scale)",
    "ICC_eta (latent)": r"$\icc_\eta$ (latent scale; her $\hat\tau_{\mathrm{abs}}$)",
    "ICC_Y absolute (manifest)": r"$\icc_Y$, absolute (manifest scale; her $\hat\varphi_{\mathrm{abs}}$)",
    "ICC_Y relative (manifest)": r"$\icc_Y$, relative",
    "Fleiss kappa (sample)": r"Fleiss kappa of the raw table (her $\hat\kappa_{\mathrm{abs}}$)",
    "mean pairwise phi (sample)": r"Mean pairwise phi of the raw table",
}
for q, ours, hers in entries:
    hers = hers.replace("-- (not reported)", "not reported").replace("--", "not reported").replace("_abs", r"$_{\mathrm{abs}}$").replace("tau$", r"$\tau$").replace("phi$", r"$\varphi$").replace("kappa$", r"$\kappa$")
    hers = re.sub(r"^(tau|phi|kappa)", lambda m: {"tau": r"$\hat\tau", "phi": r"$\hat\varphi", "kappa": r"$\hat\kappa"}[m.group(1)], hers)
    hers = hers.replace(r"$_{\mathrm{abs}}$ =", r"_{\mathrm{abs}}$ =")
    lines.append(f"{rename.get(q, q)} & {ours} & {hers} " + r"\\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
write("table_vanbelle_fhr.tex", lines)
