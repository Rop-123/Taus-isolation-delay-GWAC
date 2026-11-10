"""Builds figures/*.pdf and results/table_scenarios.tex from results/. Run after run_analysis.py."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, FIG = os.path.join(ROOT, "results"), os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "DejaVu Sans"})
COL = {7: "#4d4d4d", 5: "#e08a00", 3: "#0b6e99"}
EPS = [0.5, 0.8, 1.0]
TITLE = {0.5: "Worst: 50% effective", 0.8: "Base: 80% effective", 1.0: "Best: 100% effective"}

dr = np.load(os.path.join(RES, "draws.npz"))
res = pd.read_csv(os.path.join(RES, "scenario_results.csv"))
prcc = pd.read_csv(os.path.join(RES, "sensitivity_prcc.csv"), index_col=0).iloc[:, 0]
day0 = pd.Timestamp("2026-10-28")

# ---- Fig 2: daily incidence trajectories
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.35), sharey=True)
t = np.arange(365)
for ax, eps in zip(axs, EPS):
    for dl in [7, 5, 3]:
        inc = dr[f"inc_{eps}_{dl}"]
        lo, med, hi = np.percentile(inc, [2.5, 50, 97.5], axis=0)
        ax.fill_between(t, np.maximum(lo, 0.3), hi, color=COL[dl], alpha=0.18, lw=0)
        ax.plot(t, med, color=COL[dl], lw=1.6, label=f"{dl} days")
    ax.set_yscale("log"); ax.set_ylim(0.3, 800); ax.set_xlim(0, 365)
    ax.set_title(TITLE[eps], fontsize=8)
    ax.set_xlabel("Days from 28 Oct 2026")
axs[0].set_ylabel("New symptomatic cases per day")
axs[0].legend(title="Onset-to-isolation\ndelay", frameon=False, loc="lower left", fontsize=7, title_fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_trajectories.pdf")); plt.close(fig)

# ---- Fig 3: percent averted + PRCC
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.4), gridspec_kw={"width_ratios": [1.15, 1]})
w = 0.34
for j, dl in enumerate([5, 3]):
    sub = res[res.delay == dl].set_index("eps").loc[EPS]
    x = np.arange(3) + (j - 0.5) * w
    a.bar(x, sub.pct_averted, w, color=COL[dl], label=f"{dl} days")
    a.errorbar(x, sub.pct_averted, yerr=[sub.pct_averted - sub.pct_averted_lo, sub.pct_averted_hi - sub.pct_averted],
               fmt="none", ecolor="black", capsize=2, lw=0.9)
    for xi, v in zip(x, sub.pct_averted):
        a.text(xi, 3, f"{v:.0f}%", ha="center", color="white", fontsize=7, fontweight="bold")
a.set_xticks(range(3)); a.set_xticklabels(["50%", "80%", "100%"])
a.set_xlabel("Isolation effectiveness"); a.set_ylabel("Cases averted vs 7-day delay (%)")
a.set_ylim(0, 128); a.legend(title="New delay", frameon=False, fontsize=7, title_fontsize=7, loc="upper left", ncol=2)
a.set_title("(a) Cases averted in 12 months", fontsize=8)
p = prcc.sort_values()
b.barh(range(len(p)), p.values, color=["#b3412e" if v < 0 else "#0b6e99" for v in p.values])
b.set_yticks(range(len(p))); b.set_yticklabels([s.replace(" (transmission)", "") for s in p.index], fontsize=7)
b.axvline(0, color="black", lw=0.6); b.set_xlim(-1, 1)
b.set_xlabel("Partial rank correlation with 12-month cases")
b.set_title("(b) Global sensitivity", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig3_averted_sensitivity.pdf")); plt.close(fig)

# ---- Table (LaTeX)
def fmt(v, lo, hi, d=0):
    return f"{v:,.{d}f} ({lo:,.{d}f}--{hi:,.{d}f})"
lines = []
for eps in EPS:
    for dl in [7, 5, 3]:
        r = res[(res.eps == eps) & (res.delay == dl)].iloc[0]
        pct = "--" if dl == 7 else f"{r.pct_averted:.0f} ({r.pct_averted_lo:.0f}--{r.pct_averted_hi:.0f})"
        if dl == 7:
            dur = "no end (4\\%)" if r.p_ends_within_5y < 0.5 else f"{r.duration_d / 30.4:.0f}"
        else:
            dur = f"{r.duration_d / 30.4:.0f} ({r.duration_d_lo / 30.4:.0f}--{r.duration_d_hi / 30.4:.0f})"
            if r.p_ends_within_5y < 0.95:
                dur += f"$^\\dagger$"
        lines.append(f"{int(eps * 100)}\\% & {dl} & {r.Re:.2f} ({r.Re_lo:.2f}--{r.Re_hi:.2f}) & "
                     f"{r.cum365:,.0f} ({r.cum365_lo:,.0f}--{r.cum365_hi:,.0f}) & {pct} & {dur} \\\\")
    if eps != EPS[-1]:
        lines.append("\\midrule")
head = ("\\begin{tabular}{@{}rrcccc@{}}\\toprule\n"
        "$\\epsilon$ & Delay (d) & $R_e$ & Cases in 12 months & \\% averted vs 7 d & Duration (months)\\\\\\midrule\n")
open(os.path.join(RES, "table_scenarios.tex"), "w").write(head + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n")
print("figures and table written")
