"""Reproduces every number in the report. Run:  python src/run_analysis.py
Writes results/*.csv, results/*.json and results/draws.npz (used by make_figures.py)."""
import json, os, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import model as M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)
rng = np.random.default_rng(20261031)

DELAYS = [7, 5, 3]
EPS = [0.5, 0.8, 1.0]
EPS_LABEL = {0.5: "worst (50%)", 0.8: "base (80%)", 1.0: "best (100%)"}
HORIZON = 365
T_END = 1825
DT = 0.25
NDRAW = 1000
NSTOCH = 500

# ---------------------------------------------------------------- data -> r, inc0
d = pd.read_csv(os.path.join(ROOT, "data", "cumulative_confirmed_drc.csv"), parse_dates=["date"])
cal = d[d.used_for_calibration == "yes"].sort_values("date").reset_index(drop=True)
days = (cal.date - cal.date.iloc[0]).dt.days.values.astype(float)
inc = np.diff(cal.cum_confirmed_drc.values) / np.diff(days)       # confirmed per day per interval
mid = (days[1:] + days[:-1]) / 2
fit = stats.linregress(mid, np.log(inc))
r_hat, r_se = fit.slope, fit.stderr
inc0 = inc[-1]
# Quasi-Poisson-free but honest: with 4 intervals (2 d.f.) the SE is the regression SE.
data_summary = dict(interval_incidence_per_day=inc.round(1).tolist(), r_hat=r_hat, r_se=r_se,
                    r_ci95=[r_hat - stats.t.ppf(.975, 2) * r_se, r_hat + stats.t.ppf(.975, 2) * r_se],
                    inc0=inc0, last_date=str(cal.date.iloc[-1].date()))

# ---------------------------------------------------------------- parameter draws
def tri(lo, mode, hi, size):
    return rng.triangular(lo, mode, hi, size)

L = tri(4.0, 6.3, 10.1, NDRAW)      # latent period
D = tri(5.0, 7.0, 10.0, NDRAW)      # infectious period (also after isolation)
# only 4 intervals (2 d.f.) -> Student-t(2) draws around the fitted slope, clipped
R_ = np.clip(r_hat + r_se * rng.standard_t(2, NDRAW), -0.03, 0.03)   # current growth rate of incidence

def run_scenario(Li, Di, ri, delay, eps, tau0=1 / 7.0):
    sigma, gamma = 1 / Li, 1 / Di
    gi = gamma
    beta = M.calibrate_beta(ri, sigma, gamma, tau0, eps, gi)
    y0 = M.initial_state(beta, sigma, gamma, tau0, eps, gi, inc0)
    tau = 1.0 / delay
    t, y = M.simulate(beta, sigma, gamma, tau, eps, gi, y0, T_END)
    cum, dur, incd = M.outcomes(t, y, HORIZON)
    re = M.re_value(beta, sigma, gamma, tau, eps, gi)
    return dict(cum=cum, dur=dur, re=re, peak=incd.max(), tpeak=int(incd.argmax()), inc=incd, beta=beta)

def q(a):
    a = np.asarray(a, float)
    return np.nanmedian(a), np.nanpercentile(a, 2.5), np.nanpercentile(a, 97.5)

rows, store = [], {}
for eps in EPS:
    sims = {dl: [run_scenario(L[i], D[i], R_[i], dl, eps) for i in range(NDRAW)] for dl in DELAYS}
    for dl in DELAYS:
        s = sims[dl]
        cum = np.array([x["cum"] for x in s]); base = np.array([x["cum"] for x in sims[7]])
        dur = np.array([x["dur"] for x in s]); dur7 = np.array([x["dur"] for x in sims[7]])
        re = np.array([x["re"] for x in s]); peak = np.array([x["peak"] for x in s])
        tpk = np.array([x["tpeak"] for x in s])
        averted = base - cum
        pct = 100 * averted / base
        row = dict(eps=eps, scenario=EPS_LABEL[eps], delay=dl)
        for name, arr in [("cum365", cum), ("averted", averted), ("pct_averted", pct), ("duration_d", dur),
                          ("duration_shortened_d", dur7 - dur), ("Re", re), ("peak_inc", peak), ("t_peak", tpk)]:
            m, lo, hi = q(arr)
            row[name], row[name + "_lo"], row[name + "_hi"] = m, lo, hi
        row["p_ends_within_5y"] = float(np.mean(~np.isnan(dur)))
        row["p_Re_below_1"] = float(np.mean(re < 1))
        rows.append(row)
        store[f"cum_{eps}_{dl}"] = cum
        store[f"inc_{eps}_{dl}"] = np.array([x["inc"][:HORIZON] for x in s])
res = pd.DataFrame(rows)
res.to_csv(os.path.join(RES, "scenario_results.csv"), index=False)
np.savez_compressed(os.path.join(RES, "draws.npz"), **store)

# ---------------------------------------------------------------- stochastic check
srows = []
for eps in EPS:
    sig, gam = 1 / L[:NSTOCH], 1 / D[:NSTOCH]
    betas = np.array([M.calibrate_beta_discrete(R_[i], sig[i], gam[i], 1 / 7.0, eps, gam[i], DT) for i in range(NSTOCH)])
    y0 = np.array([M.initial_state_discrete(betas[i], sig[i], gam[i], 1 / 7.0, eps, gam[i], inc0, DT)
                   for i in range(NSTOCH)]).T
    for dl in DELAYS:
        C, ext = M.stochastic(betas, sig, gam, 1 / dl, eps, gam, [y0[k] for k in range(4)], T_END, rng, dt=DT)
        m, lo, hi = q(C) if len(C) else (np.nan,) * 3
        srows.append(dict(eps=eps, delay=dl, cum365_median=m, cum365_lo=lo, cum365_hi=hi,
                          p_extinct_5y=float(np.mean(~np.isnan(ext))),
                          ext_median_d=float(np.nanmedian(ext)) if np.any(~np.isnan(ext)) else np.nan))
pd.DataFrame(srows).to_csv(os.path.join(RES, "stochastic_results.csv"), index=False)

# ---------------------------------------------------------------- baseline-delay sensitivity
brows = []
for tau_base in [4.6, 7.0, 8.6]:
    for eps in [0.8]:
        for dl in [5, 3]:
            c, c0 = [], []
            for i in range(200):
                c0.append(run_scenario(L[i], D[i], R_[i], tau_base, eps, tau0=1 / tau_base)["cum"])
                c.append(run_scenario(L[i], D[i], R_[i], dl, eps, tau0=1 / tau_base)["cum"])
            c, c0 = np.array(c), np.array(c0)
            brows.append(dict(baseline_delay=tau_base, new_delay=dl, eps=eps, base_cum=np.median(c0),
                              cum=np.median(c), pct_averted=np.median(100 * (c0 - c) / c0)))
pd.DataFrame(brows).to_csv(os.path.join(RES, "baseline_delay_sensitivity.csv"), index=False)

# ---------------------------------------------------------------- global sensitivity (PRCC)
nS = 1000
Ls, Ds = tri(4.0, 6.3, 10.1, nS), tri(5.0, 7.0, 10.0, nS)
Rs = np.clip(r_hat + r_se * rng.standard_t(2, nS), -0.03, 0.03)
Es = rng.uniform(0.5, 1.0, nS)
Dl = rng.uniform(3.0, 7.0, nS)
out = np.array([run_scenario(Ls[i], Ds[i], Rs[i], Dl[i], Es[i])["cum"] for i in range(nS)])
X = pd.DataFrame({"latent period": Ls, "infectious period": Ds, "current growth rate (transmission)": Rs,
                  "isolation effectiveness": Es, "onset-to-isolation delay": Dl})
def prcc(X, y):
    Rk = X.rank(); yr = pd.Series(y).rank()
    res = {}
    for c in X.columns:
        others = [k for k in X.columns if k != c]
        A = np.column_stack([np.ones(len(X)), Rk[others].values])
        rx = Rk[c].values - A @ np.linalg.lstsq(A, Rk[c].values, rcond=None)[0]
        ry = yr.values - A @ np.linalg.lstsq(A, yr.values, rcond=None)[0]
        res[c] = float(np.corrcoef(rx, ry)[0, 1])
    return res
pr = prcc(X, out)
pd.Series(pr, name="PRCC").to_csv(os.path.join(RES, "sensitivity_prcc.csv"))

# ---------------------------------------------------------------- headline numbers
summary = dict(data=data_summary,
               draws=NDRAW, horizon_days=HORIZON,
               baseline_cum365_central={e: float(res[(res.eps == e) & (res.delay == 7)].cum365.iloc[0]) for e in EPS},
               prcc=pr)
with open(os.path.join(RES, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2, default=float)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
print(json.dumps(data_summary, indent=2, default=float))
print(res[["eps", "delay", "cum365", "cum365_lo", "cum365_hi", "averted", "pct_averted", "duration_d",
           "duration_d_lo", "duration_d_hi", "Re", "p_Re_below_1", "p_ends_within_5y"]].round(2))
print(pd.read_csv(os.path.join(RES, "stochastic_results.csv")).round(2))
print(pr)
print(pd.read_csv(os.path.join(RES, 'baseline_delay_sensitivity.csv')).round(2))
