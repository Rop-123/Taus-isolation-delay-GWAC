"""Model fitting / calibration diagnostics. Run:  python src/fit_diagnostics.py
Shows (1) the data and the growth-rate fit, (2) a second estimator (quasi-Poisson) as a robustness check,
(3) the calibrated transmission rate, (4) a hindcast of the model against observed cumulative counts.
Outputs: results/fit_diagnostics.json, results/hindcast.csv, figures/fig_fit.pdf"""
import json, os, sys
import numpy as np
import pandas as pd
from scipy import stats, optimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import model as M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, FIG = os.path.join(ROOT, "results"), os.path.join(ROOT, "figures")
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(7)

# ------------------------------------------------------------------ 1. data and log-linear fit
d = pd.read_csv(os.path.join(ROOT, "data", "cumulative_confirmed_drc.csv"), parse_dates=["date"])
cal = d[d.used_for_calibration == "yes"].sort_values("date").reset_index(drop=True)
t = (cal.date - cal.date.iloc[0]).dt.days.values.astype(float)       # day 0 = 2 Jun 2026
cum = cal.cum_confirmed_drc.values.astype(float)
dt_int = np.diff(t); counts = np.diff(cum); inc = counts / dt_int; mid = (t[1:] + t[:-1]) / 2
ols = stats.linregress(mid, np.log(inc))
tcrit = stats.t.ppf(0.975, len(inc) - 2)
ols_ci = (ols.slope - tcrit * ols.stderr, ols.slope + tcrit * ols.stderr)

# ------------------------------------------------------------------ 2. quasi-Poisson check (uses the counts, not the logs)
def expected(par):
    la, r = par
    a = np.exp(la)
    if abs(r) < 1e-9:
        return a * dt_int
    return a * (np.exp(r * (t[1:] - mid[0])) - np.exp(r * (t[:-1] - mid[0]))) / r
def nll(par):
    mu = expected(par)
    return np.sum(mu - counts * np.log(mu))
opt = optimize.minimize(nll, [np.log(inc.mean()), 0.0], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-10})
mu = expected(opt.x)
pearson = np.sum((counts - mu) ** 2 / mu); disp = pearson / (len(counts) - 2)
hess = np.array([optimize.approx_fprime(opt.x, lambda p, k=k: optimize.approx_fprime(p, nll, 1e-5)[k], 1e-4) for k in range(2)])
cov = np.linalg.inv(hess) * disp
qp_r, qp_se = opt.x[1], float(np.sqrt(cov[1, 1]))

# ------------------------------------------------------------------ 3. calibrated transmission rate (central parameters)
sigma, gamma = 1 / 6.3, 1 / 7.0
calib = []
for eps in [0.5, 0.8, 1.0]:
    beta = M.calibrate_beta(ols.slope, sigma, gamma, 1 / 7.0, eps, gamma)
    r_check = M.growth_rate(beta, sigma, gamma, 1 / 7.0, eps, gamma)
    calib.append(dict(eps=eps, beta=beta, R0_baseline=M.re_value(beta, sigma, gamma, 1 / 7.0, eps, gamma),
                      growth_rate_target=ols.slope, growth_rate_model=r_check))
calib = pd.DataFrame(calib)

# ------------------------------------------------------------------ 4. hindcast: start 2 Jun, compare with observed cumulative counts
NDR = 300
L = rng.triangular(4.0, 6.3, 10.1, NDR); D = rng.triangular(5.0, 7.0, 10.0, NDR)
R = np.clip(ols.slope + ols.stderr * rng.standard_t(2, NDR), -0.03, 0.03)
inc_start = inc[0]                                   # 63.7 per day in the first interval
horizon = int(t[-1])
paths = []
for i in range(NDR):
    s_, g_ = 1 / L[i], 1 / D[i]
    b = M.calibrate_beta(R[i], s_, g_, 1 / 7.0, 0.8, g_)
    y0 = M.initial_state(b, s_, g_, 1 / 7.0, 0.8, g_, inc_start)
    tt, y = M.simulate(b, s_, g_, 1 / 7.0, 0.8, g_, y0, horizon)
    paths.append(cum[0] + y[5])                      # observed cumulative at 2 Jun + model new cases
paths = np.array(paths)
rows = []
for k in range(1, len(t)):
    j = int(t[k]); lo, med, hi = np.percentile(paths[:, j], [2.5, 50, 97.5])
    rows.append(dict(date=cal.date[k].date(), observed=cum[k], model_median=med, model_lo=lo, model_hi=hi,
                     pct_error=100 * (med - cum[k]) / cum[k], inside_95UI=bool(lo <= cum[k] <= hi)))
hind = pd.DataFrame(rows); hind.to_csv(os.path.join(RES, "hindcast.csv"), index=False)
rmse_log = float(np.sqrt(np.mean((np.log(hind.model_median) - np.log(hind.observed)) ** 2)))

out = dict(interval_incidence=inc.round(2).tolist(),
           ols=dict(r=ols.slope, se=ols.stderr, ci95=ols_ci, r2=ols.rvalue ** 2, df=len(inc) - 2),
           quasi_poisson=dict(r=qp_r, se=qp_se, dispersion=disp, ci95=[qp_r - tcrit * qp_se, qp_r + tcrit * qp_se]),
           calibration=calib.round(5).to_dict("records"), hindcast_rmse_log=rmse_log,
           hindcast_all_inside=bool(hind.inside_95UI.all()))
json.dump(out, open(os.path.join(RES, "fit_diagnostics.json"), "w"), indent=2, default=float)

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4))
xs = np.linspace(mid[0], mid[-1], 100)
pred = ols.intercept + ols.slope * xs
se_line = np.sqrt(ols.stderr ** 2 * (xs - mid.mean()) ** 2 + (np.std(np.log(inc) - (ols.intercept + ols.slope * mid), ddof=2) ** 2) / len(mid))
ax[0].plot(mid, inc, "ko", label="observed")
ax[0].plot(xs, np.exp(pred), color="#0b6e99", label="log-linear fit")
ax[0].fill_between(xs, np.exp(pred - tcrit * se_line), np.exp(pred + tcrit * se_line), color="#0b6e99", alpha=0.2, lw=0)
ax[0].set_ylim(0, 130); ax[0].set_xlabel("Days from 2 Jun 2026"); ax[0].set_ylabel("Confirmed cases per day")
ax[0].set_title(f"(a) Growth-rate fit, r = {ols.slope:.4f}/d", fontsize=8); ax[0].legend(frameon=False, fontsize=7)
tt_ = np.arange(horizon + 1)
ax[1].fill_between(tt_, np.percentile(paths, 2.5, axis=0), np.percentile(paths, 97.5, axis=0), color="#e08a00", alpha=0.25, lw=0)
ax[1].plot(tt_, np.median(paths, axis=0), color="#e08a00", label="model (median, 95% UI)")
ax[1].plot(t, cum, "ko", label="observed")
ax[1].set_xlabel("Days from 2 Jun 2026"); ax[1].set_ylabel("Cumulative confirmed cases")
ax[1].set_title("(b) Hindcast from 2 June", fontsize=8); ax[1].legend(frameon=False, fontsize=7, loc="upper left")
ax[2].bar(np.arange(3), calib.beta, color="#4d4d4d", width=0.5)
ax[2].set_xticks(range(3)); ax[2].set_xticklabels(["50%", "80%", "100%"]); ax[2].set_xlabel("Isolation effectiveness")
ax[2].set_ylabel(r"Calibrated $\beta$ (per day)"); ax[2].set_title(r"(c) Calibrated transmission", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_fit.pdf")); plt.close(fig)

pd.set_option("display.width", 200)
print("Interval incidence (per day):", inc.round(1))
print(f"OLS log-linear:   r = {ols.slope:.5f}  SE = {ols.stderr:.5f}  95% CI = ({ols_ci[0]:.4f}, {ols_ci[1]:.4f})  R2 = {ols.rvalue**2:.2f}")
print(f"Quasi-Poisson:    r = {qp_r:.5f}  SE = {qp_se:.5f}  dispersion = {disp:.1f}")
print(calib.round(4).to_string(index=False))
print(hind.round(1).to_string(index=False))
print(f"Hindcast RMSE on log scale = {rmse_log:.3f}")
