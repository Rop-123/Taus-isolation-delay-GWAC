# How fast is fast enough? Onset-to-isolation delay and Bundibugyo virus disease (Team TAUS)

Deterministic SEIQR model (plus a stochastic check) quantifying how cutting the mean
symptom-onset-to-isolation delay from 7 to 5 and 3 days changes cumulative cases and outbreak
duration when isolation is 50%, 80% or 100% effective.

## Run
```bash
pip install -r requirements.txt     # numpy, scipy, pandas, matplotlib
Rscript -e 'install.packages("deSolve")'   # for R/fit_seiqr.R
./run_all.sh                        # ~2-3 min; needs pdflatex for the PDF
```
Outputs: `results/` (CSV/JSON/LaTeX table), `figures/` (PDF), `report/report.pdf`.
All seeds are fixed (`src/run_analysis.py`, seed 20261031).

## Layout
| Path | Content |
|---|---|
| `data/cumulative_confirmed_drc.csv` | Cumulative confirmed DRC cases with source for each row; rows flagged `used_for_calibration` set the current growth rate and starting incidence |
| `data/parameters.csv` | Every parameter, range, distribution and source |
| `src/model.py` | ODE model, calibration, discrete-time stochastic version |
| `src/run_analysis.py` | Scenarios (3 delays x 3 effectiveness), 1,000 parameter draws, stochastic runs, sensitivity |
| `src/fit_diagnostics.py` | Growth-rate fit, quasi-Poisson check, calibrated beta, hindcast against observed counts |
| `R/fit_seiqr.R` | **R fit** of the SEIQR model to the 2026 DRC surveillance data in `data/bvd_2026/` (onset curve with reporting-delay correction, confirmed counts by report date, changepoint in transmission, positivity check). Needs R + `deSolve` |
| `data/bvd_2026/` | Surveillance streams from epiforecasts/BVDOutbreakSize (INSP SitReps, data as of 12 Sep 2026); see its README for caveats and citation |
| `src/make_figures.py` | Figures and results table |
| `report/report.tex` | The report (reads `results/` and `figures/`) |

## Method in brief
* Compartments S, E, I (infectious, not isolated), Q (isolated), R; isolation rate tau = 1/delay;
  isolated people transmit with relative infectiousness (1 - effectiveness).
* Transmission is calibrated so that the **7-day baseline** reproduces the current growth rate of confirmed
  incidence (fitted to 4 intervals, 2 Jun-28 Oct 2026); each scenario then changes only tau.
* Uncertainty: 1,000 draws of latent period, infectious period and growth rate (Student-t, 2 d.f.);
  results are medians with 95% uncertainty intervals; cases averted are paired against the 7-day run of the same draw.
* "Cases" are confirmed symptomatic cases (under-ascertainment is not modelled). Horizon: 365 days from 28 Oct 2026.

## Caveats
Preliminary. See Limitations in the report; the main ones are the exponential-delay assumption, flat
confirmed incidence possibly reflecting testing capacity, and no explicit contact tracing or unsafe burials.
