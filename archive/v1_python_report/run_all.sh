#!/usr/bin/env bash
# Reproduces every number, figure and table in the report (about 2-3 minutes).
set -euo pipefail
cd "$(dirname "$0")"
python3 src/run_analysis.py        # model, scenarios, uncertainty, stochastic check, sensitivity -> results/
python3 src/fit_diagnostics.py     # calibration diagnostics -> results/fit_diagnostics.json, hindcast.csv, figures/fig_fit.pdf
Rscript R/fit_seiqr.R             # fit SEIQR to the BVD surveillance data in R -> results/fit_R_*.csv, figures/fig_fit_R.pdf
python3 src/make_figures.py        # figures/*.pdf and results/table_scenarios.tex
cd report && pdflatex -interaction=nonstopmode report.tex >/dev/null && pdflatex -interaction=nonstopmode report.tex >/dev/null
echo "Done: report/report.pdf"
