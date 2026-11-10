"""Export the BVD 2026 surveillance streams from epiforecasts/BVDOutbreakSize as CSV.
Usage: python export_streams.py REPO_DIR OUT_DIR"""
import sys, tomllib, subprocess, datetime as dt
import pandas as pd

REPO, OUT = sys.argv[1], sys.argv[2]
obs = tomllib.load(open(f"{REPO}/data/observations.toml", "rb"))
as_of = str(obs["as_of_date"])
commit = subprocess.check_output(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], text=True).strip()
commit_date = subprocess.check_output(["git", "-C", REPO, "log", "-1", "--format=%cd", "--date=short"], text=True).strip()

# Short, human descriptions of each dated stream (source notes come from the manifest itself).
DESC = {
    "confirmed_case_history": ("cumulative", "Cumulative laboratory-confirmed cases, national (INSP front-page 'Cumul cas confirmés')"),
    "confirmed_death_history": ("cumulative", "Cumulative deaths among confirmed cases, national"),
    "reported_case_history": ("cumulative", "Cumulative suspected cases, national, first two weeks only"),
    "death_history": ("cumulative", "Cumulative suspected deaths, national, first two weeks only"),
    "recovered_history": ("cumulative", "Cumulative recovered (discharged) patients"),
    "suspected_daily_history": ("daily", "New suspected cases notified that day (stopped 6 Aug 2026)"),
    "suspected_daily_deaths_history": ("daily", "New suspected deaths notified that day"),
    "isolation_history": ("stock", "Patients in isolation / hospitalised that day (definition changed on the break dates)"),
    "bed_capacity_history": ("stock", "Isolation bed capacity, derived as occupancy / occupancy rate"),
    "tests_received_history": ("cumulative", "Cumulative laboratory samples received (late May only)"),
    "tests_analysed_history": ("cumulative", "Cumulative laboratory samples analysed (late May only)"),
    "tests_analysed_daily_history": ("daily", "Laboratory samples analysed in the previous 24 h"),
    "treatment_admissions_history": ("daily", "Admissions to treatment / transit / isolation centres that day (table discontinued 3 Aug 2026)"),
    "treatment_deaths_history": ("daily", "Deaths in treatment centres that day"),
    "treatment_ruleout_history": ("daily", "Patients discharged as non-cases from treatment centres that day"),
    "treatment_absconded_history": ("daily", "Patients who absconded from treatment centres that day"),
    "treatment_aulit_history": ("stock", "Patients in beds in treatment centres that day"),
    "treatment_confirmed_incare_history": ("stock", "Confirmed patients in care in treatment centres that day"),
    "treatment_suspect_incare_history": ("stock", "Suspected patients in care in treatment centres that day"),
}

long_rows, meta_rows, scalar_rows = [], [], []
for key, d in obs.items():
    if not isinstance(d, dict):
        continue
    if "dates" in d and "values" in d:
        kind, desc = DESC.get(key, ("", ""))
        dates = [str(x) for x in d["dates"]]
        for dd, v in zip(dates, d["values"]):
            long_rows.append(dict(stream=key, date=dd, value=v))
        meta_rows.append(dict(stream=key, kind=kind, description=desc, first_date=dates[0], last_date=dates[-1],
                              n_obs=len(dates), last_value=d["values"][-1], source=str(d.get("source", "")).replace("\n", " ")))
    else:
        val = d.get("value", d.get("date"))
        extra = {k: v for k, v in d.items() if k not in ("value", "source")}
        scalar_rows.append(dict(name=key, value=str(val), other_fields=str(extra) if extra else "",
                                source=str(d.get("source", "")).replace("\n", " ")))

long = pd.DataFrame(long_rows)
long["date"] = pd.to_datetime(long["date"])
long = long.sort_values(["stream", "date"])
wide = long.pivot(index="date", columns="stream", values="value").sort_index()
wide = wide[[m["stream"] for m in meta_rows if m["stream"] in wide.columns]].astype("Int64")  # manifest order, integers

long.to_csv(f"{OUT}/bvd_2026_streams_long.csv", index=False, date_format="%Y-%m-%d")
wide.to_csv(f"{OUT}/bvd_2026_streams_wide.csv", date_format="%Y-%m-%d")
pd.DataFrame(meta_rows).to_csv(f"{OUT}/bvd_2026_streams_metadata.csv", index=False)
pd.DataFrame(scalar_rows).to_csv(f"{OUT}/bvd_2026_scalars.csv", index=False)

# onset curve: every report vintage, plus the latest vintage on its own
onset = pd.read_csv(f"{REPO}/data/onset_curve_scanned.csv", dtype={"sitrep": str})
onset.to_csv(f"{OUT}/bvd_2026_onset_curve_all_vintages.csv", index=False)
latest_sitrep = onset.sort_values("report_date")["sitrep"].iloc[-1]
onset[onset["sitrep"] == latest_sitrep].to_csv(f"{OUT}/bvd_2026_onset_curve_latest.csv", index=False)

# the audit trail behind the PDF-sourced streams
pd.read_csv(f"{REPO}/data/insp_sitrep_scanned.csv").to_csv(f"{OUT}/insp_sitrep_scanned.csv", index=False)

print(f"as_of={as_of} commit={commit} ({commit_date}) streams={len(meta_rows)} scalars={len(scalar_rows)} "
      f"onset_vintages={onset['sitrep'].nunique()} latest_sitrep={latest_sitrep}")
for m in meta_rows:
    print(f"  {m['stream']:38s} {m['kind']:10s} {m['first_date']}..{m['last_date']}  n={m['n_obs']:3d}  last={m['last_value']}")
