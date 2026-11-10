# Surveillance data for the 2026 DRC Bundibugyo virus outbreak

Cleaned, dated surveillance streams for the 2026 Bundibugyo virus disease (BVD)
outbreak in the Democratic Republic of the Congo, exported as CSV for the MPPR 2026
group projects.

- **Data as of:** 12 September 2026 (`as_of_date` in the upstream manifest)
- **Pulled:** 15 September 2026 from
  [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize),
  commit `6e3a0d0` (15 Sep 2026)
- **Live analysis and newer data:** <https://epiforecasts.io/BVDOutbreakSize/dev/>

Nothing here is a line list. Everything is a national aggregate read off the daily
situation reports (SitReps) of the DRC Institut National de Santé Publique (INSP),
either directly from the PDFs or via the INRB-UMIE transcriptions.

## Files

| File | What it is |
|---|---|
| `bvd_2026_streams_wide.csv` | One row per date, one column per stream. Blank means no observation that day. This is the file most groups will want. |
| `bvd_2026_streams_long.csv` | The same data in long format: `stream, date, value`. Easier for `ggplot2` / `dplyr`. |
| `bvd_2026_streams_metadata.csv` | One row per stream: whether it is a **cumulative** total, a **daily** count or a **stock** (level on the day), a one-line description, first and last date, and the upstream `source` note, which records exactly which SitRep, table and row each number came from and every correction applied. Read this before you fit anything. |
| `bvd_2026_scalars.csv` | Single-valued observations: cases exported to Uganda, the outbreak's break dates (see caveats), Ituri population and border-crossing volumes (from McCabe et al.), and the genomic estimate of the time to most recent common ancestor. Each with its source note. |
| `bvd_2026_onset_curve_latest.csv` | Confirmed cases by **symptom-onset date**, as published in the most recent SitRep (SitRep 121, 12 Sep 2026), split by outcome (`confirmed_alive`, `confirmed_dead`, `confirmed_total`). |
| `bvd_2026_onset_curve_all_vintages.csv` | The same curve as it appeared in every SitRep since it was first published (SitRep 059, 12 Jul 2026). One block per `sitrep` / `report_date`. Use it to see how much the recent tail fills in over time (the reporting triangle). |
| `insp_sitrep_scanned.csv` | The audit trail: one row per SitRep with the headline figures as read from the PDF and a free-text `notes` column recording misprints, errata and definition changes. |
| `export_streams.py` | The script that produced these files. To refresh: `git clone https://github.com/epiforecasts/BVDOutbreakSize && python export_streams.py BVDOutbreakSize .` |

## The streams

| Stream | Kind | Covers | What it measures |
|---|---|---|---|
| `confirmed_case_history` | cumulative | 14 May – 12 Sep | Laboratory-confirmed cases, national total |
| `confirmed_death_history` | cumulative | 14 May – 12 Sep | Deaths among confirmed cases |
| `recovered_history` | cumulative | 6 Jun – 12 Sep | Recovered / discharged patients |
| `isolation_history` | stock | 1 Jun – 12 Sep | Patients in isolation or hospitalised on the day |
| `bed_capacity_history` | stock | 9 Jun – 12 Sep | Isolation bed capacity, derived as occupancy ÷ reported occupancy rate |
| `tests_analysed_daily_history` | daily | 1 Jun – 12 Sep | Laboratory samples analysed in the previous 24 h |
| `suspected_daily_history` | daily | 4 Jun – 5 Aug | New suspected cases notified that day. **Stopped 6 Aug.** |
| `suspected_daily_deaths_history` | daily | 7 Jun – 11 Jul | New suspected deaths notified that day |
| `treatment_admissions_history` | daily | 13 Jun – 2 Aug | Admissions to treatment / transit / isolation centres. **Table discontinued 3 Aug.** |
| `treatment_deaths_history`, `treatment_ruleout_history`, `treatment_absconded_history` | daily | 13 Jun – 2 Aug | Deaths, non-cases discharged, and absconders from the same table |
| `treatment_aulit_history`, `treatment_confirmed_incare_history`, `treatment_suspect_incare_history` | stock | 13 Jun – 2 Aug | Patients in beds, split confirmed / suspected |
| `reported_case_history`, `death_history` | cumulative | 18 – 26 May | Suspected cases and suspected deaths, first two weeks only |
| `tests_received_history`, `tests_analysed_history` | cumulative | 23 – 28 May | Laboratory samples received and analysed, first two weeks only |

Full source notes for every stream are in `bvd_2026_streams_metadata.csv`.

## Sources

- **INSP situation reports** (Institut National de Santé Publique, DRC), the primary
  source. Read from the PDFs by the BVDOutbreakSize team and recorded in
  `insp_sitrep_scanned.csv`.
- **INRB-UMIE transcriptions** of the same SitReps:
  <https://github.com/INRB-UMIE/BDBV2026-Data>. Used to cross-check the cumulative
  confirmed-case and confirmed-death series.
- **WHO Disease Outbreak News** (2026-DON602, DON603) for the cases exported to Uganda.
- **McCabe et al.** for Ituri population and daily border-crossing volumes.
- Compiled, cleaned and documented by **epiforecasts/BVDOutbreakSize**
  (Abbott S, Sherratt K, et al.; MIT licence). Please cite the repository and the
  underlying INSP/INRB sources if you use these data.

## Caveats: read before modelling

1. **No single series covers the whole outbreak.** What the SitReps publish has
   changed several times. Daily suspected counts stop on 6 August; the
   treatment-centre patient-movement table stops on 3 August; the early suspected
   and laboratory cumulative series only cover late May. Check the first and last
   date of whatever you use.
2. **Cumulative totals, not incidence.** Most series are cumulative headline figures.
   Daily differences include backlogs and corrections, and days with no report leave
   gaps. Aggregate to weeks, or model the increments between reports, before fitting.
3. **A retrospective harmonisation step.** The jump in confirmed totals on 22 July 2026
   (SitRep 069) includes about 97 cases and 62 deaths that came from reconciling
   provincial databases with the national one. Those are not new infections on 22 July. The
   date is in `bvd_2026_scalars.csv` under `confirmed_break_dates`.
4. **"In isolation" changed meaning four times** (9 Jun, 19 Jun, 14 Jul, 11 Aug; see
   `occupancy_break_dates`). Level shifts on those days are reclassifications, not
   admissions or discharges. Bed capacity is *derived* from occupancy and the printed
   occupancy rate, so it inherits the same breaks.
5. **The onset curve is digitised from a picture.** INSP publishes cases by
   symptom-onset date only as a bar-chart image. The counts were recovered from the
   pixels and are accurate to a few percent overall, with roughly ±1–2 cases of noise
   per day. The most recent two to three weeks in any vintage are incomplete because
   of reporting delay (right truncation). Fit to the latest vintage without handling
   that and you will conclude the outbreak is ending. `bvd_2026_onset_curve_all_vintages.csv`
   lets you measure the delay.
6. **Testing grew faster than the epidemic.** Samples analysed per day rose from around
   100 in early June to over 400 by September. Confirmed counts reflect laboratory
   capacity as well as incidence. Suspected and confirmed cases are two different case
   definitions with different biases; do not mix them in one likelihood.
7. **Deaths ÷ cases is not the case-fatality ratio.** Confirmed deaths (3,475) over
   confirmed cases (7,200) comes out near 50%, but outcome follow-up is incomplete
   and who gets tested is not random. Treat any CFR you compute from these files as
   an upper bound with a caveat.
8. **What is not here.** No line list, no contact-tracing data, no age or sex
   breakdown, no sub-national series, no vaccination (there is no licensed vaccine
   for Bundibugyo virus). Movement data come from the literature, not from this
   outbreak. If your question needs any of these, simulate them and say so in the
   proposal.

Never paste these files, or any partner or surveillance data, into a public AI tool.
