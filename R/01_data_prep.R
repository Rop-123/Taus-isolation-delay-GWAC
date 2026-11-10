## 01_data_prep.R -- load the surveillance data, measure the reporting delay, build the two fitting datasets.
## Inputs : data/bvd_2026/*.csv  (epiforecasts/BVDOutbreakSize export, data as of 12 Sep 2026)
## Outputs: results/prep.rds, results/tables/data_summary.csv, results/positivity.csv
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
set.seed(20261101)

rd <- function(f, ...) read.csv(file.path(P$bvd, f), stringsAsFactors = FALSE, ...)
wide <- rd("bvd_2026_streams_wide.csv"); wide$date <- as.Date(wide$date)
lat  <- rd("bvd_2026_onset_curve_latest.csv"); lat$onset_date <- as.Date(lat$onset_date); REPORT <- as.Date(lat$report_date[1])
vint <- rd("bvd_2026_onset_curve_all_vintages.csv", colClasses = c(sitrep = "character"))
vint$onset_date <- as.Date(vint$onset_date); vint$report_date <- as.Date(vint$report_date)
sc <- rd("bvd_2026_scalars.csv")
cb <- sc[sc$name == "confirmed_break_dates", ]
break_date  <- as.Date(gsub("[^0-9-]", "", cb$value))
break_gross <- as.numeric(sub(".*'gross_cases': \\[(\\d+)\\].*", "\\1", cb$other_fields))
stopifnot(REPORT == CFG$cutoff, length(break_date) == 1, is.finite(break_gross))

## ---- 1. reporting delay (completeness of an onset day's count after L days), from all report vintages
ref <- vint[vint$sitrep == "121" & vint$onset_date <= as.Date("2026-08-15"), c("onset_date", "confirmed_total")]
mm <- merge(vint, ref, by = "onset_date", suffixes = c("", "_ref"))
mm$lag <- as.integer(mm$report_date - mm$onset_date); mm <- mm[mm$lag >= 0 & mm$lag <= 40, ]
agg <- aggregate(cbind(confirmed_total, confirmed_total_ref) ~ lag, mm, sum)
agg$ratio <- agg$confirmed_total / agg$confirmed_total_ref
wfit <- optim(c(1.5, 10), function(p) if (any(p <= 0)) 1e9 else sum((agg$ratio - pweibull(agg$lag + 0.5, p[1], p[2]))^2))
WK <- wfit$par[1]; WS <- wfit$par[2]
delay <- list(shape = WK, scale = WS, mean = WS * gamma(1 + 1 / WK),
              q50 = qweibull(.5, WK, WS) - .5, q90 = qweibull(.9, WK, WS) - .5, n_vintages = length(unique(vint$sitrep)),
              rmse = sqrt(mean((agg$ratio - pweibull(agg$lag + 0.5, WK, WS))^2)))

## ---- 2. Fit A data: confirmed cases by symptom-onset date (latest vintage), 25 May -- 6 Sep
dA <- lat[lat$onset_date >= as.Date("2026-05-25"), c("onset_date", "confirmed_total")]
dA$lag <- as.numeric(REPORT - dA$onset_date); dA$k <- as.numeric(dA$onset_date - CFG$t0)

## ---- 3. Fit B data: confirmed increments between report dates, 2 Jun -- 12 Sep, harmonisation step removed
cum <- wide[!is.na(wide$confirmed_case_history) & wide$date >= as.Date("2026-06-02"), c("date", "confirmed_case_history")]
dB <- data.frame(date = cum$date[-1], prev = cum$date[-nrow(cum)], y = diff(cum$confirmed_case_history))
dB$ndays <- as.numeric(dB$date - dB$prev)
step <- dB$y[dB$date == break_date]; retro <- step - break_gross; dB$y[dB$date == break_date] <- break_gross
dB$u <- as.numeric(dB$date - CFG$t0); dB$u0 <- as.numeric(dB$prev - CFG$t0)
tst <- wide[!is.na(wide$tests_analysed_daily_history), c("date", "tests_analysed_daily_history")]
dB$tests <- sapply(seq_len(nrow(dB)), function(j) {
  d <- tst$tests_analysed_daily_history[tst$date > dB$prev[j] & tst$date <= dB$date[j]]
  if (length(d)) mean(d) * dB$ndays[j] else NA })
dB$month <- format(dB$date, "%b")

## ---- 4. positivity diagnostic (new confirmed / samples analysed)
pos <- aggregate(cbind(y, tests) ~ month, dB[!is.na(dB$tests), ], sum)
pos$days <- aggregate(ndays ~ month, dB[!is.na(dB$tests), ], sum)$ndays
pos <- pos[match(c("Jun", "Jul", "Aug", "Sep"), pos$month), ]
pos$positivity <- 100 * pos$y / pos$tests; pos$tests_per_day <- pos$tests / pos$days
names(pos)[2:3] <- c("new_confirmed", "tests_analysed")
write.csv(pos, file.path(P$res, "positivity.csv"), row.names = FALSE)

## ---- 5. descriptive summary
last <- wide[wide$date == CFG$cutoff, ]
desc <- data.frame(
  quantity = c("confirmed cases (cumulative, 12 Sep)", "confirmed deaths (cumulative, 12 Sep)", "crude deaths/cases",
               "onset-curve cases (latest vintage)", "onset-curve share of confirmed", "report intervals in Fit B",
               "retrospective cases removed (22 Jul)", "mean new confirmed per day, 1-12 Sep"),
  value = c(last$confirmed_case_history, last$confirmed_death_history, last$confirmed_death_history / last$confirmed_case_history,
            sum(lat$confirmed_total), sum(lat$confirmed_total) / last$confirmed_case_history, nrow(dB), retro,
            (last$confirmed_case_history - wide$confirmed_case_history[wide$date == as.Date("2026-09-01")]) / 11))
write.csv(desc, file.path(P$tab, "data_summary.csv"), row.names = FALSE)

saveRDS(list(agg = agg, delay = delay, dA = dA, dB = dB, pos = pos, REPORT = REPORT, break_date = break_date,
             break_gross = break_gross, retro = retro, lat = lat, last = last, desc = desc), file.path(P$res, "prep.rds"))
cat(sprintf("Reporting delay: Weibull shape %.2f scale %.1f d, mean %.1f d, 50%% by %.0f d, 90%% by %.0f d (%d vintages)\n",
            WK, WS, delay$mean, delay$q50, delay$q90, delay$n_vintages))
cat(sprintf("Fit A: %d onset days | Fit B: %d report intervals | %d retrospective cases removed\n", nrow(dA), nrow(dB), retro))
print(desc); print(pos)
