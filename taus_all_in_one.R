## taus_all_in_one.R -- the whole analysis in ONE script (generated from R/*.R by R/make_single_script.R; do not edit by hand).
## How fast is fast enough? Onset-to-isolation delay, 2026 Bundibugyo virus disease outbreak, DRC (Team TAUS).
##
## Run from anywhere:   Rscript taus_all_in_one.R              # tests, all analysis, figures, tables and the numbers quoted in the reports
##                      Rscript taus_all_in_one.R --pdf        # additionally build the three report PDFs (needs pdflatex)
##
## CONTENTS: 00 shared model code | TESTS | 01 data prep | 02 model fitting | 03 scenarios | 04 stochastic | 05 sensitivity | 06 validation | 07 figures, tables, numbers
## Needs: R >= 4.3 + package deSolve; for the PDFs pdflatex + bibtex. Runtime about 12-15 min on one core. All seeds are fixed.
## Reads data/ ; writes results/, figures/, manuscript/numbers.tex and reports/out/*.pdf. Folder layout is explained in README.md.
## The workspace is cleared between sections (as when the separate scripts ran in separate processes); sections pass results via results/*.rds.
suppressPackageStartupMessages(library(deSolve)); T0 <- Sys.time(); DO_PDF <- "--pdf" %in% commandArgs(TRUE)
ROOT <- local({ f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)); r <- if (length(f)) dirname(normalizePath(f[1])) else normalizePath(".")
  if (!dir.exists(file.path(r, "data", "bvd_2026"))) stop("Cannot find data/bvd_2026 under ", r, " -- run from the repository root.", call. = FALSE); r })
setwd(ROOT); banner <- function(s) cat(sprintf("\n=== %s  [%.1f min elapsed]\n", s, as.numeric(difftime(Sys.time(), T0, units = "mins"))))

######################################################################################################################
## 00 | Shared model code, constants and helpers (global environment)
######################################################################################################################
## 00_functions.R -- shared model code, helpers and constants (sourced by every script)
## Model: S-E-I-Q-R with an isolated compartment Q. The infectious (not yet isolated) class is split into
## m Erlang stages (m = 1 gives the plain exponential delay used in the main analysis; m > 1 is a
## structural sensitivity analysis). Isolation occurs on leaving stage m at rate m*tau (tau = 1/mean delay).

## ------------------------------------------------------------------ paths
P <- list(data = file.path(ROOT, "data"), bvd = file.path(ROOT, "data", "bvd_2026"),
          res = file.path(ROOT, "results"), fig = file.path(ROOT, "figures"),
          tab = file.path(ROOT, "results", "tables"))
for (d in c(P$res, P$fig, P$tab)) dir.create(d, FALSE, TRUE)

## ------------------------------------------------------------------ constants (single source of truth)
CFG <- list(
  sigma0 = 1 / 6.3,       # central latent period 6.3 d            (Ebengo, Akilimali et al. 2026 model input)
  gamma0 = 1 / 7,         # central infectious period 7 d          (idem; 6.5 d in Chowell & Nishiura 2014)
  lat_tri = c(4.0, 6.3, 10.1),   # triangular (min, mode, max) latent period, days
  inf_tri = c(5.0, 7.0, 10.0),   # triangular (min, mode, max) infectious period, days
  tau0 = 1 / 7,           # baseline onset-to-isolation delay 7 d (study brief)
  delays = c(7, 5, 3),    # scenario delays, days
  eps = c(0.5, 0.8, 1.0), # isolation effectiveness: worst / base / best
  eps_base = 0.8,
  pop = 4392200,          # Ituri population (scalars file); depletion negligible
  t0 = as.Date("2026-05-02"),      # model start for the fit
  cutoff = as.Date("2026-09-12"),  # last day of the surveillance data = start of scenarios
  horizon = 365,          # 12-month horizon for cumulative cases (days)
  max_days = 1825,        # 5-year window for outbreak duration
  cfr = 0.483             # 3,475 / 7,200 confirmed deaths / cases at 12 Sep (upper bound; README caveat 7)
)

## ------------------------------------------------------------------ small helpers
rtri <- function(n, p) {                      # triangular(min, mode, max) by inverse CDF
  u <- runif(n); a <- p[1]; c <- p[2]; b <- p[3]; Fc <- (c - a) / (b - a)
  ifelse(u < Fc, a + sqrt(u * (b - a) * (c - a)), b - sqrt((1 - u) * (b - a) * (b - c)))
}
q3 <- function(x, na.rm = TRUE) quantile(x, c(.5, .025, .975), na.rm = na.rm, names = FALSE, type = 1)
fmt <- function(x, d = 0, big = TRUE) formatC(x, format = "f", digits = d, big.mark = if (big) "," else "")
mci <- function(x, d = 0, big = TRUE) { q <- q3(x); sprintf("%s (%s--%s)", fmt(q[1], d, big), fmt(q[2], d, big), fmt(q[3], d, big)) }
tex_escape_pct <- function(s) gsub("%", "\\\\%", s)

## ------------------------------------------------------------------ deterministic model
## NB: layout of the state vector for the Jacobian is (E, I_1..I_m, Q); S, R and the case counter are not needed.
jac_m <- function(beta, sigma, gamma, tau, eps, m = 1) {
  n <- m + 2; J <- matrix(0, n, n); r <- m * tau
  J[1, 1] <- -sigma; J[1, 2:(1 + m)] <- beta; J[1, n] <- beta * (1 - eps)
  J[2, 1] <- sigma
  for (k in 1:m) { J[1 + k, 1 + k] <- -(r + gamma); if (k > 1) J[1 + k, k] <- r }
  J[n, 1 + m] <- r; J[n, n] <- -gamma
  J
}
growth_rate <- function(beta, sigma, gamma, tau, eps, m = 1) max(Re(eigen(jac_m(beta, sigma, gamma, tau, eps, m), only.values = TRUE)$values))
## Reproduction number at the disease-free state: R = (beta/gamma) * (1 - eps * a^m),  a = m*tau/(m*tau + gamma)
Re_value <- function(beta, gamma, tau, eps, m = 1) { a <- m * tau / (m * tau + gamma); (beta / gamma) * (1 - eps * a^m) }
frac_isolated <- function(gamma, tau, m = 1) (m * tau / (m * tau + gamma))^m       # share of cases isolated before recovery/death

calibrate_beta <- function(r, sigma, gamma, tau0, eps, m = 1)
  uniroot(function(b) growth_rate(b, sigma, gamma, tau0, eps, m) - r, c(1e-4, 6), tol = 1e-12)$root

make_rhs <- function(m, sigma, gamma, pop) {
  function(t, y, p) {                       # p = c(beta, tau_new, eps, tau0, ramp_days)
    beta <- p[1]; eps <- p[3]
    tau <- if (p[5] > 0) p[4] + (p[2] - p[4]) * min(1, t / p[5]) else p[2]
    S <- y[1]; E <- y[2]; I <- y[3:(2 + m)]; Q <- y[3 + m]
    lam <- beta * (sum(I) + (1 - eps) * Q) / pop
    rate <- m * tau
    dI <- c(sigma * E, rate * I[-m]) - (rate + gamma) * I
    list(c(-lam * S, lam * S - sigma * E, dI, rate * I[m] - gamma * Q, gamma * (sum(I) + Q), sigma * E))
  }
}
## Initial state: dominant eigenvector of the linearised (E, I_1..I_m, Q) system at the BASELINE delay, scaled so that sigma*E = i0
init_state <- function(beta, sigma, gamma, tau0, eps, i0, m = 1, pop = CFG$pop) {
  ev <- eigen(jac_m(beta, sigma, gamma, tau0, eps, m)); v <- abs(Re(ev$vectors[, which.max(Re(ev$values))]))
  v <- v * i0 / (sigma * v[1])
  c(S = pop - sum(v), E = v[1], I = v[2:(1 + m)], Q = v[2 + m], R = 0, C = 0)
}
## Run one scenario. Returns daily symptomatic incidence (length tend) and active infections (E+I+Q) at each day.
run_det <- function(beta, sigma, gamma, tau_new, eps, i0, tau0 = CFG$tau0, m = 1, tend = CFG$max_days,
                    ramp = 0, pop = CFG$pop, y0 = NULL) {
  if (is.null(y0)) y0 <- init_state(beta, sigma, gamma, tau0, eps, i0, m, pop)
  o <- ode(y0, 0:tend, make_rhs(m, sigma, gamma, pop), c(beta, tau_new, eps, tau0, ramp),
           method = "lsoda", rtol = 1e-7, atol = 1e-7)
  nm <- m + 5
  list(inc = diff(o[, nm + 1]), active = rowSums(o[, 3:(3 + m + 1)]) , cum = o[, nm + 1])
}
## active infections = E + sum(I_k) + Q = columns 3 .. (3+m+1) of the ode matrix (col 1 = time, col 2 = S)
outcomes <- function(sim, horizon = CFG$horizon) {
  ext <- which(sim$active < 1)[1]                      # first day with < 1 exposed/infectious person
  list(cum = sim$cum[horizon + 1] - sim$cum[1], duration = if (is.na(ext)) Inf else ext - 1,
       peak = max(sim$inc[1:horizon]), tpeak = which.max(sim$inc[1:horizon]))
}

## ------------------------------------------------------------------ discrete-time stochastic model (m = 1)
disc_matrix <- function(beta, sigma, gamma, tau, eps, dt) {
  es <- exp(-sigma * dt); et <- exp(-tau * dt); eg <- exp(-gamma * dt)
  matrix(c(es, beta * dt, beta * dt * (1 - eps), 1 - es, et * eg, 0, 0, 1 - et, eg), 3, 3, byrow = TRUE)
}
calibrate_beta_disc <- function(r, sigma, gamma, tau0, eps, dt)
  uniroot(function(b) max(Re(eigen(disc_matrix(b, sigma, gamma, tau0, eps, dt), only.values = TRUE)$values)) - exp(r * dt),
          c(1e-4, 6), tol = 1e-12)$root
init_disc <- function(beta, sigma, gamma, tau0, eps, i0, dt, pop = CFG$pop) {
  ev <- eigen(disc_matrix(beta, sigma, gamma, tau0, eps, dt)); v <- abs(Re(ev$vectors[, which.max(Re(ev$values))]))
  v * i0 * dt / ((1 - exp(-sigma * dt)) * v[1])
}
## Vectorised over replicates: every argument may be a vector of length nrep. Binomial-chain (tau-leap) transitions.
run_stoch <- function(beta, sigma, gamma, tau, eps, E0, I0, Q0, tend = CFG$max_days, dt = 0.25, pop = CFG$pop,
                      horizon = CFG$horizon) {
  n <- length(beta); S <- rep(round(pop), n); E <- round(E0); I <- round(I0); Q <- round(Q0)
  C12 <- numeric(n); ext <- rep(NA_real_, n); steps <- as.integer(tend / dt)
  pE <- 1 - exp(-sigma * dt); pQ <- 1 - exp(-tau * dt); pG <- 1 - exp(-gamma * dt)
  for (k in seq_len(steps)) {
    lam <- beta * (I + (1 - eps) * Q) / pop
    nE <- rbinom(n, S, 1 - exp(-lam * dt)); nI <- rbinom(n, E, pE)
    nQ <- rbinom(n, I, pQ); rI <- rbinom(n, I - nQ, pG); rQ <- rbinom(n, Q, pG)
    S <- S - nE; E <- E + nE - nI; I <- I + nI - nQ - rI; Q <- Q + nQ - rQ
    if (k * dt <= horizon) C12 <- C12 + nI
    done <- is.na(ext) & (E + I + Q == 0); ext[done] <- k * dt
  }
  list(cum12 = C12, ext = ext)
}

## ------------------------------------------------------------------ writing numbers for LaTeX
MACROS <- new.env()
macro_reset <- function() rm(list = ls(MACROS), envir = MACROS)
macro <- function(name, value) assign(name, as.character(value), envir = MACROS)
macro_write <- function(file) {
  nm <- sort(ls(MACROS)); lines <- vapply(nm, function(n) sprintf("\\newcommand{\\%s}{%s}", n, get(n, MACROS)), "")
  writeLines(c("% AUTO-GENERATED by R/07_figures_tables.R -- do not edit", lines), file)
}

KEEP <- c(ls(all.names = TRUE), "KEEP")   # objects that survive between sections (shared model code, constants, helpers)

######################################################################################################################
## TESTS | Unit tests of the model code (stop on failure)
######################################################################################################################
banner("unit tests")
## Unit tests for the model code. Run: Rscript tests/test_model.R   (stops with an error if any check fails)
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg); cat("pass:", msg, "\n") }
sg <- CFG$sigma0; gm <- CFG$gamma0
for (m in c(1, 3, 10)) for (eps in c(0.5, 0.8, 1)) for (tau in c(1/7, 1/3)) {
  b <- calibrate_beta(0, sg, gm, tau, eps, m)                      # beta giving zero growth
  ok(abs(Re_value(b, gm, tau, eps, m) - 1) < 1e-8, sprintf("Re = 1 exactly when growth rate = 0 (m=%d, eps=%.1f, delay=%.0f)", m, eps, 1/tau))
}
ok(abs(frac_isolated(gm, 1/7, 1) - 0.5) < 1e-12, "exponential 7-d delay isolates 50% before removal")
ok(abs(frac_isolated(gm, 1/7, 1000) - exp(-1)) < 1e-3, "Erlang limit -> exp(-gamma*D)")
## numerical growth rate of the ODE equals the Jacobian eigenvalue
r <- 0.01; b <- calibrate_beta(r, sg, gm, CFG$tau0, 0.8, 1)
big <- 1e12                                   # no depletion: ODE growth must equal the eigenvalue
s <- run_det(b, sg, gm, CFG$tau0, 0.8, 100, tend = 200, pop = big, y0 = init_state(b, sg, gm, CFG$tau0, 0.8, 100, 1, big))
ok(abs(log(s$inc[200] / s$inc[100]) / 100 - r) < 1e-5, "ODE growth rate matches calibration target (no depletion)")
s <- run_det(b, sg, gm, CFG$tau0, 0.8, 100, tend = 200)
ok(log(s$inc[200] / s$inc[100]) / 100 < r, "susceptible depletion slows growth slightly (as it should)")
## paired scenario runs: faster isolation always reduces cases
c7 <- outcomes(run_det(b, sg, gm, 1/7, 0.8, 100, tend = 400))$cum; c3 <- outcomes(run_det(b, sg, gm, 1/3, 0.8, 100, tend = 400))$cum
ok(c3 < c7, "3-day delay gives fewer cases than 7-day")
## effectiveness 0 => delay irrelevant
b0 <- calibrate_beta(r, sg, gm, CFG$tau0, 0, 1)
ok(abs(outcomes(run_det(b0, sg, gm, 1/3, 0, 100, tend = 400))$cum - outcomes(run_det(b0, sg, gm, 1/7, 0, 100, tend = 400))$cum) < 1e-3 * c7, "eps = 0 makes delay irrelevant")
## stochastic mean ~ deterministic
set.seed(1); n <- 400; dt <- 0.25
bd <- calibrate_beta_disc(0.01, sg, gm, CFG$tau0, 0.8, dt); v <- init_disc(bd, sg, gm, CFG$tau0, 0.8, 100, dt)
st <- run_stoch(rep(bd, n), sg, gm, 1/7, 0.8, rep(v[1], n), rep(v[2], n), rep(v[3], n), tend = 365)
dt_cum <- outcomes(run_det(b, sg, gm, 1/7, 0.8, 100, tend = 400))$cum
ok(abs(mean(st$cum12) / dt_cum - 1) < 0.05, sprintf("stochastic mean (%.0f) within 5%% of deterministic (%.0f)", mean(st$cum12), dt_cum))
cat("ALL TESTS PASSED\n")
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 01 | Data preparation: reporting-delay curve and the two fitting datasets
######################################################################################################################
banner("01 | Data preparation: reporting-delay curve and the two fitting datasets")
## 01_data_prep.R -- load the surveillance data, measure the reporting delay, build the two fitting datasets.
## Inputs : data/bvd_2026/*.csv  (epiforecasts/BVDOutbreakSize export, data as of 12 Sep 2026)
## Outputs: results/prep.rds, results/tables/data_summary.csv, results/positivity.csv
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
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 02 | Maximum-likelihood fits, change-point profile, parameter draws, adequacy checks
######################################################################################################################
banner("02 | Maximum-likelihood fits, change-point profile, parameter draws, adequacy checks")
## 02_fit.R -- fit the SEIQR model to the 2026 DRC surveillance data (maximum likelihood, negative-binomial counts).
## Natural-history parameters are fixed at central literature values; estimated: starting incidence, transmission rate
## (optionally changing once on an estimated date) and count dispersion. Uncertainty: profile likelihood over the change
## date + multivariate-normal (Hessian) approximation for the other parameters.
## Outputs: results/fit.rds, results/fit_draws_{A,B}.csv, results/tables/fit_*.csv
set.seed(20261101); t_start <- Sys.time()
pr <- readRDS(file.path(P$res, "prep.rds")); dA <- pr$dA; dB <- pr$dB
TEND <- as.numeric(CFG$cutoff - CFG$t0); WK <- pr$delay$shape; WS <- pr$delay$scale
cc <- function(L) pweibull(L + 0.5, WK, WS)                   # share of an onset day's cases already reported after L days
pm_full <- function(n) { pm <- diff(c(0, cc(0:(n - 1)))); M <- matrix(0, n, n); for (u in 1:n) M[u, 1:u] <- pm[u:1]; M }
PM <- pm_full(TEND)

## predictions on the observation scale
predA <- function(inc) inc[dA$k] * cc(dA$lag)
predB <- function(inc) { Cc <- c(0, cumsum(as.vector(PM %*% inc))); Cc[dB$u + 1] - Cc[dB$u0 + 1] }

## piecewise-constant transmission: beta1 until day tc, beta2 afterwards; returns daily incidence on days 1..TEND
sim_fit <- function(i0, b1, b2, tc, eps = CFG$eps_base, sigma = CFG$sigma0, gamma = CFG$gamma0, tend = TEND) {
  rhs <- make_rhs(1, sigma, gamma, CFG$pop); y0 <- init_state(b1, sigma, gamma, CFG$tau0, eps, i0, 1)
  o1 <- ode(y0, 0:min(tc, tend), rhs, c(b1, CFG$tau0, eps, CFG$tau0, 0), method = "lsoda", rtol = 1e-6, atol = 1e-6)
  C <- o1[, 7]
  if (tc < tend) {
    o2 <- ode(o1[nrow(o1), -1], tc:tend, rhs, c(b2, CFG$tau0, eps, CFG$tau0, 0), method = "lsoda", rtol = 1e-6, atol = 1e-6)
    C <- c(C, o2[-1, 7])
  }
  diff(C)
}
nll_gen <- function(p4, tc, pred, y, ...) {                 # p4 = log i0, log beta1, log beta2, log k
  inc <- tryCatch(sim_fit(exp(p4[1]), exp(p4[2]), exp(p4[3]), tc, ...), error = function(e) NULL)
  if (is.null(inc) || any(!is.finite(inc))) return(1e10)
  -sum(dnbinom(y, mu = pmax(pred(inc), 1e-8), size = exp(p4[4]), log = TRUE))
}
nm2 <- function(f, start, it = 4000) { o <- optim(start, f, control = list(maxit = it, reltol = 1e-10)); optim(o$par, f, control = list(maxit = it, reltol = 1e-12)) }
fit_const  <- function(pred, y, start, ...) { f <- function(p) nll_gen(c(p[1], p[2], p[2], p[3]), TEND, pred, y, ...); o <- nm2(f, start)
  list(par4 = c(o$par[1], o$par[2], o$par[2], o$par[3]), nll = o$value, hess = optimHess(o$par, f), tc = TEND) }
fit_change <- function(pred, y, tc, start4, ...) { f <- function(p) nll_gen(p, tc, pred, y, ...); o <- nm2(f, start4)
  list(par4 = o$par, nll = o$value, hess = optimHess(o$par, f), tc = tc) }
mvn_draw <- function(mean, hess) {
  cv <- tryCatch(solve(hess), error = function(e) diag(1 / pmax(diag(hess), 1e-6))); cv <- (cv + t(cv)) / 2
  e <- eigen(cv); as.vector(mean + e$vectors %*% (sqrt(pmax(e$values, 1e-8)) * rnorm(length(mean))))
}
gr <- function(b, eps = CFG$eps_base, sigma = CFG$sigma0, gamma = CFG$gamma0) growth_rate(b, sigma, gamma, CFG$tau0, eps, 1)

analyse <- function(label, pred, y, start_const) {
  cat("\n==== ", label, " ====\n", sep = "")
  c0 <- fit_const(pred, y, start_const)
  grid <- as.numeric(seq(as.Date("2026-06-18"), as.Date("2026-08-28"), by = 5) - CFG$t0)
  prof <- list(); for (tc in grid) prof[[as.character(tc)]] <- fit_change(pred, y, tc, c0$par4)
  best <- as.numeric(names(prof))[which.min(sapply(prof, `[[`, "nll"))]
  for (tc in setdiff((best - 4):(best + 4), grid)) prof[[as.character(tc)]] <- fit_change(pred, y, tc, prof[[as.character(best)]]$par4)
  prof <- prof[order(as.numeric(names(prof)))]
  nll_tc <- sapply(prof, `[[`, "nll"); tcs <- as.numeric(names(prof)); best <- tcs[which.min(nll_tc)]
  b <- prof[[as.character(best)]]; ok <- tcs[nll_tc - min(nll_tc) <= 1.92]
  aic0 <- 2 * c0$nll + 6; aic1 <- 2 * b$nll + 10
  cat(sprintf("constant: nll %.2f r %.4f | changepoint %s: nll %.2f r %.4f -> %.4f | dAIC %.1f | tc 95%% %s to %s\n", c0$nll,
              gr(exp(c0$par4[2])), format(CFG$t0 + best), b$nll, gr(exp(b$par4[2])), gr(exp(b$par4[3])), aic0 - aic1,
              format(CFG$t0 + min(ok)), format(CFG$t0 + max(ok))))
  list(const = c0, prof = prof, tcs = tcs, nll_tc = nll_tc, best = b, aic0 = aic0, aic1 = aic1, ok = ok)
}
A <- analyse("Fit A: onset curve (reporting-delay corrected)", predA, dA$confirmed_total, c(log(25), log(0.245), log(20)))
B <- analyse("Fit B: confirmed cases by report date", predB, dB$y, c(log(35), log(0.245), log(30)))

## ---- uncertainty draws: change date ~ profile-likelihood weights; other parameters ~ MVN(Hessian) at that change date
draws_for <- function(F, pred, fit_label, nd = 1000, nband = 200) {
  w <- exp(-(F$nll_tc - min(F$nll_tc))); w <- w / sum(w)
  out <- vector("list", nd); bands <- list()
  for (i in 1:nd) {
    tc <- sample(F$tcs, 1, prob = w); fit <- F$prof[[as.character(tc)]]
    p4 <- mvn_draw(fit$par4, fit$hess); p4[2:3] <- pmin(pmax(p4[2:3], log(0.05)), log(1.5))
    inc <- tryCatch(sim_fit(exp(p4[1]), exp(p4[2]), exp(p4[3]), tc), error = function(e) NULL)
    if (is.null(inc)) next
    out[[i]] <- data.frame(fit = fit_label, tc_day = tc, logi0 = p4[1], beta1 = exp(p4[2]), beta2 = exp(p4[3]), k = exp(p4[4]),
      r1 = gr(exp(p4[2])), r2 = gr(exp(p4[3])), Re1 = Re_value(exp(p4[2]), CFG$gamma0, CFG$tau0, CFG$eps_base),
      Re2 = Re_value(exp(p4[3]), CFG$gamma0, CFG$tau0, CFG$eps_base), inc_end = inc[TEND])
    if (length(bands) < nband) bands[[length(bands) + 1]] <- list(inc = inc, mu = pred(inc))
  }
  list(draws = do.call(rbind, out), bands = bands)
}
cat("\nSampling parameter uncertainty ...\n")
DA <- draws_for(A, predA, "A"); DB <- draws_for(B, predB, "B")
write.csv(DA$draws, file.path(P$res, "fit_draws_A.csv"), row.names = FALSE); write.csv(DB$draws, file.path(P$res, "fit_draws_B.csv"), row.names = FALSE)

## ---- adequacy: 95% negative-binomial predictive coverage and residual autocorrelation (at the best fit)
adequacy <- function(F, pred, y) {
  b <- F$best; inc <- sim_fit(exp(b$par4[1]), exp(b$par4[2]), exp(b$par4[3]), b$tc); mu <- pred(inc); k <- exp(b$par4[4])
  lo <- qnbinom(.025, mu = mu, size = k); hi <- qnbinom(.975, mu = mu, size = k)
  pr_res <- (y - mu) / sqrt(mu + mu^2 / k)
  c(coverage95 = mean(y >= lo & y <= hi), resid_lag1_acf = acf(pr_res, plot = FALSE, lag.max = 1)$acf[2], resid_sd = sd(pr_res), k = k)
}
adq <- rbind(A = adequacy(A, predA, dA$confirmed_total), B = adequacy(B, predB, dB$y)); print(round(adq, 3))

## ---- sensitivity of the fitted growth rate to the fixed natural-history values (Fit B, change date fixed at best)
nh <- expand.grid(latent = c(4.0, 10.1), infectious = c(5, 10)); nh_rows <- list()
for (j in 1:nrow(nh)) {
  s <- 1 / nh$latent[j]; g <- 1 / nh$infectious[j]
  f <- function(p) nll_gen(p, B$best$tc, predB, dB$y, sigma = s, gamma = g)
  o <- nm2(f, B$best$par4, 2500)
  nh_rows[[j]] <- data.frame(latent = nh$latent[j], infectious = nh$infectious[j], nll = o$value,
    beta2 = exp(o$par[3]), r_after = growth_rate(exp(o$par[3]), s, g, CFG$tau0, CFG$eps_base, 1), r_before = growth_rate(exp(o$par[2]), s, g, CFG$tau0, CFG$eps_base, 1),
    inc_end = sim_fit(exp(o$par[1]), exp(o$par[2]), exp(o$par[3]), B$best$tc, sigma = s, gamma = g)[TEND])
}
nhtab <- do.call(rbind, nh_rows); print(round(nhtab, 4))
write.csv(nhtab, file.path(P$tab, "fit_naturalhistory_sensitivity.csv"), row.names = FALSE)

## ---- beta at each isolation effectiveness (refit at the best change date)
eps_rows <- list()
for (e in CFG$eps) for (nm in c("A", "B")) {
  F <- if (nm == "A") A else B; pred <- if (nm == "A") predA else predB; y <- if (nm == "A") dA$confirmed_total else dB$y
  fit <- fit_change(pred, y, F$best$tc, F$best$par4, eps = e)
  eps_rows[[length(eps_rows) + 1]] <- data.frame(fit = nm, eps = e, nll = fit$nll, beta_before = exp(fit$par4[2]), beta_after = exp(fit$par4[3]),
    r_before = gr(exp(fit$par4[2]), e), r_after = gr(exp(fit$par4[3]), e),
    Re_before = Re_value(exp(fit$par4[2]), CFG$gamma0, CFG$tau0, e), Re_after = Re_value(exp(fit$par4[3]), CFG$gamma0, CFG$tau0, e))
}
epstab <- do.call(rbind, eps_rows); write.csv(epstab, file.path(P$tab, "fit_by_effectiveness.csv"), row.names = FALSE); print(round(epstab[, -1], 4))

saveRDS(list(A = A, B = B, DA = DA, DB = DB, adq = adq, nhtab = nhtab, epstab = epstab, TEND = TEND), file.path(P$res, "fit.rds"))
cat(sprintf("\nFit finished in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 03 | Nine scenarios x paired draws; policy surface
######################################################################################################################
banner("03 | Nine scenarios x paired draws; policy surface")
## 03_scenarios.R -- nine intervention scenarios (delay 7/5/3 d  x  isolation effectiveness 50/80/100%) from 12 Sep 2026.
## Uncertainty: each of the n draws carries (growth rate r, starting incidence) from the fitted posterior (Fit B, joint with the change date)
## plus latent and infectious periods drawn from triangular distributions. Beta is re-calibrated per draw and per effectiveness so that the
## 7-day baseline reproduces that draw's fitted growth rate; scenarios then change ONLY the isolation rate. All scenarios share the same draws
## (paired), so cases averted are differences within a draw.
## Also computes a continuous policy surface (delay 2--10 d x effectiveness 30--100%).
## Outputs: results/scenarios.rds, results/tables/scenario_summary.csv, results/policy_surface.csv
set.seed(20261101); t_start <- Sys.time()
fit <- readRDS(file.path(P$res, "fit.rds")); D <- fit$DB$draws; n <- nrow(D)
L <- rtri(n, CFG$lat_tri); Di <- rtri(n, CFG$inf_tri); sigma <- 1 / L; gamma <- 1 / Di
i0 <- D$inc_end; r2 <- D$r2

runs <- list(); traj <- list(); k <- 0
for (eps in CFG$eps) {
  beta <- vapply(1:n, function(i) calibrate_beta(r2[i], sigma[i], gamma[i], CFG$tau0, eps, 1), 0)
  for (dl in CFG$delays) {
    tr <- matrix(NA_real_, n, CFG$horizon); out <- vector("list", n)
    for (i in 1:n) {
      s <- run_det(beta[i], sigma[i], gamma[i], 1 / dl, eps, i0[i], tend = CFG$max_days)
      o <- outcomes(s); tr[i, ] <- s$inc[1:CFG$horizon]
      out[[i]] <- data.frame(draw = i, eps = eps, delay = dl, cum12 = o$cum, duration = o$duration, peak = o$peak, tpeak = o$tpeak,
                             Re = Re_value(beta[i], gamma[i], 1 / dl, eps, 1), frac_iso = frac_isolated(gamma[i], 1 / dl, 1), beta = beta[i])
    }
    k <- k + 1; runs[[k]] <- do.call(rbind, out); traj[[sprintf("e%.1f_d%d", eps, dl)]] <- tr
    cat(sprintf("eps %.1f delay %d done (%.1f min)\n", eps, dl, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}
R <- do.call(rbind, runs)
base <- R[R$delay == 7, c("draw", "eps", "cum12")]; names(base)[3] <- "cum12_base"
R <- merge(R, base, by = c("draw", "eps")); R$averted <- R$cum12_base - R$cum12; R$pct_averted <- 100 * R$averted / R$cum12_base
R$deaths_averted <- R$averted * CFG$cfr

## ---- summary table
summ <- do.call(rbind, lapply(split(R, list(R$eps, R$delay)), function(d) {
  qa <- q3(d$averted); qp <- q3(d$pct_averted); qc <- q3(d$cum12); qr <- q3(d$Re); qd <- q3(d$deaths_averted)
  dur <- d$duration
  data.frame(eps = d$eps[1], delay = d$delay[1], cum12_med = qc[1], cum12_lo = qc[2], cum12_hi = qc[3],
    averted_med = qa[1], averted_lo = qa[2], averted_hi = qa[3], pct_med = qp[1], pct_lo = qp[2], pct_hi = qp[3],
    deaths_med = qd[1], deaths_lo = qd[2], deaths_hi = qd[3],
    Re_med = qr[1], Re_lo = qr[2], Re_hi = qr[3], p_Re_below1 = mean(d$Re < 1),
    dur_med = quantile(dur, .5, type = 1), dur_lo = quantile(dur, .025, type = 1), dur_hi = quantile(dur, .975, type = 1),
    p_end_12m = mean(dur <= 365), p_end_24m = mean(dur <= 730), p_end_60m = mean(is.finite(dur)),
    peak_med = median(d$peak), tpeak_med = median(d$tpeak), frac_iso_med = median(d$frac_iso))
}))
rownames(summ) <- NULL; summ <- summ[order(summ$eps, -summ$delay), ]
write.csv(summ, file.path(P$tab, "scenario_summary.csv"), row.names = FALSE)
print(round(summ[, c("eps", "delay", "cum12_med", "pct_med", "pct_lo", "pct_hi", "Re_med", "p_Re_below1", "p_end_12m", "p_end_24m", "p_end_60m")], 2), row.names = FALSE)
saveRDS(list(R = R, traj = traj, L = L, Di = Di, summ = summ, n = n), file.path(P$res, "scenarios.rds"))

## ---- continuous policy surface (first 100 draws; 12-month cases)
nd <- 100; dgrid <- 2:10; egrid <- seq(0.3, 1, by = 0.1); surf <- list(); kk <- 0
for (eps in egrid) {
  beta <- vapply(1:nd, function(i) calibrate_beta(r2[i], sigma[i], gamma[i], CFG$tau0, eps, 1), 0)
  c7 <- vapply(1:nd, function(i) outcomes(run_det(beta[i], sigma[i], gamma[i], CFG$tau0, eps, i0[i], tend = CFG$horizon))$cum, 0)
  for (dl in dgrid) {
    cs <- vapply(1:nd, function(i) outcomes(run_det(beta[i], sigma[i], gamma[i], 1 / dl, eps, i0[i], tend = CFG$horizon))$cum, 0)
    re <- vapply(1:nd, function(i) Re_value(beta[i], gamma[i], 1 / dl, eps, 1), 0)
    kk <- kk + 1; surf[[kk]] <- data.frame(eps = eps, delay = dl, pct_averted = median(100 * (c7 - cs) / c7), Re = median(re), p_Re_below1 = mean(re < 1))
  }
}
surf <- do.call(rbind, surf); write.csv(surf, file.path(P$res, "policy_surface.csv"), row.names = FALSE)
cat(sprintf("Scenarios finished in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 04 | Stochastic check
######################################################################################################################
banner("04 | Stochastic check")
## 04_stochastic.R -- stochastic (binomial-chain) version of the model: do the conclusions survive chance events?
## Same posterior draws as 03_scenarios.R (first n = 500). Transmission is calibrated in discrete time so the stochastic model has the
## same expected growth rate as the deterministic one. Outputs: results/stochastic.rds, results/tables/stochastic_summary.csv
set.seed(20261102); t_start <- Sys.time()
sc <- readRDS(file.path(P$res, "scenarios.rds")); fit <- readRDS(file.path(P$res, "fit.rds"))
n <- 500; DT <- 0.25; D <- fit$DB$draws[1:n, ]; sigma <- 1 / sc$L[1:n]; gamma <- 1 / sc$Di[1:n]

rows <- list()
for (eps in CFG$eps) {
  beta <- vapply(1:n, function(i) calibrate_beta_disc(D$r2[i], sigma[i], gamma[i], CFG$tau0, eps, DT), 0)
  v <- t(vapply(1:n, function(i) init_disc(beta[i], sigma[i], gamma[i], CFG$tau0, eps, D$inc_end[i], DT), numeric(3)))
  for (dl in CFG$delays) {
    s <- run_stoch(beta, sigma, gamma, 1 / dl, eps, v[, 1], v[, 2], v[, 3], dt = DT)
    d <- sc$R[sc$R$draw <= n & sc$R$eps == eps & sc$R$delay == dl, ]; d <- d[order(d$draw), ]
    rows[[length(rows) + 1]] <- data.frame(eps = eps, delay = dl, draw = 1:n, cum12 = s$cum12, ext = s$ext, det_cum12 = d$cum12, det_dur = d$duration)
    cat(sprintf("eps %.1f delay %d done (%.1f min)\n", eps, dl, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}
S <- do.call(rbind, rows); saveRDS(S, file.path(P$res, "stochastic.rds"))
ss <- do.call(rbind, lapply(split(S, list(S$eps, S$delay)), function(d) {
  q <- q3(d$cum12); qd <- q3(d$det_cum12)
  data.frame(eps = d$eps[1], delay = d$delay[1], stoch_med = q[1], stoch_lo = q[2], stoch_hi = q[3], det_med = qd[1],
             p_end_12m = mean(d$ext <= 365, na.rm = TRUE), p_end_24m = mean(d$ext <= 730, na.rm = TRUE), p_end_60m = mean(!is.na(d$ext)),
             ext_med = median(d$ext, na.rm = TRUE), det_p_end_60m = mean(is.finite(d$det_dur)))
}))
ss <- ss[order(ss$eps, -ss$delay), ]; rownames(ss) <- NULL
ss$pct_red_med <- NA; for (e in CFG$eps) { b <- ss$stoch_med[ss$eps == e & ss$delay == 7]; i <- ss$eps == e; ss$pct_red_med[i] <- 100 * (b - ss$stoch_med[i]) / b }
write.csv(ss, file.path(P$tab, "stochastic_summary.csv"), row.names = FALSE); print(round(ss, 2), row.names = FALSE)
cat(sprintf("Stochastic finished in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 05 | Global (PRCC) and structural sensitivity analyses
######################################################################################################################
banner("05 | Global (PRCC) and structural sensitivity analyses")
## 05_sensitivity.R -- (a) global sensitivity (partial rank correlation) and (b) one-at-a-time / structural variants.
## Outputs: results/tables/prcc.csv, results/tables/variants.csv
set.seed(20261103); t_start <- Sys.time()
sc <- readRDS(file.path(P$res, "scenarios.rds")); fit <- readRDS(file.path(P$res, "fit.rds"))
DB <- fit$DB$draws; DA <- fit$DA$draws; eps_b <- CFG$eps_base

## ---------------------------------------------------------------- (a) PRCC on 12-month cases
ns <- 1500; ix <- sample(nrow(DB), ns, replace = TRUE)
X <- data.frame(latent_period = rtri(ns, CFG$lat_tri), infectious_period = rtri(ns, CFG$inf_tri),
                growth_rate = DB$r2[ix], starting_incidence = DB$inc_end[ix],
                isolation_effectiveness = runif(ns, 0.5, 1), isolation_delay = runif(ns, 3, 7))
y <- vapply(1:ns, function(i) {
  s <- 1 / X$latent_period[i]; g <- 1 / X$infectious_period[i]; e <- X$isolation_effectiveness[i]
  b <- calibrate_beta(X$growth_rate[i], s, g, CFG$tau0, e, 1)
  outcomes(run_det(b, s, g, 1 / X$isolation_delay[i], e, X$starting_incidence[i], tend = CFG$horizon))$cum }, 0)
prcc <- function(X, y) {
  Rk <- apply(X, 2, rank); yr <- rank(y)
  vapply(seq_len(ncol(X)), function(j) {
    A <- cbind(1, Rk[, -j, drop = FALSE]); rx <- Rk[, j] - A %*% qr.solve(A, Rk[, j]); ry <- yr - A %*% qr.solve(A, yr); cor(rx, ry) }, 0)
}
pe <- prcc(X, y); B <- 300
bs <- replicate(B, { i <- sample(ns, replace = TRUE); prcc(X[i, ], y[i]) })
pt <- data.frame(parameter = names(X), prcc = pe, lo = apply(bs, 1, quantile, .025), hi = apply(bs, 1, quantile, .975))
pt <- pt[order(-abs(pt$prcc)), ]; write.csv(pt, file.path(P$tab, "prcc.csv"), row.names = FALSE); print(round(pt[, -1], 3))
cat(sprintf("PRCC done (%.1f min)\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

## ---------------------------------------------------------------- (b) variants (isolation effectiveness 80%)
nd <- 300; idx <- 1:nd
variants <- list(
  list(name = "Reference (all uncertainty propagated)", group = "Reference"),
  list(name = "Latent period 4.0 d", group = "Natural history", L = 4.0),
  list(name = "Latent period 10.1 d", group = "Natural history", L = 10.1),
  list(name = "Infectious period 5 d", group = "Natural history", D = 5),
  list(name = "Infectious period 10 d", group = "Natural history", D = 10),
  list(name = "Growth rate at 2.5th percentile", group = "Calibration", r = quantile(DB$r2, .025)),
  list(name = "Growth rate at 97.5th percentile", group = "Calibration", r = quantile(DB$r2, .975)),
  list(name = "Starting incidence at 2.5th percentile", group = "Calibration", i0 = quantile(DB$inc_end, .025)),
  list(name = "Starting incidence at 97.5th percentile", group = "Calibration", i0 = quantile(DB$inc_end, .975)),
  list(name = "Calibrated to onset curve (Fit A) instead of report counts", group = "Calibration", src = "A"),
  list(name = "True baseline delay 4.6 d (onset-to-sample, preprint)", group = "Baseline delay", tau0 = 1 / 4.6),
  list(name = "True baseline delay 8.6 d (onset-to-result, preprint)", group = "Baseline delay", tau0 = 1 / 8.6),
  list(name = "Erlang-3 delay (less variable)", group = "Structure", m = 3),
  list(name = "Erlang-10 delay (nearly fixed delay)", group = "Structure", m = 10),
  list(name = "Delay = mean time to isolation among those isolated", group = "Structure", cond = TRUE),
  list(name = "Population at risk 10 million", group = "Structure", pop = 1e7),
  list(name = "Interventions ramp up linearly over 30 d", group = "Implementation", ramp = 30),
  list(name = "Interventions ramp up linearly over 90 d", group = "Implementation", ramp = 90))
g <- function(v, key, default = NULL) { if (key %in% names(v)) v[[key]] else default }   # exact lookup (no $ partial matching)
run_variant <- function(v) {
  src <- if (!is.null(g(v, "src"))) DA else DB
  L <- if (!is.null(g(v, "L"))) rep(v[["L"]], nd) else sc$L[idx]; Dd <- if (!is.null(g(v, "D"))) rep(v[["D"]], nd) else sc$Di[idx]
  sigma <- 1 / L; gamma <- 1 / Dd
  r <- if (!is.null(g(v, "r"))) rep(v[["r"]], nd) else src$r2[idx]; i0 <- if (!is.null(g(v, "i0"))) rep(v[["i0"]], nd) else src$inc_end[idx]
  m <- g(v, "m", 1); ramp <- g(v, "ramp", 0); pop <- g(v, "pop", CFG$pop)
  tau0 <- rep(g(v, "tau0", CFG$tau0), nd)
  cond <- isTRUE(g(v, "cond", FALSE))
  if (cond) tau0 <- pmax(1 / 7 - gamma, 1e-6)           # delay defined among those isolated: tau = 1/D - gamma
  tauf <- function(dl, i) if (cond) max(1 / dl - gamma[i], 1e-6) else 1 / dl
  res <- sapply(1:nd, function(i) {
    b <- calibrate_beta(r[i], sigma[i], gamma[i], tau0[i], eps_b, m)
    base <- outcomes(run_det(b, sigma[i], gamma[i], tau0[i], eps_b, i0[i], tau0 = tau0[i], m = m, tend = CFG$horizon, pop = pop))$cum
    o5 <- outcomes(run_det(b, sigma[i], gamma[i], tauf(5, i), eps_b, i0[i], tau0 = tau0[i], m = m, tend = CFG$horizon, ramp = ramp, pop = pop))$cum
    o3 <- outcomes(run_det(b, sigma[i], gamma[i], tauf(3, i), eps_b, i0[i], tau0 = tau0[i], m = m, tend = CFG$horizon, ramp = ramp, pop = pop))$cum
    c(100 * (base - o5) / base, 100 * (base - o3) / base, Re_value(b, gamma[i], tauf(3, i), eps_b, m), frac_isolated(gamma[i], tau0[i], m), base)
  })
  q5 <- q3(res[1, ]); q3_ <- q3(res[2, ])
  data.frame(group = v$group, variant = v$name, pct5 = q5[1], pct5_lo = q5[2], pct5_hi = q5[3], pct3 = q3_[1], pct3_lo = q3_[2], pct3_hi = q3_[3],
             Re3 = median(res[3, ]), frac_iso_base = median(res[4, ]), base_cases = median(res[5, ]))
}
V <- do.call(rbind, lapply(variants, function(v) { cat(" variant:", v$name, "\n"); run_variant(v) }))
write.csv(V, file.path(P$tab, "variants.csv"), row.names = FALSE); print(round(V[, -(1:2)], 1), row.names = FALSE)
cat(sprintf("Sensitivity finished in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 06 | Out-of-sample validation
######################################################################################################################
banner("06 | Out-of-sample validation")
## 06_validation.R -- out-of-sample check. The model is fitted to data up to 12 Sep 2026 ONLY. We forecast confirmed cases reported up to
## 27 Sep and 28 Oct 2026 (baseline: nothing changes) and compare with later official cumulative counts (UKHSA 27 Sep; CDC 28 Oct), which
## were never used in fitting. Benchmark: constant incidence equal to the mean of the last 14 days.
## Outputs: results/validation.rds, results/tables/validation.csv
set.seed(20261104)
pr <- readRDS(file.path(P$res, "prep.rds")); fit <- readRDS(file.path(P$res, "fit.rds")); DB <- fit$DB$draws
ext <- read.csv(file.path(P$data, "cumulative_confirmed_drc.csv"), stringsAsFactors = FALSE); ext$date <- as.Date(ext$date)
TEND <- fit$TEND; HOR <- 46; TT <- TEND + HOR; WK <- pr$delay$shape; WS <- pr$delay$scale
cc <- function(L) pweibull(L + 0.5, WK, WS); pm <- diff(c(0, cc(0:(TT - 1)))); PM <- matrix(0, TT, TT); for (u in 1:TT) PM[u, 1:u] <- pm[u:1]

cum0 <- pr$last$confirmed_case_history                                    # 7,200 at 12 Sep
obs <- ext[ext$date > CFG$cutoff & ext$date <= as.Date("2026-10-28"), ]
obs$u <- as.numeric(obs$date - CFG$t0); obs$observed <- obs$cum_confirmed_drc - cum0; obs$ndays <- as.numeric(obs$date - CFG$cutoff)

nd <- 500; ix <- sample(nrow(DB), nd, replace = TRUE); E <- matrix(NA, nd, HOR)           # expected reported cases per day after cutoff
sg <- CFG$sigma0; gm <- CFG$gamma0; rhs <- make_rhs(1, sg, gm, CFG$pop)
for (j in 1:nd) {
  d <- DB[ix[j], ]; y0 <- init_state(d$beta1, sg, gm, CFG$tau0, CFG$eps_base, exp(d$logi0), 1)
  p <- function(b) c(b, CFG$tau0, CFG$eps_base, CFG$tau0, 0)
  o1 <- ode(y0, 0:d$tc_day, rhs, p(d$beta1), method = "lsoda", rtol = 1e-6, atol = 1e-6); C <- o1[, 7]
  o2 <- ode(o1[nrow(o1), -1], d$tc_day:TT, rhs, p(d$beta2), method = "lsoda", rtol = 1e-6, atol = 1e-6); C <- c(C, o2[-1, 7])
  E[j, ] <- as.vector(PM %*% diff(C))[(TEND + 1):TT]
}
fc <- function(days) { cs <- rowSums(E[, 1:days, drop = FALSE]); k <- DB$k[ix]
  pred <- rnbinom(nd, mu = cs, size = k * days)                                 # observation noise: independent NB days
  list(mu = quantile(cs, c(.5, .025, .975)), pred = quantile(pred, c(.5, .025, .975))) }
wide <- read.csv(file.path(P$bvd, "bvd_2026_streams_wide.csv")); wide$date <- as.Date(wide$date)
last14 <- (cum0 - wide$confirmed_case_history[wide$date == CFG$cutoff - 14]) / 14       # mean new confirmed/day over the last 14 days
rows <- list(); for (i in seq_len(nrow(obs))) {
  f <- fc(obs$ndays[i]); nv <- last14 * obs$ndays[i]
  rows[[i]] <- data.frame(date = obs$date[i], days_ahead = obs$ndays[i], observed = obs$observed[i], model_median = f$mu[1], model_lo = f$mu[2], model_hi = f$mu[3],
    pred_lo = f$pred[2], pred_hi = f$pred[3], inside_param_interval = obs$observed[i] >= f$mu[2] & obs$observed[i] <= f$mu[3],
    inside_pred_interval = obs$observed[i] >= f$pred[2] & obs$observed[i] <= f$pred[3], naive = nv,
    err_model_pct = 100 * (f$mu[1] - obs$observed[i]) / obs$observed[i], err_naive_pct = 100 * (nv - obs$observed[i]) / obs$observed[i],
    cum_model = cum0 + f$mu[1], cum_observed = cum0 + obs$observed[i]) }
V <- do.call(rbind, rows); write.csv(V, file.path(P$tab, "validation.csv"), row.names = FALSE); print(round(V[, -1], 1), row.names = FALSE)
saveRDS(list(V = V, E = E, obs = obs, cum0 = cum0), file.path(P$res, "validation.rds"))
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 07 | Figures, LaTeX tables and in-text number macros
######################################################################################################################
banner("07 | Figures, LaTeX tables and in-text number macros")
## 07_figures_tables.R -- every figure, table and in-text number of the manuscript is produced here from results/*.
## Outputs: figures/fig*.pdf, results/tables/*.tex, manuscript/numbers.tex (LaTeX macros used in the text)
pr <- readRDS(file.path(P$res, "prep.rds")); fit <- readRDS(file.path(P$res, "fit.rds")); sc <- readRDS(file.path(P$res, "scenarios.rds"))
st <- readRDS(file.path(P$res, "stochastic.rds")); va <- readRDS(file.path(P$res, "validation.rds"))
summ <- sc$summ; R <- sc$R; traj <- sc$traj
prcc <- read.csv(file.path(P$tab, "prcc.csv")); vtab <- read.csv(file.path(P$tab, "variants.csv"))
ss <- read.csv(file.path(P$tab, "stochastic_summary.csv")); surf <- read.csv(file.path(P$res, "policy_surface.csv"))
V <- va$V; nh <- fit$nhtab; eps_tab <- fit$epstab
COL <- c("7" = "#4d4d4d", "5" = "#e08a00", "3" = "#0b6e99"); ORG <- "#e08a00"; BLU <- "#0b6e99"
LTY <- c("0.5" = 3, "0.8" = 1, "1" = 2)
bandq <- function(M) apply(M, 2, quantile, c(.025, .5, .975), na.rm = TRUE)
numw <- c("0.5" = "Fifty", "0.8" = "Eighty", "1" = "Hundred", "7" = "Seven", "5" = "Five", "3" = "Three")
em <- function(e) numw[as.character(e)]; dm <- function(d) numw[as.character(d)]
d0 <- function(x, d = 0) gsub("^-", "$-$", formatC(x, format = "f", digits = d, big.mark = ","))
ci <- function(m, lo, hi, d = 0, sep = " to ") sprintf("%s (%s%s%s)", d0(m, d), d0(lo, d), sep, d0(hi, d))
pcs <- function(x, d = 0) paste0(d0(x, d), "\\%")
dstr <- function(date) format(date, "%-d %B")
months <- function(days) days / 30.44

## =================================================================== Figure 2: fit
pdf(file.path(P$fig, "fig2_fit.pdf"), width = 7.4, height = 6.4)
par(mfrow = c(2, 2), mar = c(3.7, 4, 2.1, 0.8), cex = 0.68, mgp = c(2.3, 0.7, 0))
WK <- pr$delay$shape; WS <- pr$delay$scale; cc <- function(L) pweibull(L + 0.5, WK, WS)
plot(pr$agg$lag, pr$agg$ratio, pch = 16, ylim = c(0, 1.1), xlab = "Days between symptom onset and report", ylab = "Share of final count already reported",
     main = "(a) Reporting delay (55 report vintages)"); lines(0:40, cc(0:40), col = BLU, lwd = 2); abline(h = 1, lty = 3)
legend("bottomright", bty = "n", lwd = 2, col = BLU, legend = sprintf("Weibull fit, mean %.1f d", pr$delay$mean))
dA <- pr$dA; dB <- pr$dB; xA <- dA$onset_date; bA <- bandq(do.call(rbind, lapply(fit$DA$bands, `[[`, "mu"))); nA <- bandq(do.call(rbind, lapply(fit$DA$bands, function(z) z$inc[dA$k])))
plot(xA, dA$confirmed_total, type = "n", ylim = c(0, max(dA$confirmed_total, bA) * 1.08), xlab = "Symptom-onset date", ylab = "Confirmed cases per day", main = "(b) Fit A: cases by onset date")
polygon(c(xA, rev(xA)), c(bA[1, ], rev(bA[3, ])), col = adjustcolor(ORG, .25), border = NA); lines(xA, nA[2, ], col = BLU, lwd = 1.6, lty = 2); lines(xA, bA[2, ], col = ORG, lwd = 2)
points(xA, dA$confirmed_total, pch = ifelse(dA$lag >= 21, 16, 1), cex = 0.6)
legend("topleft", bty = "n", cex = 0.85, pch = c(16, 1, NA, NA), lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 1.6), col = c("black", "black", ORG, BLU),
       legend = c("observed", "observed, last 3 weeks (incomplete)", "model, as reported", "model, true onset incidence"))
bB <- bandq(do.call(rbind, lapply(fit$DB$bands, `[[`, "mu")) ) / rep(dB$ndays, each = 3) -> bBq
bB <- apply(do.call(rbind, lapply(fit$DB$bands, `[[`, "mu")), 2, quantile, c(.025, .5, .975)); bB <- sweep(bB, 2, dB$ndays, "/")
plot(dB$date, dB$y / dB$ndays, pch = 16, cex = 0.6, ylim = c(0, max(dB$y / dB$ndays, bB) * 1.05), xlab = "Report date", ylab = "Confirmed cases per day", main = "(c) Fit B: cases by report date")
polygon(c(dB$date, rev(dB$date)), c(bB[1, ], rev(bB[3, ])), col = adjustcolor(ORG, .25), border = NA); lines(dB$date, bB[2, ], col = ORG, lwd = 2); points(dB$date, dB$y / dB$ndays, pch = 16, cex = 0.6)
abline(v = CFG$t0 + median(fit$DB$draws$tc_day), lty = 3, col = "grey40")
A <- fit$A; B <- fit$B
plot(CFG$t0 + A$tcs, -(A$nll_tc - min(A$nll_tc)), type = "b", pch = 16, cex = 0.6, col = BLU, lwd = 1.5, ylim = c(-8, 0.3), xlab = "Date transmission changes", ylab = "Profile log-likelihood (vs best)", main = "(d) Timing of the change in transmission")
lines(CFG$t0 + B$tcs, -(B$nll_tc - min(B$nll_tc)), type = "b", pch = 17, cex = 0.6, col = ORG, lwd = 1.5); abline(h = -1.92, lty = 2)
legend("topright", bty = "n", pch = c(16, 17), col = c(BLU, ORG), lty = 1, legend = c("Fit A (onset)", "Fit B (report)")); dev.off()

## =================================================================== Figure 3: incidence trajectories
pdf(file.path(P$fig, "fig3_trajectories.pdf"), width = 7.4, height = 2.9)
par(mfrow = c(1, 3), mar = c(3.8, 3.9, 2, 0.5), cex = 0.68, mgp = c(2.3, 0.7, 0), oma = c(0, 0, 0, 0)); xs <- CFG$cutoff + 1:CFG$horizon
for (e in CFG$eps) {
  plot(xs, rep(1, CFG$horizon), type = "n", log = "y", ylim = c(0.1, 600), yaxt = "n", xlab = "Date", ylab = if (e == 0.5) "New symptomatic cases per day" else "",
       main = sprintf("Isolation %.0f%% effective", 100 * e)); axis(2, c(0.1, 1, 10, 100), c("0.1", "1", "10", "100"), las = 1)
  for (dl in CFG$delays) { q <- bandq(traj[[sprintf("e%.1f_d%d", e, dl)]]); polygon(c(xs, rev(xs)), c(pmax(q[1, ], 0.05), rev(pmax(q[3, ], 0.05))), col = adjustcolor(COL[as.character(dl)], .2), border = NA)
    lines(xs, pmax(q[2, ], 0.05), col = COL[as.character(dl)], lwd = 2) }
  abline(h = 1, lty = 3, col = "grey50")
  if (e == 0.5) legend("bottomleft", bty = "n", lwd = 2, col = COL, title = "Delay to isolation", legend = c("7 days", "5 days", "3 days"), cex = 0.9) }
dev.off()

## =================================================================== Figure 4: cases averted + time to end
pdf(file.path(P$fig, "fig4_averted_duration.pdf"), width = 7.4, height = 3.3)
par(mfrow = c(1, 2), mar = c(3.8, 4, 2, 0.8), cex = 0.68, mgp = c(2.3, 0.7, 0))
xb <- barplot(rbind(summ$pct_med[summ$delay == 5], summ$pct_med[summ$delay == 3]), beside = TRUE, col = COL[c("5", "3")], ylim = c(0, 112), names.arg = c("50%", "80%", "100%"),
              xlab = "Isolation effectiveness", ylab = "Cases averted in 12 months (%)", main = "(a) Reduction relative to a 7-day delay")
arrows(xb, rbind(summ$pct_lo[summ$delay == 5], summ$pct_lo[summ$delay == 3]), xb, rbind(summ$pct_hi[summ$delay == 5], summ$pct_hi[summ$delay == 3]), angle = 90, code = 3, length = 0.03)
legend("topleft", fill = COL[c("5", "3")], bty = "n", legend = c("5 days", "3 days"), horiz = TRUE, cex = 0.9)
par(mar = c(3.8, 4, 2, 5.6)); plot(NA, xlim = c(0, 60), ylim = c(0, 1), xlab = "Months from 12 September 2026", ylab = "Probability outbreak has ended", main = "(b) Time to end of transmission")
for (e in CFG$eps) for (dl in CFG$delays) { d <- R$duration[R$eps == e & R$delay == dl] / 30.44; xx <- seq(0, 60, by = 0.25); lines(xx, vapply(xx, function(x) mean(d <= x), 0), col = COL[as.character(dl)], lty = LTY[as.character(e)], lwd = 1.6) }
legend("right", inset = c(-0.27, 0), xpd = NA, bty = "n", cex = 0.85, lwd = 1.6, col = c(COL, "black", "black", "black"), lty = c(1, 1, 1, LTY), legend = c("7 d", "5 d", "3 d", "50%", "80%", "100%")); dev.off()

## =================================================================== Figure 5: policy surface
pdf(file.path(P$fig, "fig5_policy_surface.pdf"), width = 5.6, height = 4.2)
par(mar = c(4, 4, 2, 1), cex = 0.75, mgp = c(2.3, 0.7, 0))
dg <- sort(unique(surf$delay)); eg <- sort(unique(surf$eps))
Z <- matrix(surf$pct_averted, length(dg), length(eg), dimnames = list(dg, eg), byrow = FALSE)[, ]
Z <- sapply(eg, function(e) surf$pct_averted[surf$eps == e][order(surf$delay[surf$eps == e])]); Zp <- sapply(eg, function(e) surf$p_Re_below1[surf$eps == e][order(surf$delay[surf$eps == e])])
image(dg, eg, pmax(pmin(Z, 100), -100), col = hcl.colors(40, "Blue-Red 3", rev = TRUE), zlim = c(-100, 100), xlab = "Mean onset-to-isolation delay (days)", ylab = "Isolation effectiveness",
      main = "Cases averted in 12 months vs a 7-day delay (%)")
contour(dg, eg, Z, levels = c(25, 50, 75, 90), add = TRUE, labcex = 0.8, col = "black"); contour(dg, eg, Zp, levels = 0.5, add = TRUE, col = "white", lwd = 2.5, drawlabels = FALSE)
points(rep(c(7, 5, 3), 3), rep(CFG$eps, each = 3), pch = 21, bg = "white", cex = 1.1)
mtext("blue: fewer cases than at 7 days; red: more cases", side = 1, line = 3.1, cex = 0.7)
legend("bottomleft", bty = "n", lwd = 2.5, col = "white", legend = expression(paste("P(", R[e] < 1, ") = 0.5")), text.col = "white", cex = 0.9); dev.off()

## =================================================================== Figure 6: sensitivity
pdf(file.path(P$fig, "fig6_sensitivity.pdf"), width = 7.4, height = 5.6)
layout(matrix(1:2, 1), widths = c(1, 1.55)); par(cex = 0.66, mgp = c(2.2, 0.7, 0))
lab <- c(latent_period = "Latent period", infectious_period = "Infectious period", growth_rate = "Baseline growth rate", starting_incidence = "Starting incidence",
         isolation_effectiveness = "Isolation effectiveness", isolation_delay = "Isolation delay"); pp <- prcc[order(abs(prcc$prcc)), ]
par(mar = c(4, 9.5, 2, 0.6)); plot(pp$prcc, seq_len(nrow(pp)), xlim = c(-1, 1), pch = 16, yaxt = "n", ylab = "", xlab = "PRCC with 12-month cases", ylim = c(0.5, nrow(pp) + 0.5), main = "(a) Global sensitivity")
segments(pp$lo, seq_len(nrow(pp)), pp$hi, seq_len(nrow(pp))); axis(2, seq_len(nrow(pp)), lab[pp$parameter], las = 1); abline(v = 0, lty = 3)
short <- c("Reference", "Latent 4.0 d", "Latent 10.1 d", "Infectious 5 d", "Infectious 10 d", "Growth rate, 2.5th pct", "Growth rate, 97.5th pct", "Start incidence, 2.5th pct", "Start incidence, 97.5th pct",
           "Fitted to onset curve (A)", "Baseline delay 4.6 d", "Baseline delay 8.6 d", "Erlang-3 delay", "Erlang-10 delay", "Delay among isolated only", "Population 10 million", "Ramp-up 30 d", "Ramp-up 90 d")
par(mar = c(4, 10.5, 2, 0.8)); n <- nrow(vtab); y <- n:1
plot(NA, xlim = c(-110, 100), ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "", xlab = "Cases averted in 12 months (%), isolation 80% effective", main = "(b) Assumptions and structure")
abline(v = c(0, 50, 100), lty = 3, col = "grey70"); abline(h = y[1] - 0.5, col = "grey85")
segments(vtab$pct5_lo, y + 0.17, vtab$pct5_hi, y + 0.17, col = COL["5"]); points(vtab$pct5, y + 0.17, pch = 16, col = COL["5"])
segments(vtab$pct3_lo, y - 0.17, vtab$pct3_hi, y - 0.17, col = COL["3"]); points(vtab$pct3, y - 0.17, pch = 16, col = COL["3"])
axis(2, y, short, las = 1); legend("bottomleft", bty = "n", pch = 16, col = COL[c("5", "3")], legend = c("5-day delay", "3-day delay"), cex = 0.9); dev.off()

## =================================================================== Figure 7: stochastic vs deterministic
pdf(file.path(P$fig, "fig7_stochastic.pdf"), width = 7.4, height = 3.3)
par(mar = c(4.2, 4.3, 2, 0.8), cex = 0.68, mgp = c(2.3, 0.7, 0))
st$grp <- factor(sprintf("%d%%|%dd", round(100 * st$eps), st$delay), levels = unlist(lapply(CFG$eps, function(e) sprintf("%d%%|%dd", round(100 * e), CFG$delays))))
bp <- boxplot(cum12 ~ grp, st, log = "y", col = rep(COL[c("7", "5", "3")], 3), outline = FALSE, xaxt = "n", yaxt = "n", xlab = "Isolation effectiveness", ylab = "Symptomatic cases in 12 months",
              main = "Stochastic model (boxes: median, IQR, 95% range) vs deterministic median (diamonds)", range = 0)
axis(1, c(2, 5, 8), c("50%", "80%", "100%")); axis(2, c(1e3, 3e3, 1e4, 3e4, 1e5), c("1,000", "3,000", "10,000", "30,000", "100,000"), las = 1); abline(v = c(3.5, 6.5), lty = 3, col = "grey70")
dm_ <- aggregate(det_cum12 ~ grp, st, median); points(seq_len(nrow(dm_)), dm_$det_cum12[match(levels(st$grp), dm_$grp)], pch = 23, bg = "white", cex = 1.1)
legend("bottomleft", fill = COL[c("7", "5", "3")], bty = "n", legend = c("7 d", "5 d", "3 d"), horiz = TRUE, cex = 0.9); dev.off()

## =================================================================== Figure 8: out-of-sample validation
pdf(file.path(P$fig, "fig8_validation.pdf"), width = 5.4, height = 3.8)
par(mar = c(4, 4.2, 2, 0.8), cex = 0.72, mgp = c(2.3, 0.7, 0))
cs <- t(apply(va$E, 1, cumsum)); q <- apply(cs, 2, quantile, c(.025, .5, .975)); xd <- 1:ncol(cs)
plot(NA, xlim = c(0, 46), ylim = c(0, max(q, va$obs$observed) * 1.05), xlab = "Days after 12 September 2026", ylab = "Cumulative new confirmed cases", main = "Forecast fitted to data up to 12 Sep vs later counts")
polygon(c(xd, rev(xd)), c(q[1, ], rev(q[3, ])), col = adjustcolor(ORG, .25), border = NA); lines(xd, q[2, ], col = ORG, lwd = 2)
last14 <- V$naive[2] / V$days_ahead[2]; lines(c(0, 46), c(0, 46 * last14), lty = 2, col = "grey40")
points(c(0, va$obs$ndays), c(0, va$obs$observed), pch = 16, cex = 1.1); lines(c(0, va$obs$ndays), c(0, va$obs$observed), lwd = 1)
legend("topleft", bty = "n", cex = 0.85, lwd = c(2, 1, 1), lty = c(1, 2, 1), pch = c(NA, NA, 16), col = c(ORG, "grey40", "black"), legend = c("model median and 95% band", "constant last-14-day rate", "observed (official counts)")); dev.off()

## =================================================================== LaTeX tables
rowsep <- "\\\\"
wtab <- function(file, lines) writeLines(lines, file.path(P$tab, file))
# --- main scenario table
durtxt <- function(r) { if (!is.finite(r$dur_med)) return("$>$60")
  hi <- if (is.finite(r$dur_hi)) sprintf("%.0f", months(r$dur_hi)) else "$>$60"; sprintf("%.0f (%.0f to %s)", months(r$dur_med), months(r$dur_lo), hi) }
L <- c("\\begin{tabular}{@{}rrcrrrrr@{}}", "\\toprule",
  "Isolation & Delay & $R_e$ & Cases in 12 months & Cases averted (\\%) & Deaths averted$^\\dagger$ & Months to end & P(ended $\\le$24 mo) \\\\",
  "effectiveness & (days) & median (95\\% UI) & median (95\\% UI) & median (95\\% UI) & median & median (95\\% UI) &  \\\\", "\\midrule")
for (e in CFG$eps) { for (dl in CFG$delays) { r <- summ[summ$eps == e & summ$delay == dl, ]
  L <- c(L, sprintf("%s & %d & %s & %s & %s & %s & %s & %s \\\\", if (dl == 7) sprintf("%.0f\\%%", 100 * e) else "", dl, ci(r$Re_med, r$Re_lo, r$Re_hi, 2, "--"),
    ci(r$cum12_med, r$cum12_lo, r$cum12_hi, 0, "--"), if (dl == 7) "--" else ci(r$pct_med, r$pct_lo, r$pct_hi, 0, "--"), if (dl == 7) "--" else d0(r$deaths_med, 0), durtxt(r), d0(r$p_end_24m, 2))) }
  if (e != 1) L <- c(L, "\\midrule") }
wtab("scenarios.tex", c(L, "\\bottomrule", "\\end{tabular}"))
# --- fit table
fr <- function(D, key) { x <- D[[key]]; q <- quantile(x, c(.5, .025, .975)); q }
mk <- function(lbl, D, F) { r1 <- fr(D, "r1"); r2 <- fr(D, "r2"); R1 <- fr(D, "Re1"); R2 <- fr(D, "Re2"); ie <- fr(D, "inc_end"); tc <- quantile(D$tc_day, .5)
  c(lbl, ci(r1[1], r1[2], r1[3], 4, "--"), ci(r2[1], r2[2], r2[3], 4, " to "), ci(R1[1], R1[2], R1[3], 2, "--"), ci(R2[1], R2[2], R2[3], 2, "--"),
    sprintf("%s (%s--%s)", dstr(CFG$t0 + tc), dstr(CFG$t0 + min(F$ok)), dstr(CFG$t0 + max(F$ok))), sprintf("%.1f", F$aic0 - F$aic1), ci(ie[1], ie[2], ie[3], 0, "--")) }
fm <- rbind(mk("A: onset curve", fit$DA$draws, fit$A), mk("B: report date", fit$DB$draws, fit$B))
wtab("fit.tex", c("\\begin{tabular}{@{}lccccccc@{}}", "\\toprule",
  "Fit & Growth rate before (/d) & Growth rate after (/d) & $R_e$ before & $R_e$ after & Change date (95\\% interval) & $\\Delta$AIC & Incidence 12 Sep (/d) \\\\", "\\midrule",
  paste0(apply(fm, 1, paste, collapse = " & "), " \\\\"), "\\bottomrule", "\\end{tabular}"))
# --- variants table
Lv <- c("\\begin{tabular}{@{}lccc@{}}", "\\toprule", "Variant & 5-day delay (\\%) & 3-day delay (\\%) & $R_e$ at 3 d \\\\", "\\midrule"); grp <- ""
for (i in seq_len(nrow(vtab))) { r <- vtab[i, ]; if (r$group != grp && r$group != "Reference") Lv <- c(Lv, sprintf("\\multicolumn{4}{@{}l}{\\textit{%s}} \\\\", r$group)); grp <- r$group
  Lv <- c(Lv, sprintf("%s%s & %s & %s & %.2f \\\\", if (r$group == "Reference") "" else "\\quad ", gsub("%", "\\\\%", r$variant), ci(r$pct5, r$pct5_lo, r$pct5_hi, 0, " to "), ci(r$pct3, r$pct3_lo, r$pct3_hi, 0, " to "), r$Re3)) }
wtab("variants.tex", c(Lv, "\\bottomrule", "\\end{tabular}"))
# --- stochastic table
Ls <- c("\\begin{tabular}{@{}rrcccccc@{}}", "\\toprule", "Isolation & Delay & Stochastic cases & Deterministic & P(ended $\\le$12 mo) & P(ended $\\le$24 mo) & P(ended $\\le$60 mo) & Deterministic \\\\",
        "effectiveness & (days) & median (95\\% range) & median & & & & P($\\le$60 mo) \\\\", "\\midrule")
for (e in CFG$eps) { for (dl in CFG$delays) { r <- ss[ss$eps == e & ss$delay == dl, ]
  Ls <- c(Ls, sprintf("%s & %d & %s & %s & %s & %s & %s & %s \\\\", if (dl == 7) sprintf("%.0f\\%%", 100 * e) else "", dl, ci(r$stoch_med, r$stoch_lo, r$stoch_hi, 0, "--"), d0(r$det_med, 0), d0(r$p_end_12m, 2), d0(r$p_end_24m, 2), d0(r$p_end_60m, 2), d0(r$det_p_end_60m, 2))) }
  if (e != 1) Ls <- c(Ls, "\\midrule") }
wtab("stochastic.tex", c(Ls, "\\bottomrule", "\\end{tabular}"))
# --- validation table
Lval <- c("\\begin{tabular}{@{}rcccccc@{}}", "\\toprule", "Horizon & Observed new & Model median & Model 95\\% UI & Predictive 95\\% interval & Model error & Constant-rate error \\\\", "\\midrule")
for (i in 1:nrow(V)) { r <- V[i, ]; Lval <- c(Lval, sprintf("%d d (%s) & %s & %s & %s--%s & %s--%s & %s\\%% & %s\\%% \\\\", r$days_ahead, dstr(r$date), d0(r$observed), d0(r$model_median), d0(r$model_lo), d0(r$model_hi), d0(r$pred_lo), d0(r$pred_hi), d0(r$err_model_pct, 0), d0(r$err_naive_pct, 0))) }
wtab("validation.tex", c(Lval, "\\bottomrule", "\\end{tabular}"))
# --- supplementary tables: PRCC, natural-history sensitivity, by-effectiveness fit, positivity, parameters
wtab("prcc.tex", c("\\begin{tabular}{@{}lc@{}}", "\\toprule", "Parameter (range sampled) & PRCC (95\\% bootstrap CI) \\\\", "\\midrule",
  vapply(seq_len(nrow(prcc)), function(i) sprintf("%s & %s \\\\", lab[prcc$parameter[i]], sprintf("%.2f (%.2f to %.2f)", prcc$prcc[i], prcc$lo[i], prcc$hi[i])), ""), "\\bottomrule", "\\end{tabular}"))
wtab("naturalhistory.tex", c("\\begin{tabular}{@{}cccccc@{}}", "\\toprule", "Latent (d) & Infectious (d) & Neg.\\ log-lik. & $\\beta$ after & Growth rate after (/d) & Incidence 12 Sep (/d) \\\\", "\\midrule",
  vapply(seq_len(nrow(nh)), function(i) sprintf("%.1f & %.0f & %.2f & %.3f & %.4f & %.1f \\\\", nh$latent[i], nh$infectious[i], nh$nll[i], nh$beta2[i], nh$r_after[i], nh$inc_end[i]), ""), "\\bottomrule", "\\end{tabular}"))
wtab("byeffectiveness.tex", c("\\begin{tabular}{@{}clccccc@{}}", "\\toprule", "Fit & Isolation eff. & $\\beta$ before & $\\beta$ after & Growth after (/d) & $R_e$ before & $R_e$ after \\\\", "\\midrule",
  vapply(seq_len(nrow(eps_tab)), function(i) { r <- eps_tab[i, ]; sprintf("%s & %.0f\\%% & %.3f & %.3f & %.4f & %.2f & %.2f \\\\", r$fit, 100 * r$eps, r$beta_before, r$beta_after, r$r_after, r$Re_before, r$Re_after) }, ""), "\\bottomrule", "\\end{tabular}"))
pos <- pr$pos
wtab("positivity.tex", c("\\begin{tabular}{@{}lcccc@{}}", "\\toprule", "Month (2026) & New confirmed cases & Samples analysed & Samples per day & Positivity (\\%) \\\\", "\\midrule",
  vapply(seq_len(nrow(pos)), function(i) sprintf("%s & %s & %s & %.0f & %.1f \\\\", pos$month[i], d0(pos$new_confirmed[i]), d0(pos$tests_analysed[i]), pos$tests_per_day[i], pos$positivity[i]), ""), "\\bottomrule", "\\end{tabular}"))

## =================================================================== in-text numbers (LaTeX macros)
macro_reset(); m <- macro; last <- pr$last; S <- function(e, d) summ[summ$eps == e & summ$delay == d, ]
m("nConfirmed", d0(last$confirmed_case_history)); m("nDeaths", d0(last$confirmed_death_history)); m("crudeCFR", pcs(100 * last$confirmed_death_history / last$confirmed_case_history, 0))
m("nOnsetCases", d0(sum(pr$lat$confirmed_total))); m("shareOnset", pcs(100 * sum(pr$lat$confirmed_total) / last$confirmed_case_history, 0)); m("nRetro", d0(pr$retro)); m("nIntervals", nrow(pr$dB)); m("nOnsetDays", nrow(pr$dA))
m("delayMean", sprintf("%.1f", pr$delay$mean)); m("delayHalf", sprintf("%.0f", pr$delay$q50)); m("delayNinety", sprintf("%.0f", pr$delay$q90)); m("nVintages", pr$delay$n_vintages)
m("posJun", sprintf("%.1f", pos$positivity[1])); m("posSep", sprintf("%.1f", pos$positivity[4])); m("testsJun", sprintf("%.0f", pos$tests_per_day[1])); m("testsSep", sprintf("%.0f", pos$tests_per_day[4]))
m("meanSep", sprintf("%.0f", pr$desc$value[8]))
for (nm in c("A", "B")) { D <- fit[[paste0("D", nm)]]$draws; F <- fit[[nm]]; pre <- paste0("fit", nm)
  m(paste0(pre, "rBefore"), ci(median(D$r1), quantile(D$r1, .025), quantile(D$r1, .975), 3, " to ")); m(paste0(pre, "rAfter"), ci(median(D$r2), quantile(D$r2, .025), quantile(D$r2, .975), 3, " to "))
  m(paste0(pre, "ReBefore"), ci(median(D$Re1), quantile(D$Re1, .025), quantile(D$Re1, .975), 2, " to ")); m(paste0(pre, "ReAfter"), ci(median(D$Re2), quantile(D$Re2, .025), quantile(D$Re2, .975), 2, " to "))
  m(paste0(pre, "ReBeforeMed"), sprintf("%.2f", median(D$Re1))); m(paste0(pre, "ReAfterMed"), sprintf("%.2f", median(D$Re2)))
  m(paste0(pre, "tc"), dstr(CFG$t0 + median(D$tc_day))); m(paste0(pre, "tcRange"), sprintf("%s to %s", dstr(CFG$t0 + min(F$ok)), dstr(CFG$t0 + max(F$ok)))); m(paste0(pre, "dAIC"), sprintf("%.1f", F$aic0 - F$aic1))
  m(paste0(pre, "incEnd"), ci(median(D$inc_end), quantile(D$inc_end, .025), quantile(D$inc_end, .975), 0, " to ")); m(paste0(pre, "incEndMed"), sprintf("%.0f", median(D$inc_end)))
  m(paste0(pre, "pGrow"), sprintf("%.2f", mean(D$r2 > 0))); m(paste0(pre, "cover"), pcs(100 * fit$adq[nm, "coverage95"], 0)); m(paste0(pre, "acf"), sprintf("%.2f", fit$adq[nm, "resid_lag1_acf"])) }
m("nhRMin", d0(min(nh$r_after), 4)); m("nhRMax", d0(max(nh$r_after), 4)); m("nhIncMin", sprintf("%.0f", min(nh$inc_end))); m("nhIncMax", sprintf("%.0f", max(nh$inc_end)))
m("baseCases", ci(S(.8, 7)$cum12_med, S(.8, 7)$cum12_lo, S(.8, 7)$cum12_hi, 0, " to ")); m("baseCasesMed", d0(S(.8, 7)$cum12_med)); m("pReBaseBelow", sprintf("%.2f", S(.8, 7)$p_Re_below1))
for (e in CFG$eps) for (dl in c(5, 3)) { r <- S(e, dl); k <- paste0(em(e), dm(dl))
  m(paste0("pct", k), ci(r$pct_med, r$pct_lo, r$pct_hi, 0, " to ")); m(paste0("pctMed", k), pcs(r$pct_med, 0)); m(paste0("Re", k), ci(r$Re_med, r$Re_lo, r$Re_hi, 2, " to ")); m(paste0("ReMed", k), sprintf("%.2f", r$Re_med))
  m(paste0("avert", k), ci(r$averted_med, r$averted_lo, r$averted_hi, 0, " to ")); m(paste0("averted", k), d0(r$averted_med)); m(paste0("deaths", k), d0(r$deaths_med))
  m(paste0("pEndYearTwo", k), sprintf("%.2f", r$p_end_24m)); m(paste0("pEndYearOne", k), sprintf("%.2f", r$p_end_12m)); m(paste0("pEndYearFive", k), sprintf("%.2f", r$p_end_60m)); m(paste0("pReBelow", k), sprintf("%.2f", r$p_Re_below1))
  m(paste0("dur", k), if (is.finite(r$dur_med)) sprintf("%.0f", months(r$dur_med)) else "$>$60"); m(paste0("durLo", k), if (is.finite(r$dur_lo)) sprintf("%.0f", months(r$dur_lo)) else "$>$60")
  m(paste0("durHi", k), if (is.finite(r$dur_hi)) sprintf("%.0f", months(r$dur_hi)) else "$>$60"); m(paste0("fracIso", k), pcs(100 * r$frac_iso_med, 0)) }
m("fracIsoBase", pcs(100 * S(.8, 7)$frac_iso_med, 0)); m("pEndBaseSixty", sprintf("%.2f", S(.8, 7)$p_end_60m)); m("pEndBaseTwentyFour", sprintf("%.2f", S(.8, 7)$p_end_24m))
sx <- function(e, d, col) ss[ss$eps == e & ss$delay == d, col]
m("stochBase", d0(sx(.8, 7, "stoch_med"))); m("stochFive", d0(sx(.8, 5, "stoch_med"))); m("stochThree", d0(sx(.8, 3, "stoch_med"))); m("stochRedFive", pcs(sx(.8, 5, "pct_red_med"), 0)); m("stochRedThree", pcs(sx(.8, 3, "pct_red_med"), 0))
ratio <- ss$stoch_med / ss$det_med; m("stochRatioMin", sprintf("%.2f", min(ratio))); m("stochRatioMax", sprintf("%.2f", max(ratio)))
m("stochEndThree", sprintf("%.2f", sx(.8, 3, "p_end_12m"))); m("stochEndFive", sprintf("%.2f", sx(.8, 5, "p_end_24m"))); m("stochEndBase", sprintf("%.2f", sx(.8, 7, "p_end_60m")))
for (i in seq_len(nrow(prcc))) m(paste0("prcc", gsub(" ", "", tools::toTitleCase(gsub("_", " ", prcc$parameter[i])))), d0(prcc$prcc[i], 2))
vv <- function(i, col) vtab[i, col]; m("varFiveMin", sprintf("%.0f", min(vtab$pct5[-11]))); m("varFiveMax", sprintf("%.0f", max(vtab$pct5[-11]))); m("varThreeMin", sprintf("%.0f", min(vtab$pct3))); m("varThreeMax", sprintf("%.0f", max(vtab$pct3)))
m("varBaseFourFive", ci(vv(11, "pct5"), vv(11, "pct5_lo"), vv(11, "pct5_hi"), 0, " to ")); m("varBaseFourThree", ci(vv(11, "pct3"), vv(11, "pct3_lo"), vv(11, "pct3_hi"), 0, " to "))
m("varErlangThreeFive", sprintf("%.0f", vv(13, "pct5"))); m("varErlangTenFive", sprintf("%.0f", vv(14, "pct5"))); m("varErlangTenThree", sprintf("%.0f", vv(14, "pct3")))
m("varCondFive", sprintf("%.0f", vv(15, "pct5"))); m("varCondThree", sprintf("%.0f", vv(15, "pct3"))); m("varCondIso", pcs(100 * vv(15, "frac_iso_base"), 0))
m("varRampThirtyFive", sprintf("%.0f", vv(17, "pct5"))); m("varRampThirtyThree", sprintf("%.0f", vv(17, "pct3"))); m("varRampNinetyFive", sprintf("%.0f", vv(18, "pct5"))); m("varRampNinetyThree", sprintf("%.0f", vv(18, "pct3"))); m("varRampNinetyUiLo", sprintf("%.0f", vv(18, "pct5_lo")))
m("varFitAFive", sprintf("%.0f", vv(10, "pct5"))); m("varFitAThree", sprintf("%.0f", vv(10, "pct3"))); m("varFitABase", d0(vv(10, "base_cases")))
m("varIncLoFive", sprintf("%.0f", vv(8, "pct5"))); m("varIncHiFive", sprintf("%.0f", vv(9, "pct5"))); m("varGrowthLoFive", sprintf("%.0f", vv(6, "pct5"))); m("varGrowthHiFive", sprintf("%.0f", vv(7, "pct5"))); m("varGrowthHiThree", sprintf("%.0f", vv(7, "pct3"))); m("varGrowthLoThree", sprintf("%.0f", vv(6, "pct3")))
m("varGrowthLoBase", d0(vv(6, "base_cases"))); m("varGrowthHiBase", d0(vv(7, "base_cases"))); m("varPopFive", sprintf("%.0f", vv(16, "pct5")))
vo <- V[V$days_ahead == max(V$days_ahead), ]; vs <- V[V$days_ahead == min(V$days_ahead), ]
m("valObsLong", d0(vo$observed)); m("valModLong", d0(vo$model_median)); m("valModLongCI", ci(vo$model_median, vo$model_lo, vo$model_hi, 0, " to ")); m("valErrLong", pcs(vo$err_model_pct, 0)); m("valNaiveErrLong", pcs(vo$err_naive_pct, 0))
m("valObsShort", d0(vs$observed)); m("valModShort", d0(vs$model_median)); m("valErrShort", pcs(vs$err_model_pct, 0)); m("valRateObsShort", sprintf("%.0f", vs$observed / vs$days_ahead)); m("valRateObsLong", sprintf("%.0f", vo$observed / vo$days_ahead))
m("valPredLongLo", d0(vo$pred_lo)); m("valPredLongHi", d0(vo$pred_hi)); m("valInsideLong", if (vo$inside_param_interval) "inside" else "outside"); m("valInsideShort", if (vs$inside_param_interval) "inside" else "outside")
m("valCumObs", d0(vo$cum_observed)); m("valCumMod", d0(vo$cum_model))
m("cfrAssumed", pcs(100 * CFG$cfr, 0)); m("popAtRisk", d0(CFG$pop)); m("nDraws", d0(sc$n)); m("nStoch", "500"); m("nPRCC", "1,500")
## critical baseline reproduction number: the largest Re at the 7-day baseline for which the new delay still gives Re < 1
aa <- function(d) (1 / d) / (1 / d + CFG$gamma0); crit <- function(e, d) (1 - e * aa(7)) / (1 - e * aa(d))
Lc <- c("\\begin{tabular}{@{}rcc@{}}", "\\toprule", "Isolation effectiveness & 5-day delay & 3-day delay \\\\", "\\midrule",
        vapply(CFG$eps, function(e) sprintf("%.0f\\%% & %.2f & %.2f \\\\", 100 * e, crit(e, 5), crit(e, 3)), ""), "\\bottomrule", "\\end{tabular}")
wtab("critre.tex", Lc)
for (e in CFG$eps) for (dl in c(5, 3)) m(paste0("critRe", em(e), dm(dl)), sprintf("%.2f", crit(e, dl)))
m("rVersion", sub("R version ", "", R.version.string)); m("desolveVersion", as.character(packageVersion("deSolve"))); m("rYear", R.version$year)
m("nVariants", nrow(vtab) - 1); m("nVarDraws", "300"); m("aIsoFive", pcs(100 * aa(5), 0)); m("aIsoThree", pcs(100 * aa(3), 0)); m("aIsoSeven", pcs(100 * aa(7), 0))
macro_write(file.path(ROOT, "manuscript", "numbers.tex"))
cat(sprintf("Wrote %d macros, 7 figures, 12 tables\n", length(ls(MACROS))))
rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds

######################################################################################################################
## 08 | Optional: build the three report PDFs from the LaTeX sources (only with --pdf)
######################################################################################################################
if (DO_PDF && nzchar(Sys.which("pdflatex"))) local({
banner("building PDFs")
root <- ROOT
## 08_build_reports.R -- builds the three report PDFs from the generated numbers/tables/figures (run R/run_all.R first,
## or at least scripts 01-07). Needs only base R + pdflatex/bibtex (no knitr/Quarto).
##   Rscript R/08_build_reports.R
## Outputs (reports/out/):
##   1_short_report_3pages.pdf     course report: <= 3 pages of body text, references on a separate extra page
##   2_model_fitting_report.pdf    data, reporting-delay correction, likelihood, estimation, diagnostics, fitting code
##   3_full_manuscript.pdf         full journal-style manuscript with R-code appendix
need <- c("manuscript/numbers.tex", "results/tables/scenarios.tex", "results/tables/fit.tex", "figures/fig3_trajectories.pdf")
if (!all(file.exists(file.path(root, need)))) stop("Run Rscript R/run_all.R --no-pdf first (generated tables/numbers are missing).", call. = FALSE)
if (!nzchar(Sys.which("pdflatex"))) stop("pdflatex not found", call. = FALSE)
out <- file.path(root, "reports", "out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
build <- function(dir, tex) {
  old <- setwd(file.path(root, dir)); on.exit(setwd(old)); base <- sub("\\.tex$", "", tex)
  run <- function(cmd, args) system2(cmd, args, stdout = FALSE, stderr = FALSE)
  tx <- function() run("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", base))
  if (tx() != 0) stop("LaTeX failed for ", tex, " (see ", file.path(dir, paste0(base, ".log")), ")", call. = FALSE)
  run("bibtex", base); tx(); if (tx() != 0) stop("LaTeX failed (final pass) for ", tex, call. = FALSE)
  log <- readLines(paste0(base, ".log"), warn = FALSE)
  if (any(grepl("undefined (references|citations)|Citation .* undefined|Reference .* undefined", log))) warning("Undefined references/citations in ", tex)
  file.path(root, dir, paste0(base, ".pdf"))
}
npages <- function(p) { if (nzchar(Sys.which("pdfinfo"))) as.integer(sub("Pages:\\s+", "", grep("^Pages:", system2("pdfinfo", p, stdout = TRUE), value = TRUE))) else NA }
jobs <- list(c("reports", "short_report.tex", "1_short_report_3pages.pdf"),
             c("reports", "fit_report.tex", "2_model_fitting_report.pdf"),
             c("manuscript", "main.tex", "3_full_manuscript.pdf"))
for (j in jobs) {
  p <- build(j[1], j[2]); file.copy(p, file.path(out, j[3]), overwrite = TRUE)
  cat(sprintf("%-32s %s pages\n", j[3], npages(file.path(out, j[3]))))
}
## page-limit check for the short report: body (everything before the references page) must be <= 3 pages
sp <- file.path(out, jobs[[1]][3]); n <- npages(sp)
if (!is.na(n) && nzchar(Sys.which("pdftotext"))) {
  last <- system2("pdftotext", c("-f", n, "-l", n, sp, "-"), stdout = TRUE)
  body <- if (any(grepl("^\\[1\\]", last))) n - 1 else n
  cat(sprintf("Short report body: %d page(s) (limit 3) -> %s\n", body, if (body <= 3) "OK" else "TOO LONG"))
  if (body > 3) stop("Short report exceeds three pages", call. = FALSE)
}
cat("PDFs written to reports/out/\n")
})
writeLines(capture.output(sessionInfo()), file.path(ROOT, "results", "sessionInfo.txt"))
banner("ALL DONE")
