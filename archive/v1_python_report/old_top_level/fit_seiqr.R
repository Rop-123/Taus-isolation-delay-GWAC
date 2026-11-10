## Fit the SEIQR model to the 2026 DRC BVD surveillance data (base R + deSolve).
## Run from anywhere:  Rscript R/fit_seiqr.R      (about 3-6 minutes)
##
## Data (data/bvd_2026/, epiforecasts/BVDOutbreakSize export, data as of 12 Sep 2026):
##   * onset curve, latest vintage  -> Fit A: confirmed cases by SYMPTOM-ONSET date. The last ~3 weeks are
##     incomplete (right truncation); the completeness curve is measured from the 55 earlier vintages.
##   * cumulative confirmed cases   -> Fit B: confirmed cases by REPORT date, with the 22 Jul harmonisation
##     step (369 reported vs 97 new) removed, 2 Jun - 12 Sep.
##   * tests analysed per day       -> positivity diagnostic (is "flat" incidence a testing artefact?)
## Model: S-E-I-Q-R, isolation rate tau = 1/7 per day (baseline), isolated people transmit with (1-eps).
## Natural history is FIXED at the central literature values; only transmission (beta), the starting
## incidence and the count dispersion are estimated. Optionally one change in beta (changepoint) is allowed.
suppressPackageStartupMessages(library(deSolve))
set.seed(20261101)

args <- commandArgs(trailingOnly = FALSE)
f <- sub("--file=", "", args[grep("--file=", args)])
root <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
dat <- file.path(root, "data", "bvd_2026"); res <- file.path(root, "results"); fig <- file.path(root, "figures")
dir.create(res, FALSE, TRUE); dir.create(fig, FALSE, TRUE)
sink(file.path(res, "fit_R_log.txt"), split = TRUE)
t_start <- Sys.time()

## ------------------------------------------------------------------ fixed quantities
SIGMA <- 1 / 6.3          # latent period 6.3 d
GAMMA <- 1 / 7            # infectious period 7 d (also while isolated)
TAU0  <- 1 / 7            # baseline onset-to-isolation delay 7 d
EPS   <- 0.8              # isolation effectiveness (base case); 0.5 and 1.0 refitted below
POP   <- 4392200          # Ituri population (scalars file); depletion is negligible either way
T0    <- as.Date("2026-05-02")      # model start (before first detection; gives history for delay convolution)
END   <- as.Date("2026-09-12")      # data cut-off
TEND  <- as.numeric(END - T0)

## ------------------------------------------------------------------ data
wide <- read.csv(file.path(dat, "bvd_2026_streams_wide.csv"), stringsAsFactors = FALSE)
wide$date <- as.Date(wide$date)
lat <- read.csv(file.path(dat, "bvd_2026_onset_curve_latest.csv"), stringsAsFactors = FALSE)
lat$onset_date <- as.Date(lat$onset_date); REPORT <- as.Date(lat$report_date[1])
vint <- read.csv(file.path(dat, "bvd_2026_onset_curve_all_vintages.csv"), stringsAsFactors = FALSE,
                 colClasses = c(sitrep = "character"))
vint$onset_date <- as.Date(vint$onset_date); vint$report_date <- as.Date(vint$report_date)
sc <- read.csv(file.path(dat, "bvd_2026_scalars.csv"), stringsAsFactors = FALSE)
cb <- sc[sc$name == "confirmed_break_dates", ]
break_date  <- as.Date(gsub("[^0-9-]", "", cb$value))
break_gross <- as.numeric(sub(".*'gross_cases': \\[(\\d+)\\].*", "\\1", cb$other_fields))
cat("Harmonisation break:", format(break_date), " printed new cases that day:", break_gross, "\n")

## ------------------------------------------------------------------ 1. reporting-delay (completeness) curve
ref <- vint[vint$sitrep == "121" & vint$onset_date <= as.Date("2026-08-15"), c("onset_date", "confirmed_total")]
m <- merge(vint, ref, by = "onset_date", suffixes = c("", "_ref"))
m$lag <- as.integer(m$report_date - m$onset_date); m <- m[m$lag >= 0 & m$lag <= 40, ]
agg <- aggregate(cbind(confirmed_total, confirmed_total_ref) ~ lag, m, sum)
agg$ratio <- agg$confirmed_total / agg$confirmed_total_ref
wfit <- optim(c(1.5, 10), function(p) {
  if (any(p <= 0)) return(1e9)
  sum((agg$ratio - pweibull(agg$lag + 0.5, p[1], p[2]))^2) }, method = "Nelder-Mead")
WK <- wfit$par[1]; WS <- wfit$par[2]
cc <- function(L) pweibull(L + 0.5, WK, WS)                   # share of an onset day's cases visible after L days
pm <- diff(c(0, cc(0:(TEND - 1))))                             # delay pmf, lags 0..TEND-1
PM <- matrix(0, TEND, TEND); for (u in 1:TEND) PM[u, 1:u] <- pm[u:1]
cat(sprintf("Reporting delay (onset -> in the SitRep onset curve): Weibull shape %.2f scale %.1f d, mean %.1f d; 50%% by %.0f d, 90%% by %.0f d\n",
            WK, WS, WS * gamma(1 + 1 / WK), qweibull(.5, WK, WS) - .5, qweibull(.9, WK, WS) - .5))

## Fit A data: onset curve of the latest vintage
dA <- lat[lat$onset_date >= as.Date("2026-05-25"), ]
yA <- dA$confirmed_total; kA <- as.numeric(dA$onset_date - T0); lagA <- as.numeric(REPORT - dA$onset_date)
predA <- function(inc) inc[kA] * cc(lagA)

## Fit B data: confirmed increments between report dates, break-corrected
cum <- wide[!is.na(wide$confirmed_case_history) & wide$date >= as.Date("2026-06-02"), c("date", "confirmed_case_history")]
dB <- data.frame(date = cum$date[-1], prev = cum$date[-nrow(cum)], y = diff(cum$confirmed_case_history))
dB$ndays <- as.numeric(dB$date - dB$prev)
step <- dB$y[dB$date == break_date]; retro <- step - break_gross
dB$y[dB$date == break_date] <- break_gross
cat(sprintf("Fit B: %d report intervals; removed %d retrospective cases (%d reported - %d new) on %s\n",
            nrow(dB), retro, step, break_gross, format(break_date)))
uB <- as.numeric(dB$date - T0); u0B <- as.numeric(dB$prev - T0)
predB <- function(inc) { Cc <- c(0, cumsum(as.vector(PM %*% inc))); Cc[uB + 1] - Cc[u0B + 1] }

## ------------------------------------------------------------------ 2. SEIQR model
jac <- function(beta, eps = EPS, tau = TAU0)
  matrix(c(-SIGMA, beta, beta * (1 - eps), SIGMA, -(tau + GAMMA), 0, 0, tau, -GAMMA), 3, 3, byrow = TRUE)
rhs <- function(t, y, p) {
  lam <- p["beta"] * (y[3] + (1 - p["eps"]) * y[4]) / POP
  list(c(-lam * y[1], lam * y[1] - SIGMA * y[2], SIGMA * y[2] - (p["tau"] + GAMMA) * y[3],
         p["tau"] * y[3] - GAMMA * y[4], GAMMA * (y[3] + y[4]), SIGMA * y[2]))
}
growth_rate <- function(beta, eps = EPS, tau = TAU0) max(Re(eigen(jac(beta, eps, tau))$values))
Re_value <- function(beta, eps = EPS, tau = TAU0) beta * (1 / (tau + GAMMA) + tau / (tau + GAMMA) * (1 - eps) / GAMMA)

## daily symptomatic incidence (sigma*E integrated per day), dates T0+1 ... END; beta changes from b1 to b2 on day tc
sim <- function(i0, b1, b2, tc, eps = EPS, tau = TAU0) {
  ev <- eigen(jac(b1, eps, tau)); v <- abs(Re(ev$vectors[, which.max(Re(ev$values))]))
  v <- v * i0 / (SIGMA * v[1])
  y0 <- c(S = POP - sum(v), E = v[1], I = v[2], Q = v[3], R = 0, C = 0)
  o1 <- ode(y0, 0:min(tc, TEND), rhs, c(beta = b1, tau = tau, eps = eps), method = "lsoda", rtol = 1e-6, atol = 1e-6)
  C <- o1[, "C"]
  if (tc < TEND) {
    y1 <- o1[nrow(o1), -1]
    o2 <- ode(y1, tc:TEND, rhs, c(beta = b2, tau = tau, eps = eps), method = "lsoda", rtol = 1e-6, atol = 1e-6)
    C <- c(C, o2[-1, "C"])
  }
  diff(C)
}

## ------------------------------------------------------------------ 3. likelihood (negative binomial) and fitting
nll_gen <- function(p4, tc, pred, y) {                       # p4 = log i0, log beta1, log beta2, log k
  inc <- tryCatch(sim(exp(p4[1]), exp(p4[2]), exp(p4[3]), tc), error = function(e) NULL)
  if (is.null(inc) || any(!is.finite(inc))) return(1e10)
  mu <- pmax(pred(inc), 1e-8)
  -sum(dnbinom(y, mu = mu, size = exp(p4[4]), log = TRUE))
}
fit_const <- function(pred, y, start) {
  f <- function(p) nll_gen(c(p[1], p[2], p[2], p[3]), TEND, pred, y)
  o <- optim(start, f, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-10))
  o <- optim(o$par, f, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-12))
  list(par4 = c(o$par[1], o$par[2], o$par[2], o$par[3]), nll = o$value, hess = optimHess(o$par, f), tc = TEND, npar = 3)
}
fit_change <- function(pred, y, tc, start4) {
  f <- function(p) nll_gen(p, tc, pred, y)
  o <- optim(start4, f, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10))
  o <- optim(o$par, f, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
  list(par4 = o$par, nll = o$value, hess = optimHess(o$par, f), tc = tc, npar = 4)
}
mvn_draw <- function(mean, hess) {
  cv <- tryCatch(solve(hess), error = function(e) diag(1 / pmax(diag(hess), 1e-6)))
  cv <- (cv + t(cv)) / 2; e <- eigen(cv); vals <- pmax(e$values, 1e-8)
  as.vector(mean + e$vectors %*% (sqrt(vals) * rnorm(length(mean))))
}

analyse <- function(label, pred, y, start_const) {
  cat("\n==== ", label, " ====\n", sep = "")
  c0 <- fit_const(pred, y, start_const)
  cat(sprintf("constant beta : nll %.2f  beta %.4f  r %.4f/d  k %.1f\n", c0$nll, exp(c0$par4[2]),
              growth_rate(exp(c0$par4[2])), exp(c0$par4[4])))
  grid <- as.numeric(seq(as.Date("2026-06-18"), as.Date("2026-08-28"), by = 5) - T0)
  prof <- list()
  for (tc in grid) prof[[as.character(tc)]] <- fit_change(pred, y, tc, c0$par4)
  best <- as.numeric(names(prof))[which.min(sapply(prof, `[[`, "nll"))]
  for (tc in setdiff((best - 4):(best + 4), grid)) prof[[as.character(tc)]] <- fit_change(pred, y, tc, prof[[as.character(best)]]$par4)
  prof <- prof[order(as.numeric(names(prof)))]
  nll_tc <- sapply(prof, `[[`, "nll"); tcs <- as.numeric(names(prof)); best <- tcs[which.min(nll_tc)]
  b <- prof[[as.character(best)]]
  cat(sprintf("changepoint   : nll %.2f  change on %s  beta1 %.4f -> beta2 %.4f  r %.4f -> %.4f /d\n", b$nll,
              format(T0 + best), exp(b$par4[2]), exp(b$par4[3]), growth_rate(exp(b$par4[2])), growth_rate(exp(b$par4[3]))))
  aic0 <- 2 * c0$nll + 2 * 3; aic1 <- 2 * b$nll + 2 * 5          # changepoint date counts as a parameter
  cat(sprintf("AIC constant %.1f  AIC changepoint %.1f  (delta %.1f)\n", aic0, aic1, aic0 - aic1))
  ok <- tcs[nll_tc - min(nll_tc) <= 1.92]
  cat("tc 95% profile interval:", format(T0 + min(ok)), "to", format(T0 + max(ok)), "\n")
  list(const = c0, prof = prof, tcs = tcs, nll_tc = nll_tc, best = b, aic0 = aic0, aic1 = aic1, ok = ok)
}

A <- analyse("Fit A: onset curve (truncation-corrected)", predA, yA, c(log(25), log(0.245), log(20)))
B <- analyse("Fit B: confirmed increments by report date", predB, dB$y, c(log(35), log(0.245), log(30)))

## ------------------------------------------------------------------ 4. uncertainty: draws over (changepoint, parameters)
draws_for <- function(F, pred, nd = 1000, nband = 200) {
  w <- exp(-(F$nll_tc - min(F$nll_tc))); w <- w / sum(w)
  out <- matrix(NA, nd, 9); bands <- list(); cnt <- 0
  colnames(out) <- c("tc_day", "beta1", "beta2", "r1", "r2", "Re1", "Re2", "inc_end", "k")
  for (i in 1:nd) {
    tc <- sample(F$tcs, 1, prob = w); fit <- F$prof[[as.character(tc)]]
    p4 <- mvn_draw(fit$par4, fit$hess)
    p4[2:3] <- pmin(pmax(p4[2:3], log(0.05)), log(1.5))
    inc <- tryCatch(sim(exp(p4[1]), exp(p4[2]), exp(p4[3]), tc), error = function(e) NULL)
    if (is.null(inc)) next
    out[i, ] <- c(tc, exp(p4[2]), exp(p4[3]), growth_rate(exp(p4[2])), growth_rate(exp(p4[3])),
                  Re_value(exp(p4[2])), Re_value(exp(p4[3])), inc[TEND], exp(p4[4]))
    if (cnt < nband) { cnt <- cnt + 1; bands[[cnt]] <- list(inc = inc, mu = pred(inc)) }
  }
  list(draws = as.data.frame(out[complete.cases(out), ]), bands = bands)
}
cat("\nDrawing parameter uncertainty ...\n")
dA_ <- draws_for(A, predA); dB_ <- draws_for(B, predB)

qs <- function(x, d = 4) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), median(x), quantile(x, .025), quantile(x, .975))
summ <- function(F, D, label) {
  d <- D$draws
  data.frame(dataset = label, n_obs = NA, nll_const = F$const$nll, nll_change = F$best$nll,
             AIC_const = F$aic0, AIC_change = F$aic1,
             change_date = format(T0 + median(d$tc_day)),
             change_date_95 = paste(format(T0 + min(F$ok)), "to", format(T0 + max(F$ok))),
             r_before = qs(d$r1), r_after = qs(d$r2), Re_before = qs(d$Re1, 2), Re_after = qs(d$Re2, 2),
             beta_before = qs(d$beta1), beta_after = qs(d$beta2), incidence_12Sep = qs(d$inc_end, 1),
             p_growth_after_change = mean(d$r2 > 0), stringsAsFactors = FALSE)
}
tabA <- summ(A, dA_, "A: onset curve (truncation-corrected)"); tabA$n_obs <- length(yA)
tabB <- summ(B, dB_, "B: confirmed by report date (break-corrected)"); tabB$n_obs <- nrow(dB)
tab <- rbind(tabA, tabB)
write.csv(tab, file.path(res, "fit_R_parameters.csv"), row.names = FALSE)
write.csv(cbind(fit = "A", dA_$draws), file.path(res, "fit_R_draws_A.csv"), row.names = FALSE)
write.csv(cbind(fit = "B", dB_$draws), file.path(res, "fit_R_draws_B.csv"), row.names = FALSE)
print(t(tab), quote = FALSE)

## ------------------------------------------------------------------ 5. beta for each isolation effectiveness (refit at best changepoint)
cat("\nBeta for each isolation effectiveness (refit):\n")
eps_rows <- list()
for (e in c(0.5, 0.8, 1.0)) {
  EPS <- e
  for (nm in c("A", "B")) {
    F <- if (nm == "A") A else B; pred <- if (nm == "A") predA else predB; y <- if (nm == "A") yA else dB$y
    fit <- fit_change(pred, y, F$best$tc, F$best$par4)
    eps_rows[[length(eps_rows) + 1]] <- data.frame(fit = nm, eps = e, nll = fit$nll,
      beta_before = exp(fit$par4[2]), beta_after = exp(fit$par4[3]),
      r_before = growth_rate(exp(fit$par4[2]), e), r_after = growth_rate(exp(fit$par4[3]), e),
      Re_before = Re_value(exp(fit$par4[2]), e), Re_after = Re_value(exp(fit$par4[3]), e))
  }
}
EPS <- 0.8
epstab <- do.call(rbind, eps_rows); write.csv(epstab, file.path(res, "fit_R_by_effectiveness.csv"), row.names = FALSE)
print(cbind(epstab[1:2], round(epstab[-(1:2)], 4)), row.names = FALSE)

## ------------------------------------------------------------------ 6. testing diagnostic: positivity by month
tst <- wide[!is.na(wide$tests_analysed_daily_history), c("date", "tests_analysed_daily_history")]
dB$tests <- sapply(seq_len(nrow(dB)), function(j) {
  days <- tst$tests_analysed_daily_history[tst$date > dB$prev[j] & tst$date <= dB$date[j]]
  if (length(days)) mean(days) * dB$ndays[j] else NA })
dB$month <- format(dB$date, "%b")
pos <- aggregate(cbind(y, tests) ~ month, dB[!is.na(dB$tests), ], sum)
pos$positivity_pct <- 100 * pos$y / pos$tests
pos$tests_per_day <- pos$tests / aggregate(ndays ~ month, dB[!is.na(dB$tests), ], sum)$ndays
pos <- pos[match(c("Jun", "Jul", "Aug", "Sep"), pos$month), ]
names(pos)[2:3] <- c("new_confirmed", "tests_analysed")
write.csv(pos, file.path(res, "fit_R_positivity.csv"), row.names = FALSE)
cat("\nPositivity (new confirmed / tests analysed):\n"); print(round(pos[, -1], 1), row.names = FALSE)

## ------------------------------------------------------------------ 7. figure
q3 <- function(lst, key, idx = NULL) {
  M <- do.call(rbind, lapply(lst, function(z) if (is.null(idx)) z[[key]] else z[[key]][idx]))
  apply(M, 2, quantile, c(.025, .5, .975))
}
col1 <- "#0b6e99"; col2 <- "#e08a00"
pdf(file.path(fig, "fig_fit_R.pdf"), width = 7.4, height = 6.3)
par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2, 0.8), cex = 0.68, mgp = c(2.3, 0.7, 0))
## (a) completeness
plot(agg$lag, agg$ratio, pch = 16, ylim = c(0, 1.1), xlab = "Days between onset and report", ylab = "Share of final count reported",
     main = "(a) Reporting delay measured from 55 vintages")
lines(0:40, cc(0:40), col = col1, lwd = 2); abline(h = 1, lty = 3)
legend("bottomright", bty = "n", legend = sprintf("Weibull fit, mean %.1f d", WS * gamma(1 + 1 / WK)), lty = 1, col = col1, lwd = 2)
## (b) fit A
bandA <- q3(dA_$bands, "mu"); nowA <- q3(dA_$bands, "inc", kA)
xA <- dA$onset_date
plot(xA, yA, type = "n", ylim = c(0, max(yA, bandA) * 1.05), xlab = "Symptom-onset date", ylab = "Confirmed cases per day",
     main = "(b) Fit A: onset curve, truncation-corrected")
polygon(c(xA, rev(xA)), c(bandA[1, ], rev(bandA[3, ])), col = adjustcolor(col2, .25), border = NA)
lines(xA, nowA[2, ], col = col1, lwd = 1.6, lty = 2)
lines(xA, bandA[2, ], col = col2, lwd = 2)
points(xA, yA, pch = ifelse(lagA >= 21, 16, 1), cex = 0.6)
abline(v = T0 + median(dA_$draws$tc_day), lty = 3, col = "grey40")
legend("topleft", bty = "n", cex = 0.9, pch = c(16, 1, NA, NA), lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 1.6),
       col = c("black", "black", col2, col1), legend = c("observed (complete)", "observed (last 3 weeks, incomplete)", "model, as it would be reported", "model, true onset incidence"))
## (c) fit B
bandB <- q3(dB_$bands, "mu") / rep(dB$ndays, each = 3)
plot(dB$date, dB$y / dB$ndays, pch = 16, cex = 0.6, ylim = c(0, max(dB$y / dB$ndays, bandB) * 1.05),
     xlab = "Report date", ylab = "Confirmed cases per day", main = "(c) Fit B: confirmed cases by report date")
polygon(c(dB$date, rev(dB$date)), c(bandB[1, ], rev(bandB[3, ])), col = adjustcolor(col2, .25), border = NA)
lines(dB$date, bandB[2, ], col = col2, lwd = 2)
points(dB$date, dB$y / dB$ndays, pch = 16, cex = 0.6)
abline(v = T0 + median(dB_$draws$tc_day), lty = 3, col = "grey40")
text(break_date, max(dB$y / dB$ndays) * 1.02, "22 Jul step\nremoved", cex = 0.8, pos = 4)
## (d) profile over changepoint
rng <- range(c(A$nll_tc - min(A$nll_tc), B$nll_tc - min(B$nll_tc)))
plot(T0 + A$tcs, -(A$nll_tc - min(A$nll_tc)), type = "b", pch = 16, cex = 0.6, col = col1, lwd = 1.5, ylim = c(-8, 0.3),
     xlab = "Date transmission changes", ylab = "Profile log-likelihood (relative to best)", main = "(d) When did transmission change?")
lines(T0 + B$tcs, -(B$nll_tc - min(B$nll_tc)), type = "b", pch = 17, cex = 0.6, col = col2, lwd = 1.5)
abline(h = -1.92, lty = 2); text(T0 + A$tcs[4], -1.92, "95% limit", pos = 1, cex = 0.8)
legend("topright", bty = "n", pch = c(16, 17), col = c(col1, col2), lty = 1, legend = c("Fit A (onset)", "Fit B (report)"))
dev.off()

cat(sprintf("\nFinished in %.1f minutes\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
sink()
