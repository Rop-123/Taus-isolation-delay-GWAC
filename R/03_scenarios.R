## 03_scenarios.R -- nine intervention scenarios (delay 7/5/3 d  x  isolation effectiveness 50/80/100%) from 12 Sep 2026.
## Uncertainty: each of the n draws carries (growth rate r, starting incidence) from the fitted posterior (Fit B, joint with the change date)
## plus latent and infectious periods drawn from triangular distributions. Beta is re-calibrated per draw and per effectiveness so that the
## 7-day baseline reproduces that draw's fitted growth rate; scenarios then change ONLY the isolation rate. All scenarios share the same draws
## (paired), so cases averted are differences within a draw.
## Also computes a continuous policy surface (delay 2--10 d x effectiveness 30--100%).
## Outputs: results/scenarios.rds, results/tables/scenario_summary.csv, results/policy_surface.csv
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
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
