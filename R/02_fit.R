## 02_fit.R -- fit the SEIQR model to the 2026 DRC surveillance data (maximum likelihood, negative-binomial counts).
## Natural-history parameters are fixed at central literature values; estimated: starting incidence, transmission rate
## (optionally changing once on an estimated date) and count dispersion. Uncertainty: profile likelihood over the change
## date + multivariate-normal (Hessian) approximation for the other parameters.
## Outputs: results/fit.rds, results/fit_draws_{A,B}.csv, results/tables/fit_*.csv
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
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
