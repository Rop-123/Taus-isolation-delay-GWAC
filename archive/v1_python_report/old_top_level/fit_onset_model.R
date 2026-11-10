# Fit the SEIQR model to the BVD 2026 symptom-onset curve (INSP SitReps, via epiforecasts/BVDOutbreakSize).
# Run from the repository root:  Rscript R/fit_onset_model.R      (needs only base R + deSolve)
#
# What is fitted
#   Data  : confirmed cases by symptom-onset date, composite of all report vintages (latest vintage per onset date).
#   Model : S-E-I-Q-R, isolation rate tau = 1/delay, residual transmission among isolated (1 - eps).
#           Fixed at baseline: latent 6.3 d, infectious 7 d, delay 7 d, eps 0.8 (sensitivity below).
#   Fitted: transmission beta(t) (constant / one change point / two change points), initial onset incidence,
#           negative-binomial dispersion.
#   Obs.  : onset counts in the latest vintage are right-truncated; the completeness curve F(lag) is estimated
#           from the report vintages and applied to the expected counts.
suppressPackageStartupMessages(library(deSolve))
set.seed(20261031)
stopifnot(file.exists("data/bvd_streams/bvd_2026_onset_curve_all_vintages.csv"))
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

# ------------------------------------------------------------------ settings
T0        <- as.Date("2026-05-18")            # day 0 of the fitted window
LATEST    <- as.Date("2026-09-12")            # report date of the latest vintage (SitRep 121)
EXCL_LAG  <- 10                               # drop onset dates fewer than 10 days before LATEST (F < 0.75)
N_POP     <- 4392200                          # Ituri population (McCabe et al., from scalars file)
W_LOGIT   <- 3                                # width (days) of the logistic change in beta
BASE      <- list(sigma = 1/6.3, gamma = 1/7, tau = 1/7, eps = 0.8)   # baseline natural history / isolation
BASE$gamma_iso <- BASE$gamma
LAG_REPORT <- 9                               # onset-to-result delay (days), used only to compare with report-date counts

# ------------------------------------------------------------------ 1. data
vint <- read.csv("data/bvd_streams/bvd_2026_onset_curve_all_vintages.csv", colClasses = c(sitrep = "character"))
vint$report_date <- as.Date(vint$report_date); vint$onset_date <- as.Date(vint$onset_date)
M <- tapply(vint$confirmed_total, list(as.character(vint$onset_date), as.character(vint$report_date)), sum)
onset_dates <- as.Date(rownames(M)); report_dates <- as.Date(colnames(M))

# 1a. completeness F(lag): ratio of the count at lag k to the (near-final) count >= 35 days later
rat <- NULL
for (i in seq_along(onset_dates)) {
  d <- onset_dates[i]; row <- M[i, ]
  ok_ref <- which(as.numeric(report_dates - d) >= 35 & !is.na(row))
  if (!length(ok_ref)) next
  ref <- ok_ref[which.min(as.numeric(report_dates[ok_ref] - d))]
  base <- row[ref]; if (base < 10) next
  for (j in seq_along(report_dates)) {
    k <- as.numeric(report_dates[j] - d)
    if (k >= 0 && k <= 34 && !is.na(row[j])) rat <- rbind(rat, c(k, row[j] / base, base))
  }
}
rat <- as.data.frame(rat); names(rat) <- c("k", "ratio", "base")
sse <- function(p) sum(rat$base * (rat$ratio - pgamma(rat$k + 1, exp(p[1]), exp(p[2])))^2)
fo <- optim(c(log(8), log(0.6)), sse)
F_shape <- exp(fo$par[1]); F_rate <- exp(fo$par[2])
Fk <- function(k) ifelse(k >= 35, 1, pgamma(k + 1, F_shape, F_rate))
comp_tab <- data.frame(lag = 0:34, empirical_median = tapply(rat$ratio, rat$k, median)[as.character(0:34)],
                       fitted = Fk(0:34))
write.csv(comp_tab, "results/R_fit_completeness.csv", row.names = FALSE)

# 1b. composite onset curve: latest vintage that contains each onset date
curve <- do.call(rbind, lapply(seq_along(onset_dates), function(i) {
  j <- max(which(!is.na(M[i, ]))); data.frame(onset = onset_dates[i], y = M[i, j], lag = as.numeric(report_dates[j] - onset_dates[i]))
}))
curve$F <- Fk(curve$lag)
curve$day <- as.numeric(curve$onset - T0)
fitd <- subset(curve, day >= 0 & lag >= EXCL_LAG)
stopifnot(!anyDuplicated(fitd$day), all(diff(fitd$day) == 1))   # no gaps in the fitted window
TMAX_FIT <- max(fitd$day)
cat(sprintf("Fit window: %s to %s (%d days); completeness F: gamma(shape %.2f, rate %.2f)\n",
            format(T0), format(T0 + TMAX_FIT), nrow(fitd), F_shape, F_rate))

# ------------------------------------------------------------------ 2. model
init_state <- function(b0, inc0, cfg) {
  J <- matrix(c(-cfg$sigma, b0, b0 * (1 - cfg$eps),
                cfg$sigma, -(cfg$tau + cfg$gamma), 0,
                0, cfg$tau, -cfg$gamma_iso), 3, 3, byrow = TRUE)
  ev <- eigen(J); v <- abs(Re(ev$vectors[, which.max(Re(ev$values))]))
  v <- v * inc0 / (cfg$sigma * v[1])
  c(S = N_POP - sum(v), E = v[1], I = v[2], Q = v[3], R = 0, C = 0)
}
beta_fun <- function(t, b, tc) { x <- b[1]; for (k in seq_along(tc)) x <- x + (b[k + 1] - b[k]) * plogis((t - tc[k]) / W_LOGIT); x }
rhs <- function(t, y, p) {
  with(as.list(y), {
    be <- beta_fun(t, p$b, p$tc); lam <- be * I / N_POP + be * (1 - p$cfg$eps) * Q / N_POP
    list(c(-lam * S, lam * S - p$cfg$sigma * E, p$cfg$sigma * E - (p$cfg$tau + p$cfg$gamma) * I,
           p$cfg$tau * I - p$cfg$gamma_iso * Q, p$cfg$gamma * I + p$cfg$gamma_iso * Q, p$cfg$sigma * E))
  })
}
# unconstrained parameter vector -> list(b, tc, inc0, k)
unpack <- function(th, K) {
  b <- exp(th[1:(K + 1)]); tc <- numeric(0); i <- K + 1
  if (K >= 1) { tc[1] <- 5 + 60 * plogis(th[i + 1]); i <- i + 1 }
  if (K >= 2) { tc[2] <- tc[1] + 8 + (112 - tc[1] - 8) * plogis(th[i + 1]); i <- i + 1 }
  list(b = b, tc = tc, inc0 = exp(th[i + 1]), k = exp(th[i + 2]))
}
sim_inc <- function(pp, cfg, tmax) {
  y0 <- init_state(pp$b[1], pp$inc0, cfg)
  out <- ode(y0, 0:(tmax + 1), rhs, list(b = pp$b, tc = pp$tc, cfg = cfg), method = "lsoda", rtol = 1e-7, atol = 1e-6)
  diff(out[, "C"])                                   # incidence on days 0..tmax
}
neglogpost <- function(th, K, cfg) {
  pp <- unpack(th, K)
  if (pp$k > 1e4 || any(pp$b > 5)) return(1e10)
  mu <- sim_inc(pp, cfg, TMAX_FIT)[fitd$day + 1] * fitd$F
  ll <- sum(dnbinom(round(fitd$y), mu = pmax(mu, 1e-9), size = pp$k, log = TRUE))
  lp <- sum(dnorm(th[1:(K + 1)], log(0.25), 1, log = TRUE)) + dnorm(log(pp$inc0), log(30), 1.5, log = TRUE) +
        dnorm(log(pp$k), log(20), 1.5, log = TRUE) + (if (K > 0) sum(dnorm(th[(K + 2):(2 * K + 1)], 0, 2, log = TRUE)) else 0)
  if (!is.finite(ll)) return(1e10)
  -(ll + lp)
}
loglik_only <- function(th, K, cfg) {
  pp <- unpack(th, K); mu <- sim_inc(pp, cfg, TMAX_FIT)[fitd$day + 1] * fitd$F
  sum(dnbinom(round(fitd$y), mu = pmax(mu, 1e-9), size = pp$k, log = TRUE))
}
fit_model <- function(K, cfg, starts) {
  best <- NULL
  for (s in starts) {
    o <- try(optim(s, neglogpost, K = K, cfg = cfg, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10)), silent = TRUE)
    if (inherits(o, "try-error")) next
    o <- optim(o$par, neglogpost, K = K, cfg = cfg, method = "BFGS", control = list(maxit = 500, reltol = 1e-12))
    if (is.null(best) || o$value < best$value) best <- o
  }
  best
}
starts_for <- function(K) {
  lk <- log(20); li <- log(30)
  if (K == 0) return(list(c(log(0.3), li, lk), c(log(0.25), li, lk)))
  if (K == 1) return(lapply(expand.grid(b1 = c(0.2, 0.3), u = c(-1, 0, 1)) |> split(1:6), function(r) c(log(0.35), log(r$b1), r$u, li, lk)))
  lapply(expand.grid(b2 = c(0.2, 0.3), u1 = c(-1, 0), u2 = c(-1, 1)) |> split(1:8), function(r) c(log(0.35), log(0.25), log(r$b2), r$u1, r$u2, li, lk))
}

# ------------------------------------------------------------------ 3. fit and compare (baseline natural history)
fits <- lapply(0:2, function(K) fit_model(K, BASE, starts_for(K)))
cmp <- do.call(rbind, lapply(0:2, function(K) {
  th <- fits[[K + 1]]$par; npar <- length(th); ll <- loglik_only(th, K, BASE)
  data.frame(change_points = K, n_par = npar, loglik = ll, AIC = -2 * ll + 2 * npar, BIC = -2 * ll + log(nrow(fitd)) * npar)
}))
cmp$dAIC <- cmp$AIC - min(cmp$AIC)
print(round(cmp, 1)); write.csv(cmp, "results/R_fit_model_comparison.csv", row.names = FALSE)
Kbest <- cmp$change_points[which.min(cmp$AIC)]
cat("Best model by AIC: change points =", Kbest, "\n")

# ------------------------------------------------------------------ 4. MCMC (adaptive random walk) for the best model
hess <- optimHess(fits[[Kbest + 1]]$par, neglogpost, K = Kbest, cfg = BASE, control = list(ndeps = rep(1e-3, length(fits[[Kbest + 1]]$par))))
Sig <- tryCatch(solve(hess), error = function(e) diag(0.01, nrow(hess)))
if (any(diag(Sig) <= 0) || any(!is.finite(Sig))) Sig <- diag(0.01, nrow(hess))
Sig <- (Sig + t(Sig)) / 2 + diag(1e-8, nrow(Sig))
run_chain <- function(start, n = 15000, K, cfg) {
  d <- length(start); chain <- matrix(NA, n, d); cur <- start; lp <- -neglogpost(cur, K, cfg); acc <- 0
  cov <- Sig; L <- t(chol(cov * 2.38^2 / d))
  for (i in 1:n) {
    if (i > 500 && i %% 500 == 0) { cm <- cov(chain[max(1, i - 4000):(i - 1), , drop = FALSE]); cov <- 0.8 * cm + 0.2 * Sig; L <- t(chol(cov * 2.38^2 / d + diag(1e-10, d))) }
    prop <- cur + as.numeric(L %*% rnorm(d)); lpp <- -neglogpost(prop, K, cfg)
    if (log(runif(1)) < lpp - lp) { cur <- prop; lp <- lpp; acc <- acc + 1 }
    chain[i, ] <- cur
  }
  list(chain = chain, acc = acc / n)
}
th_hat <- fits[[Kbest + 1]]$par
starts_mc <- lapply(1:3, function(i) th_hat + as.numeric(t(chol(Sig)) %*% rnorm(length(th_hat))) * 0.5)
chains <- lapply(starts_mc, run_chain, K = Kbest, cfg = BASE)
burn <- 5000; keep <- seq(burn + 1, 15000, by = 10)
draws <- do.call(rbind, lapply(chains, function(ch) ch$chain[keep, , drop = FALSE]))
rhat <- function(chs) { m <- length(chs); n <- nrow(chs[[1]]); apply(sapply(1:ncol(chs[[1]]), function(j) {
  x <- sapply(chs, function(ch) ch[, j]); W <- mean(apply(x, 2, var)); B <- n * var(colMeans(x)); sqrt(((n - 1) / n * W + B / n) / W) }), 2, identity) }
rh <- rhat(lapply(chains, function(ch) ch$chain[keep, , drop = FALSE]))
ess1 <- function(x) { a <- acf(x, plot = FALSE, lag.max = 200)$acf[-1]; cut <- which(a < 0.05)[1]; if (is.na(cut)) cut <- 200; length(x) / (1 + 2 * sum(a[seq_len(cut - 1)])) }
ess <- apply(draws, 2, ess1)
cat(sprintf("MCMC: acceptance %s; max R-hat %.3f; min ESS %.0f\n", paste(round(sapply(chains, `[[`, "acc"), 2), collapse = "/"), max(rh), min(ess)))

# ------------------------------------------------------------------ 5. posterior summaries
pars <- t(apply(draws, 1, function(th) { pp <- unpack(th, Kbest); unlist(c(setNames(pp$b, paste0("beta", seq_along(pp$b))),
                  setNames(pp$tc, paste0("tc", seq_along(pp$tc))), inc0 = pp$inc0, k = pp$k)) }))
RF <- function(b, cfg = BASE) b * (1 / (cfg$tau + cfg$gamma) + cfg$tau / (cfg$tau + cfg$gamma) * (1 - cfg$eps) / cfg$gamma_iso)
summ <- function(x) c(median = median(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975)))
ptab <- as.data.frame(t(apply(pars, 2, summ))); ptab$parameter <- rownames(ptab)
for (nm in grep("^beta", colnames(pars), value = TRUE)) {
  s <- summ(RF(pars[, nm])); ptab <- rbind(ptab, data.frame(median = s[1], lo = s[2], hi = s[3], parameter = paste0("R_", sub("beta", "", nm), " (=RF*beta)")))
}
if (Kbest >= 1) for (j in 1:Kbest) { d_ <- T0 + pars[, paste0("tc", j)]; ptab <- rbind(ptab, data.frame(median = as.numeric(median(d_)), lo = as.numeric(quantile(d_, .025)), hi = as.numeric(quantile(d_, .975)), parameter = paste0("tc", j, "_date_days_since_epoch"))) }
ptab <- ptab[, c("parameter", "median", "lo", "hi")]; rownames(ptab) <- NULL
write.csv(ptab, "results/R_fit_params.csv", row.names = FALSE); print(ptab, digits = 4)

# ------------------------------------------------------------------ 6. posterior predictive and R_t(t)
TMAX_PRED <- as.numeric(as.Date("2026-10-28") - T0)
sel <- sample(nrow(draws), 400)
inc_pp <- sapply(sel, function(i) sim_inc(unpack(draws[i, ], Kbest), BASE, TMAX_PRED))          # (TMAX_PRED+1) x 400
rt_pp <- sapply(sel, function(i) { pp <- unpack(draws[i, ], Kbest); RF(beta_fun(0:TMAX_PRED, pp$b, pp$tc)) })
qs <- function(m) t(apply(m, 1, quantile, c(.025, .5, .975)))
rt_q <- qs(rt_pp); rt_df <- data.frame(date = T0 + 0:TMAX_PRED, Rt_lo = rt_q[, 1], Rt_med = rt_q[, 2], Rt_hi = rt_q[, 3])
write.csv(rt_df, "results/R_fit_Rt.csv", row.names = FALSE)
cat(sprintf("Model-implied R at end of fit window (%s): %.2f (%.2f-%.2f)\n", format(T0 + TMAX_FIT), rt_df$Rt_med[TMAX_FIT + 1], rt_df$Rt_lo[TMAX_FIT + 1], rt_df$Rt_hi[TMAX_FIT + 1]))
cat(sprintf("Model-implied R at %s: %.2f (%.2f-%.2f)\n", "2026-10-28", rt_df$Rt_med[TMAX_PRED + 1], rt_df$Rt_lo[TMAX_PRED + 1], rt_df$Rt_hi[TMAX_PRED + 1]))

# goodness of fit: Pearson residuals and coverage of the 95% posterior-predictive interval for the counts
mu_pp <- inc_pp[fitd$day + 1, ] * fitd$F
kk <- pars[sel, "k"]
ypp <- sapply(seq_along(sel), function(j) rnbinom(nrow(fitd), mu = mu_pp[, j], size = kk[j]))
pi_lo <- apply(ypp, 1, quantile, .025); pi_hi <- apply(ypp, 1, quantile, .975)
cover <- mean(fitd$y >= pi_lo & fitd$y <= pi_hi)
mu_med <- apply(mu_pp, 1, median); kmed <- median(kk)
pres <- (fitd$y - mu_med) / sqrt(mu_med + mu_med^2 / kmed)
gof <- data.frame(coverage_95PI = cover, pearson_sd = sd(pres), pearson_lag1_autocorr = cor(pres[-1], pres[-length(pres)]),
                  rmse_log1p = sqrt(mean((log1p(fitd$y) - log1p(mu_med))^2)), n = nrow(fitd))
write.csv(gof, "results/R_fit_gof.csv", row.names = FALSE); print(round(gof, 3))

# ------------------------------------------------------------------ 7. out-of-sample check against report-date confirmed totals
w <- read.csv("data/bvd_streams/bvd_2026_streams_wide.csv"); w$date <- as.Date(w$date)
cc <- function(d) w$confirmed_case_history[w$date == as.Date(d)]
onset_complete <- subset(curve, onset >= as.Date("2026-06-01") & onset <= as.Date("2026-08-15"))
rho <- sum(onset_complete$y) / (cc("2026-08-24") - cc("2026-06-10") - 272)       # 272 = retrospective harmonisation on 22 Jul
ext <- read.csv("data/cumulative_confirmed_drc.csv"); ext$date <- as.Date(ext$date)
day_of <- function(d) as.numeric(as.Date(d) - T0)
val <- do.call(rbind, lapply(c("2026-09-27", "2026-10-28"), function(D) {
  from <- day_of("2026-09-12") - LAG_REPORT; to <- day_of(D) - LAG_REPORT - 1
  pred <- cc("2026-09-12") + colSums(inc_pp[(from + 1):(to + 1), , drop = FALSE]) / rho
  data.frame(date = D, observed = ext$cum_confirmed_drc[ext$date == as.Date(D)], pred_median = median(pred), pred_lo = unname(quantile(pred, .025)), pred_hi = unname(quantile(pred, .975)))
}))
val$pct_error <- 100 * (val$pred_median - val$observed) / val$observed
val$inside_95 <- val$observed >= val$pred_lo & val$observed <= val$pred_hi
write.csv(val, "results/R_fit_validation.csv", row.names = FALSE)
cat(sprintf("Onset curve captures about %.0f%% of confirmed cases (rho = %.2f)\n", 100 * rho, rho)); print(round(val[, -1], 1)); print(val$date)

# ------------------------------------------------------------------ 8. sensitivity of the fit to the fixed isolation settings
alt <- list(list(lab = "baseline (delay 7 d, eps 0.8)", tau = 1/7, eps = 0.8), list(lab = "eps 0.5", tau = 1/7, eps = 0.5),
            list(lab = "eps 1.0", tau = 1/7, eps = 1.0), list(lab = "delay 4.6 d", tau = 1/4.6, eps = 0.8), list(lab = "delay 8.6 d", tau = 1/8.6, eps = 0.8),
            list(lab = "latent 10.1 d", tau = 1/7, eps = 0.8, sigma = 1/10.1), list(lab = "infectious 10 d", tau = 1/7, eps = 0.8, gamma = 1/10))
sens <- do.call(rbind, lapply(alt, function(a) {
  cfg <- BASE; for (n in setdiff(names(a), "lab")) cfg[[n]] <- a[[n]]; cfg$gamma_iso <- cfg$gamma
  f <- fit_model(Kbest, cfg, list(th_hat)); pp <- unpack(f$par, Kbest); ll <- loglik_only(f$par, Kbest, cfg)
  bT <- beta_fun(TMAX_FIT, pp$b, pp$tc)
  data.frame(setting = a$lab, loglik = ll, beta_start = pp$b[1], beta_end = bT, R_start = RF(pp$b[1], cfg), R_end = RF(bT, cfg))
}))
write.csv(sens, "results/R_fit_sensitivity.csv", row.names = FALSE); print(sens, digits = 3)

# ------------------------------------------------------------------ 9. figure
pdf("figures/fig_R_fit.pdf", width = 7.4, height = 5.6); par(mfrow = c(2, 2), mar = c(4, 4, 2.2, 0.8), cex = 0.72, mgp = c(2.4, 0.7, 0))
dts <- T0 + 0:TMAX_PRED
mu_all <- inc_pp[1:(TMAX_FIT + 1 + 14), ] * Fk(as.numeric(LATEST - (T0 + 0:(TMAX_FIT + 14))))
x <- T0 + 0:(TMAX_FIT + 14); qm <- qs(mu_all)
plot(curve$onset[curve$day >= 0], curve$y[curve$day >= 0], pch = 16, col = ifelse(curve$lag[curve$day >= 0] >= EXCL_LAG, "black", "grey70"), cex = 0.6,
     xlab = "Symptom-onset date", ylab = "Confirmed cases by onset date", main = "(a) Fit to onset curve (SitRep 121)", xlim = range(x), ylim = c(0, max(curve$y) * 1.05))
polygon(c(x, rev(x)), c(qm[, 1], rev(qm[, 3])), col = rgb(0.9, 0.54, 0, 0.3), border = NA); lines(x, qm[, 2], col = "#e08a00", lwd = 2)
tr <- qs(inc_pp[1:(TMAX_FIT + 15), ]); lines(x, tr[, 2], col = "#0b6e99", lty = 2, lwd = 1.5)
abline(v = LATEST - EXCL_LAG, lty = 3); legend("topleft", bty = "n", cex = 0.8, pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 2), col = c("black", "grey70", "#e08a00", "#0b6e99"),
       legend = c("fitted data", "excluded (truncated)", "expected observed (95% CrI)", "implied true onsets"))
plot(rt_df$date, rt_df$Rt_med, type = "n", ylim = c(0, max(rt_df$Rt_hi, 3)), xlab = "Date", ylab = "Model-implied reproduction number", main = "(b) R(t) from fitted beta(t)")
polygon(c(rt_df$date, rev(rt_df$date)), c(rt_df$Rt_lo, rev(rt_df$Rt_hi)), col = rgb(0.04, 0.43, 0.6, 0.25), border = NA); lines(rt_df$date, rt_df$Rt_med, col = "#0b6e99", lwd = 2); abline(h = 1, lty = 2)
abline(v = LATEST, lty = 3); text(LATEST, 0.15, "end of fit data", pos = 4, cex = 0.8)
plot(fitd$onset, pres, pch = 16, cex = 0.6, xlab = "Symptom-onset date", ylab = "Pearson residual", main = "(c) Residuals", ylim = c(-4, 4)); abline(h = c(-2, 0, 2), lty = c(3, 1, 3))
ob <- ext[ext$date >= as.Date("2026-09-01"), ]; ob <- rbind(data.frame(date = as.Date("2026-09-12"), cum_confirmed_drc = cc("2026-09-12")), ob[, c("date", "cum_confirmed_drc")])
dd <- seq(as.Date("2026-09-13"), as.Date("2026-10-28"), by = "day")
pc <- sapply(dd, function(D) { from <- day_of("2026-09-12") - LAG_REPORT; to <- day_of(D) - LAG_REPORT - 1; cc("2026-09-12") + colSums(inc_pp[(from + 1):(to + 1), , drop = FALSE]) / rho })
pcq <- t(apply(pc, 1, quantile, c(.025, .5, .975)))
plot(dd, pcq[, 2], type = "n", ylim = c(7000, max(pcq[, 3], 10500)), xlab = "Report date", ylab = "Cumulative confirmed DRC cases", main = "(d) Out-of-sample check")
polygon(c(dd, rev(dd)), c(pcq[, 1], rev(pcq[, 3])), col = rgb(0.9, 0.54, 0, 0.3), border = NA); lines(dd, pcq[, 2], col = "#e08a00", lwd = 2)
points(ob$date, ob$cum_confirmed_drc, pch = 16); legend("topleft", bty = "n", cex = 0.8, pch = c(16, NA), lty = c(NA, 1), col = c("black", "#e08a00"), legend = c("observed (CDC, UKHSA, SitRep)", "predicted (95% CrI)"))
dev.off()
cat("Done. Outputs in results/R_fit_*.csv and figures/fig_R_fit.pdf\n")
