## 04_stochastic.R -- stochastic (binomial-chain) version of the model: do the conclusions survive chance events?
## Same posterior draws as 03_scenarios.R (first n = 500). Transmission is calibrated in discrete time so the stochastic model has the
## same expected growth rate as the deterministic one. Outputs: results/stochastic.rds, results/tables/stochastic_summary.csv
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
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
