## 05_sensitivity.R -- (a) global sensitivity (partial rank correlation) and (b) one-at-a-time / structural variants.
## Outputs: results/tables/prcc.csv, results/tables/variants.csv
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
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
