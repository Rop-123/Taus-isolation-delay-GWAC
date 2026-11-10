## 06_validation.R -- out-of-sample check. The model is fitted to data up to 12 Sep 2026 ONLY. We forecast confirmed cases reported up to
## 27 Sep and 28 Oct 2026 (baseline: nothing changes) and compare with later official cumulative counts (UKHSA 27 Sep; CDC 28 Oct), which
## were never used in fitting. Benchmark: constant incidence equal to the mean of the last 14 days.
## Outputs: results/validation.rds, results/tables/validation.csv
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
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
