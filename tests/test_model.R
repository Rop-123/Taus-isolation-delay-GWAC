## Unit tests for the model code. Run: Rscript tests/test_model.R   (stops with an error if any check fails)
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "..", "R", "00_functions.R"))
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
