## 00_functions.R -- shared model code, helpers and constants (sourced by every script)
## Model: S-E-I-Q-R with an isolated compartment Q. The infectious (not yet isolated) class is split into
## m Erlang stages (m = 1 gives the plain exponential delay used in the main analysis; m > 1 is a
## structural sensitivity analysis). Isolation occurs on leaving stage m at rate m*tau (tau = 1/mean delay).
suppressPackageStartupMessages(library(deSolve))

## ------------------------------------------------------------------ paths
find_root <- function() {
  a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(".")
}
ROOT <- find_root()
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
