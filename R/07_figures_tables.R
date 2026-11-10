## 07_figures_tables.R -- every figure, table and in-text number of the manuscript is produced here from results/*.
## Outputs: figures/fig*.pdf, results/tables/*.tex, manuscript/numbers.tex (LaTeX macros used in the text)
find_root_ <- function() { a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); normalizePath(file.path(dirname(f[1]), "..")) }
source(file.path(find_root_(), "R", "00_functions.R"))
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
