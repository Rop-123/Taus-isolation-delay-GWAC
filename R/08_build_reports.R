## 08_build_reports.R -- builds the three report PDFs from the generated numbers/tables/figures (run R/run_all.R first,
## or at least scripts 01-07). Needs only base R + pdflatex/bibtex (no knitr/Quarto).
##   Rscript R/08_build_reports.R
## Outputs (reports/out/):
##   1_short_report_3pages.pdf     course report: <= 3 pages of body text, references on a separate extra page
##   2_model_fitting_report.pdf    data, reporting-delay correction, likelihood, estimation, diagnostics, fitting code
##   3_full_manuscript.pdf         full journal-style manuscript with R-code appendix
a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); root <- normalizePath(file.path(dirname(f[1]), ".."))
need <- c("manuscript/numbers.tex", "results/tables/scenarios.tex", "results/tables/fit.tex", "figures/fig3_trajectories.pdf")
if (!all(file.exists(file.path(root, need)))) stop("Run Rscript R/run_all.R --no-pdf first (generated tables/numbers are missing).", call. = FALSE)
if (!nzchar(Sys.which("pdflatex"))) stop("pdflatex not found", call. = FALSE)
out <- file.path(root, "reports", "out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
build <- function(dir, tex) {
  old <- setwd(file.path(root, dir)); on.exit(setwd(old)); base <- sub("\\.tex$", "", tex)
  run <- function(cmd, args) system2(cmd, args, stdout = FALSE, stderr = FALSE)
  tx <- function() run("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", base))
  if (tx() != 0) stop("LaTeX failed for ", tex, " (see ", file.path(dir, paste0(base, ".log")), ")", call. = FALSE)
  run("bibtex", base); tx(); if (tx() != 0) stop("LaTeX failed (final pass) for ", tex, call. = FALSE)
  log <- readLines(paste0(base, ".log"), warn = FALSE)
  if (any(grepl("undefined (references|citations)|Citation .* undefined|Reference .* undefined", log))) warning("Undefined references/citations in ", tex)
  file.path(root, dir, paste0(base, ".pdf"))
}
npages <- function(p) { if (nzchar(Sys.which("pdfinfo"))) as.integer(sub("Pages:\\s+", "", grep("^Pages:", system2("pdfinfo", p, stdout = TRUE), value = TRUE))) else NA }
jobs <- list(c("reports", "short_report.tex", "1_short_report_3pages.pdf"),
             c("reports", "fit_report.tex", "2_model_fitting_report.pdf"),
             c("manuscript", "main.tex", "3_full_manuscript.pdf"))
for (j in jobs) {
  p <- build(j[1], j[2]); file.copy(p, file.path(out, j[3]), overwrite = TRUE)
  cat(sprintf("%-32s %s pages\n", j[3], npages(file.path(out, j[3]))))
}
## page-limit check for the short report: body (everything before the references page) must be <= 3 pages
sp <- file.path(out, jobs[[1]][3]); n <- npages(sp)
if (!is.na(n) && nzchar(Sys.which("pdftotext"))) {
  last <- system2("pdftotext", c("-f", n, "-l", n, sp, "-"), stdout = TRUE)
  body <- if (any(grepl("^\\[1\\]", last))) n - 1 else n
  cat(sprintf("Short report body: %d page(s) (limit 3) -> %s\n", body, if (body <= 3) "OK" else "TOO LONG"))
  if (body > 3) stop("Short report exceeds three pages", call. = FALSE)
}
cat("PDFs written to reports/out/\n")
