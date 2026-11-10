## run_all.R -- reproduces every number, table and figure, then (optionally) builds the manuscript PDF.
## Usage (from the repository root):   Rscript R/run_all.R            # analysis + the three PDFs (reports/out/)
##                                      Rscript R/run_all.R --no-pdf   # analysis only
## Runtime: about 12-15 minutes on one core (02_fit ~2 min, 03_scenarios ~5 min, 05_sensitivity ~2 min, rest < 1 min each).
a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); root <- normalizePath(file.path(dirname(f[1]), ".."))
setwd(root); rs <- file.path(R.home("bin"), "Rscript")
run <- function(s) { cat(sprintf("\n=== %s (%s)\n", s, format(Sys.time(), "%H:%M:%S"))); st <- system2(rs, s); if (st != 0) stop("failed: ", s, call. = FALSE) }
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
run("tests/test_model.R")
for (s in c("01_data_prep", "02_fit", "03_scenarios", "04_stochastic", "05_sensitivity", "06_validation", "07_figures_tables")) run(file.path("R", paste0(s, ".R")))
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
if (!"--no-pdf" %in% commandArgs(TRUE) && nzchar(Sys.which("pdflatex"))) run("R/08_build_reports.R")   # three PDFs -> reports/out/
cat("\nAll done.\n")
