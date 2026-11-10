## make_single_script.R -- merges tests + scripts 00-08 into ONE self-contained script (taus_all_in_one.R) at the repository root.
## The modular files in R/ remain the source of truth; re-run this after editing them:   Rscript R/make_single_script.R
a <- commandArgs(FALSE); f <- sub("--file=", "", a[grep("--file=", a)]); root <- normalizePath(file.path(dirname(f[1]), ".."))
rd <- function(p) readLines(file.path(root, p), warn = FALSE)
drop <- function(x, pats) x[!Reduce(`|`, lapply(pats, function(p) grepl(p, x, fixed = TRUE)))]
sec <- function(title, code) c("", strrep("#", 118), paste0("## ", title), strrep("#", 118), code, "rm(list = setdiff(ls(all.names = TRUE), KEEP))   # clean workspace; sections share results only through results/*.rds")
# --- 00: shared functions, defined in the global environment; root comes from the single ROOT defined at the top
f00 <- rd("R/00_functions.R"); i1 <- grep("^find_root <- function", f00); i2 <- grep("^ROOT <- find_root", f00)
f00 <- f00[-(i1:i2)]; f00 <- f00[!grepl("^suppressPackageStartupMessages\\(library\\(deSolve\\)\\)", f00)]
mods <- list(c("01_data_prep.R", "01 | Data preparation: reporting-delay curve and the two fitting datasets"),
             c("02_fit.R", "02 | Maximum-likelihood fits, change-point profile, parameter draws, adequacy checks"),
             c("03_scenarios.R", "03 | Nine scenarios x paired draws; policy surface"),
             c("04_stochastic.R", "04 | Stochastic check"),
             c("05_sensitivity.R", "05 | Global (PRCC) and structural sensitivity analyses"),
             c("06_validation.R", "06 | Out-of-sample validation"),
             c("07_figures_tables.R", "07 | Figures, LaTeX tables and in-text number macros"))
clean <- function(x) drop(x, c("find_root_ <- function", "source(file.path(find_root_()"))
tst <- drop(rd("tests/test_model.R"), "source(file.path(dirname(normalizePath(sub(")
bld <- rd("R/08_build_reports.R"); bld <- bld[!grepl("^a <- commandArgs\\(FALSE\\)", bld)]; bld <- c("root <- ROOT", bld)
hdr <- c(
"## taus_all_in_one.R -- the whole analysis in ONE script (generated from R/*.R by R/make_single_script.R; do not edit by hand).",
"## How fast is fast enough? Onset-to-isolation delay, 2026 Bundibugyo virus disease outbreak, DRC (Team TAUS).",
"##",
"## Run from anywhere:   Rscript taus_all_in_one.R              # tests, all analysis, figures, tables and the numbers quoted in the reports",
"##                      Rscript taus_all_in_one.R --pdf        # additionally build the three report PDFs (needs pdflatex)",
"##",
"## CONTENTS: 00 shared model code | TESTS | 01 data prep | 02 model fitting | 03 scenarios | 04 stochastic | 05 sensitivity | 06 validation | 07 figures, tables, numbers",
"## Needs: R >= 4.3 + package deSolve; for the PDFs pdflatex + bibtex. Runtime about 12-15 min on one core. All seeds are fixed.",
"## Reads data/ ; writes results/, figures/, manuscript/numbers.tex and reports/out/*.pdf. Folder layout is explained in README.md.",
"## The workspace is cleared between sections (as when the separate scripts ran in separate processes); sections pass results via results/*.rds.",
"suppressPackageStartupMessages(library(deSolve)); T0 <- Sys.time(); DO_PDF <- \"--pdf\" %in% commandArgs(TRUE)",
"ROOT <- local({ f <- sub(\"--file=\", \"\", grep(\"--file=\", commandArgs(FALSE), value = TRUE)); r <- if (length(f)) dirname(normalizePath(f[1])) else normalizePath(\".\")",
"  if (!dir.exists(file.path(r, \"data\", \"bvd_2026\"))) stop(\"Cannot find data/bvd_2026 under \", r, \" -- run from the repository root.\", call. = FALSE); r })",
"setwd(ROOT); banner <- function(s) cat(sprintf(\"\\n=== %s  [%.1f min elapsed]\\n\", s, as.numeric(difftime(Sys.time(), T0, units = \"mins\"))))")
out <- hdr
out <- c(out, "", strrep("#", 118), "## 00 | Shared model code, constants and helpers (global environment)", strrep("#", 118), f00)
out <- c(out, "", "KEEP <- c(ls(all.names = TRUE), \"KEEP\")   # objects that survive between sections (shared model code, constants, helpers)")
out <- c(out, sec("TESTS | Unit tests of the model code (stop on failure)", c("banner(\"unit tests\")", tst)))
for (m in mods) out <- c(out, sec(m[2], c(sprintf("banner(\"%s\")", m[2]), clean(rd(file.path("R", m[1]))))))
out <- c(out, "", strrep("#", 118), "## 08 | Optional: build the three report PDFs from the LaTeX sources (only with --pdf)", strrep("#", 118),
         "if (DO_PDF && nzchar(Sys.which(\"pdflatex\"))) local({", "banner(\"building PDFs\")", bld, "})",
         "writeLines(capture.output(sessionInfo()), file.path(ROOT, \"results\", \"sessionInfo.txt\"))",
         "banner(\"ALL DONE\")")
writeLines(out, file.path(root, "taus_all_in_one.R")); cat("Wrote taus_all_in_one.R:", length(out), "lines\n")
