# Run the six simulation scenarios in parallel, one background R process
# per scenario. Works on Windows, macOS and Linux.
#
# Usage, from the package root (the folder that contains `analysis/`):
#   source("analysis/run-simulation-parallel.R")
#
# Each process writes its replications to analysis/simulation-checkpoints/
# and its progress log to analysis/simulation-logs/scenario-<k>.txt. No new
# replication starts after STOP_AT. When every process has stopped, run
#   Rscript analysis/simulation-study.R summary
# (or, inside R: system2(rscript, c("analysis/simulation-study.R", "summary")))
# to write simulation-draws.csv and simulation-performance.csv.

STOP_AT <- "2026-10-12 10:00"      # local time; adjust to deadline
SCENARIOS <- 1:6                   # 1-3 MSNBurr, 4-6 MSNBurr-IIa (alpha 0.5, 1, 3)

n_cores <- parallel::detectCores(logical = FALSE)
message("Physical cores: ", n_cores, "; scenarios: ", length(SCENARIOS))
if (!is.na(n_cores) && n_cores < length(SCENARIOS)) {
  message("Fewer cores than scenarios: processes will share cores and each ",
          "will run slower. Consider SCENARIOS <- 1:", n_cores,
          " here and the rest on another computer.")
}

rscript <- file.path(R.home("bin"), "Rscript")
dir.create("analysis/simulation-logs", showWarnings = FALSE, recursive = TRUE)

for (k in SCENARIOS) {
  log_file <- sprintf("analysis/simulation-logs/scenario-%d.txt", k)
  system2(rscript,
          args = c("analysis/simulation-study.R", k, shQuote(STOP_AT)),
          stdout = sub("[.]txt$", "-output.txt", log_file),
          stderr = log_file, wait = FALSE)
  message("Started scenario ", k, " -> ", log_file)
}

message("All processes started. Watch progress with, for example:\n",
        "  readLines('analysis/simulation-logs/scenario-1.txt')\n",
        "or count finished replications with:\n",
        "  table(sub('-r[0-9]+[.]rds$', '', ",
        "list.files('analysis/simulation-checkpoints')))")
