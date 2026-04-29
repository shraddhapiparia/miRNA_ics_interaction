# run_all.R
# Run the analysis pipeline in order using source()

out_dir <- "results"

# Input data directory checks
required_dirs <- c("data/CAMP", "data/CRA")
for (d in required_dirs) {
  if (!dir.exists(d)) stop("Required input directory missing: ", d)
}

# Create required output directories
output_dirs <- c(
  out_dir,
  file.path(out_dir, "qc"),
  file.path(out_dir, "figures"),
  file.path(out_dir, "models"),
  file.path(out_dir, "qtl"),
  file.path(out_dir, "mediation"),
  file.path(out_dir, "enrichment")
)
for (d in output_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Script list (run in order)
scripts <- c(
  "analysis/01_mirna_ics_interaction_analysis.R",
  "analysis/02_mirna_cis_miR_QTL.R",
  "analysis/03_miR_Targets.R",
  "analysis/04_mediation_rs36047_miR584.R"
)

# Existence checks
for (f in scripts) {
  if (!file.exists(f)) stop("Script not found: ", f)
}

message("Starting analysis pipeline...")

for (f in scripts) {
  message("[START] ", f, " — ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  tryCatch(
    source(f),
    error = function(e) stop("Error in ", f, ": ", conditionMessage(e))
  )
  message("[END]   ", f, " — ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

message("Pipeline complete.")
