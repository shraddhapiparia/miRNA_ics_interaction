# run_all.R
# Render the analysis notebooks in order

rmd_files <- c(
  "analysis/01_mirna_ics_interaction_analysis.Rmd",
  "analysis/02_mirna_cis_miR_QTL.Rmd",
  "analysis/03_miR_Targets.Rmd",
  "analysis/04_mediation_rs36047_miR584.Rmd"
)

required_pkgs <- c("rmarkdown")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

for (f in rmd_files) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}

message("Starting analysis pipeline...")
for (f in rmd_files) {
  message("Rendering: ", f)
  rmarkdown::render(
    input = f,
    output_format = "html_document",
    clean = TRUE,
    envir = new.env(parent = globalenv())
  )
}

message("Pipeline complete.")