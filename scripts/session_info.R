# Save R session information for reproducibility

output_file <- "results/session_info.txt"

dir.create("results", showWarnings = FALSE)

sink(output_file)

cat("Session information\n")
cat("====================\n\n")

print(sessionInfo())

sink()

cat("Session info saved to:", output_file, "\n")
