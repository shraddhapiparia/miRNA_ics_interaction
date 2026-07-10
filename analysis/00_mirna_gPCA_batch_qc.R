set.seed(42)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
  library(preprocessCore)
  library(ggplot2)
})

# --- Package checks ---
required_pkgs <- c("gPCA", "preprocessCore", "dplyr", "readr", "tibble")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: BiocManager::install(c('gPCA', 'preprocessCore'))")
}
library(gPCA)

# -----------------------------
# Inputs (mirror 01)
# -----------------------------
camp_counts_path <- "data/CAMP/CAMP.492.COUNT.csv"
camp_pheno_path  <- "data/CAMP/camp_pheno_2022-08-04.tsv"
camp_link_path   <- "data/CAMP/shortened_phenotype_files_ver5.csv"

cra_counts_path     <- "data/CRA/2022-08-16_CRA_1159_exceRpt_miRNA_ReadCountscopy.csv"
cra_link_path       <- "data/CRA/cra_shortened_phenotype_files10-11.csv"
cra_hosp_pheno_path <- "data/CRA/COS_TRIO_pheno_1165.csv"
cra_er_pheno_path   <- "data/CRA/PhenotypeInfo.csv"

out_dir <- if (exists("out_dir")) out_dir else "results"
qc_dir  <- file.path(out_dir, "qc", "gPCA")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

min_count        <- 5
min_prop_samples <- 0.5
nperm            <- 1000

batch_candidates <- c("seq_run", "sequencing_run", "batch", "flowcell", "lane")

# -----------------------------
# Helpers (identical logic to 01)
# -----------------------------
filter_mirnas <- function(count_df, min_count = 5, min_prop_samples = 0.5) {
  mat <- as.matrix(count_df[, -1, drop = FALSE])
  rownames(mat) <- count_df[[1]]
  storage.mode(mat) <- "numeric"
  mat[rowSums(mat > min_count) >= ncol(mat) * min_prop_samples, , drop = FALSE]
}

norm_quantile_log <- function(mat) {
  qn <- preprocessCore::normalize.quantiles(as.matrix(mat))
  rownames(qn) <- rownames(mat); colnames(qn) <- colnames(mat)
  log2(qn + 1)
}

find_existing_file <- function(paths) paths[file.exists(paths)][1]

find_batch_col <- function(pheno, candidates, cohort_name) {
  hit_idx <- which(tolower(names(pheno)) %in% tolower(candidates))
  if (length(hit_idx) == 0) {
    stop(
      cohort_name, ": no batch column found among {",
      paste(candidates, collapse = ", "), "}.\n",
      "Available phenotype columns:\n  ",
      paste(names(pheno), collapse = "\n  ")
    )
  }
  names(pheno)[hit_idx[1]]
}

# -----------------------------
# Load data
# -----------------------------
load_camp_qc <- function() {
  counts <- read.csv(camp_counts_path, check.names = FALSE)
  pheno  <- read.table(camp_pheno_path, sep = "\t", header = TRUE, check.names = FALSE)
  link   <- read.csv(camp_link_path, check.names = FALSE)

  mat   <- filter_mirnas(counts, min_count, min_prop_samples)
  pheno <- pheno[pheno$TG != 2, , drop = FALSE] |>
    merge(link, by.x = "camp", by.y = "CAMP_ID", all = FALSE)
  rownames(pheno) <- pheno$`row.names(CAMP_492_COUNT.Trans)`

  common  <- intersect(colnames(mat), rownames(pheno))
  log_mat <- norm_quantile_log(mat[, common, drop = FALSE])
  list(log = log_mat, pheno = pheno[common, , drop = FALSE])
}

load_cra_qc <- function() {
  counts <- read.csv(cra_counts_path, check.names = FALSE)
  link   <- read.csv(cra_link_path, check.names = FALSE)

  # Use whichever pheno file exists; only need it for batch column lookup
  pheno_path <- find_existing_file(c(cra_hosp_pheno_path, cra_er_pheno_path))
  if (is.na(pheno_path)) stop("No CRA phenotype file found")
  pheno <- read.csv(pheno_path, check.names = FALSE)
  pheno <- merge(pheno, link, by.x = "S_SUBJECTID", by.y = "ST.ID")

  mat <- filter_mirnas(counts, min_count, min_prop_samples)
  pheno$sample.id <- gsub("-", ".", pheno$sample.id)
  rownames(pheno)  <- pheno$sample.id

  common  <- intersect(colnames(mat), rownames(pheno))
  log_mat <- norm_quantile_log(mat[, common, drop = FALSE])
  list(log = log_mat, pheno = pheno[common, , drop = FALSE])
}

# -----------------------------
# gPCA runner + plots
# -----------------------------
run_gpca <- function(log_mat, pheno, cohort_name, prefix) {
  batch_col <- find_batch_col(pheno, batch_candidates, cohort_name)
  batch     <- pheno[[batch_col]]
  X         <- t(log_mat)   # samples x miRNAs

  cat("\n--- gPCA:", cohort_name, "---\n")
  cat("  Batch column  :", batch_col, "\n")
  cat("  n_samples     :", nrow(X), "\n")
  cat("  n_features    :", ncol(X), "\n")
  cat("  Batch levels  :", paste(sort(unique(as.character(batch))), collapse = ", "), "\n")


  set.seed(42)
  out <- gPCA::gPCA.batchdetect(
    x      = X,
    batch  = as.factor(batch),
    filt   = NULL,
    nperm  = nperm,
    center = FALSE,
    scaleY = FALSE,
    seed   = 42
  )

  cum_g       <- out$cumulative.var.g
  n_guided_pc <- length(out$cumulative.var.g)

  # --- Summary table ---
  perm_df <- data.frame(delta_star = out$delta.p)

  summary_row <- tibble(
    cohort              = cohort_name,
    batch_column        = batch_col,
    n_samples           = nrow(X),
    n_features          = ncol(X),
    n_batch_levels      = length(unique(batch)),
    delta               = out$delta,
    p_value             = out$p.val,
    varPCu1             = out$varPCu1,
    varPCg1             = out$varPCg1,
    cumulative_var_gPC1 = out$cumulative.var.g[1],
    cumulative_var_gPC2 = if (length(out$cumulative.var.g) >= 2) out$cumulative.var.g[2] else NA_real_,
    cumulative_var_gPC3 = if (length(out$cumulative.var.g) >= 3) out$cumulative.var.g[3] else NA_real_
  )
  write_csv(summary_row, file.path(qc_dir, paste0(prefix, "_gPCA_summary.csv")))
  cat("  delta  :", round(out$delta, 4), "\n")
  cat("  p-value:", round(out$p.val,  4), "\n")

  batch_fac <- factor(batch)

  # --- Plot 1: permutation null distribution ---
  pct_lbl <- function(x, d) paste0(round(x * 100, d), "%")

  p_perm <- ggplot(perm_df, aes(x = delta_star)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_vline(xintercept = out$delta, color = "firebrick", linewidth = 1) +
    annotate("text",
      x = out$delta, y = Inf,
      label = paste0(" observed\n delta = ", round(out$delta, 3),
                     "\n p = ", round(out$p.val, 3)),
      hjust = -0.05, vjust = 1.3, color = "firebrick", size = 3.5
    ) +
    labs(
      title = paste0(cohort_name, ": gPCA permutation null distribution"),
      x = "Permutation delta", y = "Count"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_permutation_delta.png")), p_perm, width = 7, height = 4, dpi = 300)
  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_permutation_delta.pdf")), p_perm, width = 7, height = 4)

  # --- Plot 2: unguided PC1 vs PC2 ---
  n_u <- ncol(out$PCu)
  pc_u_df <- as.data.frame(out$PCu[, seq_len(min(2L, n_u)), drop = FALSE])
  if (ncol(pc_u_df) == 1L) pc_u_df$PC2 <- 0
  names(pc_u_df) <- c("PC1", "PC2")
  pc_u_df$batch <- batch_fac

  p_pcu <- ggplot(pc_u_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2.2, alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0(cohort_name, ": unguided PCA by batch"),
    x     = "PC1",
    y     = "PC2",
    color = batch_col
  ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_unguided_PC1_PC2.png")), p_pcu, width = 7, height = 5, dpi = 300)
  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_unguided_PC1_PC2.pdf")), p_pcu, width = 7, height = 5)

  # --- Plot 3: guided PC1 vs PC2 ---
  n_g <- ncol(out$PCg)
  pc_g_df <- as.data.frame(out$PCg[, seq_len(min(2L, n_g)), drop = FALSE])
  if (ncol(pc_g_df) == 1L) pc_g_df$gPC2 <- 0
  names(pc_g_df) <- c("gPC1", "gPC2")
  pc_g_df$batch <- batch_fac

  p_pcg <- ggplot(pc_g_df, aes(x = gPC1, y = gPC2, color = batch)) +
    geom_point(size = 2.2, alpha = 0.85) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(cohort_name, ": guided PCA by batch"),
      x     = "gPC1",
      y     = if (n_g >= 2) "gPC2" else "gPC2 (n/a — only 1 guided PC)",
      color = batch_col
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_guided_PC1_PC2.png")), p_pcg, width = 7, height = 5, dpi = 300)
  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_guided_PC1_PC2.pdf")), p_pcg, width = 7, height = 5)

  # --- Plot 4: guided PC cumulative variance ---
  n_show <- min(10L, n_guided_pc)
  cum_df <- tibble(
    gPC     = factor(seq_len(n_show)),
    var_cum = cumsum(out$varPCg)[seq_len(n_show)]
  )

  p_cum <- ggplot(cum_df, aes(x = gPC, y = var_cum, group = 1)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      title = paste0(cohort_name, ": guided PC cumulative variance explained"),
      x     = "Guided PC",
      y     = "Cumulative proportion of variance"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_guided_cumulative_variance.png")), p_cum, width = 7, height = 4, dpi = 300)
  ggsave(file.path(qc_dir, paste0(prefix, "_gPCA_guided_cumulative_variance.pdf")), p_cum, width = 7, height = 4)

  invisible(list(out = out, summary = summary_row))
}

# -----------------------------
# Load + run
# -----------------------------
camp <- load_camp_qc()
cra  <- load_cra_qc()

camp_res <- run_gpca(camp$log, camp$pheno, "CAMP", "camp")
cra_res  <- run_gpca(cra$log,  cra$pheno,  "CRA",  "cra")

# --- Combined summary ---
combined_summary <- bind_rows(camp_res$summary, cra_res$summary)
write_csv(combined_summary, file.path(qc_dir, "gPCA_summary_combined.csv"))

cat("\n--- gPCA QC complete ---\n")
cat("Outputs written to:", qc_dir, "\n\n")
print(combined_summary, width = 120)
