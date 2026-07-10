# ============================================================
# Unified GC-overlap Reactome enrichment for miR-584-5p and miR-205-5p
# Outputs: results/enrichment_unified/
# ============================================================

if (!exists("out_dir")) out_dir <- "results"

enrich_dir <- file.path(out_dir, "enrichment")
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(multiMiR)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

# -----------------------------
# Helper functions
# -----------------------------

clean_filename <- function(x) {
  x |>
    stringr::str_replace_all("hsa-", "") |>
    stringr::str_replace_all("-", "_")
}

plot_reactome_dot <- function(plot_df, title_text) {
  ggplot(
    plot_df,
    aes(
      x = neg_log10_padj,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(alpha = 0.85) +
    scale_color_continuous(trans = "reverse") +
    labs(
      title = title_text,
      x = expression(-log[10]("BH-adjusted p-value")),
      y = NULL,
      color = "p.adjust",
      size = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
}

run_mir_gc_reactome <- function(
  mir_name,
  top_n_pathways = 15,
  count_cutoff = 2,
  fdr_cutoff = 0.05
) {

  mir_short <- clean_filename(mir_name)

  cat("\n============================================================\n")
  cat("Running:", mir_name, "\n")
  cat("============================================================\n")

  # -----------------------------
  # 1. Experimentally validated targets
  # -----------------------------

  res_valid <- get_multimir(
    mirna = mir_name,
    table = "validated",
    summary = FALSE
  )

  df_valid <- as_tibble(res_valid@data)

  validated_evidence_table <- df_valid %>%
    filter(!is.na(target_symbol)) %>%
    filter(tolower(database) %in% c("mirtarbase", "tarbase", "mirecords")) %>%
    distinct()

  validated_targets <- validated_evidence_table %>%
    pull(target_symbol) %>%
    unique() %>%
    sort()

  cat("Total validated unique targets:", length(validated_targets), "\n")

  write.csv(
    validated_evidence_table,
    file.path(enrich_dir, paste0(mir_short, "_validated_evidence_table.csv")),
    row.names = FALSE
  )

  write.table(
    validated_targets,
    file.path(enrich_dir, paste0(mir_short, "_validated_targets.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  validated_source_summary <- validated_evidence_table %>%
    count(database, sort = TRUE)

  write.csv(
    validated_source_summary,
    file.path(enrich_dir, paste0(mir_short, "_validated_source_summary.csv")),
    row.names = FALSE
  )

  if (length(validated_targets) == 0) {
    warning("No validated targets found for ", mir_name)
    return(NULL)
  }

  # -----------------------------
  # 2. SYMBOL to ENTREZ mapping
  # -----------------------------

  gene_df <- bitr(
    validated_targets,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, ENTREZID)

  cat("Mapped validated target genes:", nrow(gene_df), "\n")

  write.csv(
    gene_df,
    file.path(enrich_dir, paste0(mir_short, "_validated_targets_symbol_to_entrez.csv")),
    row.names = FALSE
  )

  unmapped_validated_targets <- setdiff(validated_targets, gene_df$SYMBOL)

  write.table(
    unmapped_validated_targets,
    file.path(enrich_dir, paste0(mir_short, "_unmapped_validated_targets.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  if (nrow(gene_df) == 0) {
    warning("No validated targets mapped to Entrez IDs for ", mir_name)
    return(NULL)
  }

  # -----------------------------
  # 3. GO BP lookup for response to glucocorticoid
  # -----------------------------

  go_result_all <- enrichGO(
    gene = unique(gene_df$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    readable = TRUE,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )

  go_all_df <- as.data.frame(go_result_all) %>%
    arrange(p.adjust)

  write.csv(
    go_all_df,
    file.path(enrich_dir, paste0(mir_short, "_GO_BP_all_validated_targets.csv")),
    row.names = FALSE
  )

  gc_hit <- go_all_df %>%
    filter(ID == "GO:0051384")

  write.csv(
    gc_hit,
    file.path(enrich_dir, paste0(mir_short, "_glucocorticoid_GO0051384_overlap.csv")),
    row.names = FALSE
  )

  if (nrow(gc_hit) == 0) {
    warning("No GO:0051384 response-to-glucocorticoid hit found for ", mir_name)
    return(NULL)
  }

  gc_overlap_genes <- strsplit(gc_hit$geneID[1], "/")[[1]]

  cat("GO:0051384 overlap genes:", length(gc_overlap_genes), "\n")
  print(gc_overlap_genes)

  write.table(
    gc_overlap_genes,
    file.path(enrich_dir, paste0(mir_short, "_glucocorticoid_overlap_genes.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # -----------------------------
  # 4. Map GC-overlap genes to Entrez
  # -----------------------------

  gc_gene_df <- bitr(
    gc_overlap_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, ENTREZID)

  cat("Mapped GC-overlap genes:", nrow(gc_gene_df), "\n")

  write.csv(
    gc_gene_df,
    file.path(enrich_dir, paste0(mir_short, "_glucocorticoid_overlap_symbol_to_entrez.csv")),
    row.names = FALSE
  )

  if (nrow(gc_gene_df) == 0) {
    warning("No GC-overlap genes mapped to Entrez IDs for ", mir_name)
    return(NULL)
  }

  # -----------------------------
  # 5. Reactome enrichment on GC-overlap genes
  # -----------------------------

  gc_reactome_result <- enrichPathway(
    gene = unique(gc_gene_df$ENTREZID),
    organism = "human",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = TRUE
  )

  gc_reactome_df <- as.data.frame(gc_reactome_result) %>%
    filter(Count >= count_cutoff) %>%
    arrange(p.adjust)

  write.csv(
    gc_reactome_df,
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_genes_reactome.csv")),
    row.names = FALSE
  )

  gc_reactome_sig_df <- gc_reactome_df %>%
    filter(p.adjust < fdr_cutoff)

  write.csv(
    gc_reactome_sig_df,
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_genes_reactome_FDR05.csv")),
    row.names = FALSE
  )

  cat("Reactome terms with Count >=", count_cutoff, ":", nrow(gc_reactome_df), "\n")
  cat("Reactome terms with FDR <", fdr_cutoff, ":", nrow(gc_reactome_sig_df), "\n")

  if (nrow(gc_reactome_sig_df) == 0) {
    warning("No significant GC-overlap Reactome pathways for ", mir_name)
    return(list(
      mir_name = mir_name,
      validated_targets = validated_targets,
      gc_overlap_genes = gc_overlap_genes,
      gc_reactome_df = gc_reactome_df,
      gc_reactome_sig_df = gc_reactome_sig_df
    ))
  }

  # -----------------------------
  # 6. Plot ALL significant GC-overlap Reactome pathways
  # -----------------------------

  all_plot_df <- gc_reactome_sig_df %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      Description = fct_reorder(Description, neg_log10_padj)
    )

  p_all <- plot_reactome_dot(
    all_plot_df,
    paste0("Reactome pathways for GC-overlap ", mir_name, " targets")
  )

  ggsave(
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_reactome_all_pathways_FDR05.png")),
    plot = p_all,
    width = 10,
    height = max(6, 0.25 * nrow(all_plot_df)),
    dpi = 600,
    limitsize = FALSE
  )

  ggsave(
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_reactome_all_pathways_FDR05.pdf")),
    plot = p_all,
    width = 10,
    height = max(6, 0.25 * nrow(all_plot_df)),
    limitsize = FALSE
  )

  # -----------------------------
  # 7. Plot TOP N significant GC-overlap Reactome pathways
  # Same rule for both miRNAs
  # -----------------------------

  top_plot_df <- gc_reactome_sig_df %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n_pathways) %>%
  mutate(
    neg_log10_padj = -log10(p.adjust),
    Description = stringr::str_wrap(Description, width = 50),
    Description = forcats::fct_reorder(Description, neg_log10_padj)
  )

  write.csv(
    top_plot_df,
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_reactome_top", top_n_pathways, "_pathways.csv")),
    row.names = FALSE
  )

  p_top <- plot_reactome_dot(
    top_plot_df,
    paste0("Top ", top_n_pathways, " Reactome pathways for GC-overlap ", mir_name, " targets")
  )

  ggsave(
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_reactome_top", top_n_pathways, "_pathways.png")),
    plot = p_top,
    width = 12,
    height = max(6, 0.45 * nrow(top_plot_df)),
    dpi = 600,
    limitsize = FALSE
  )

  ggsave(
    file.path(enrich_dir, paste0(mir_short, "_GC_overlap_reactome_top", top_n_pathways, "_pathways.pdf")),
    plot = p_top,
    width = 12,
    height = max(6, 0.45 * nrow(top_plot_df)),
    limitsize = FALSE
  )

  return(list(
    mir_name = mir_name,
    validated_targets = validated_targets,
    gc_overlap_genes = gc_overlap_genes,
    gc_gene_df = gc_gene_df,
    gc_reactome_df = gc_reactome_df,
    gc_reactome_sig_df = gc_reactome_sig_df,
    top_plot_df = top_plot_df
  ))
}

# -----------------------------
# Run both miRNAs
# -----------------------------

mir_list <- c(
  "hsa-miR-584-5p",
  "hsa-miR-205-5p"
)

results_list <- lapply(
  mir_list,
  run_mir_gc_reactome,
  top_n_pathways = 15,
  count_cutoff = 2,
  fdr_cutoff = 0.05
)

names(results_list) <- mir_list

# -----------------------------
# Combined summary table
# -----------------------------

summary_df <- bind_rows(lapply(results_list, function(x) {
  if (is.null(x)) return(NULL)

  tibble(
    miRNA = x$mir_name,
    n_validated_targets = length(x$validated_targets),
    n_GC_overlap_genes = length(x$gc_overlap_genes),
    n_GC_overlap_reactome_terms = nrow(x$gc_reactome_df),
    n_GC_overlap_reactome_FDR05 = nrow(x$gc_reactome_sig_df)
  )
}))

write.csv(
  summary_df,
  file.path(enrich_dir, "GC_overlap_reactome_summary_both_miRNAs.csv"),
  row.names = FALSE
)

print(summary_df)

cat("\nUnified GC-overlap Reactome analysis complete. Outputs saved to:", enrich_dir, "\n")