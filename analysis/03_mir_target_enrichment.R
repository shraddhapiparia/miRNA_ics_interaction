# Target search and pathway enrichment for hsa-miR-584-5p
# Outputs are written to: results/enrichment/

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
})

mir_name <- "hsa-miR-584-5p"
target_gene_check <- "SH3TC2"

# 1. Experimentally validated targets
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

write.table(
  validated_targets,
  file.path(enrich_dir, "mir584_validated_targets.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.csv(
  validated_evidence_table,
  file.path(enrich_dir, "mir584_validated_evidence_table.csv"),
  row.names = FALSE
)

validated_source_summary <- validated_evidence_table %>%
  count(database, sort = TRUE)

write.csv(
  validated_source_summary,
  file.path(enrich_dir, "mir584_validated_source_summary.csv"),
  row.names = FALSE
)

# 2. High-confidence miRTarBase-only sensitivity target list
validated_targets_strict <- df_valid %>%
  filter(!is.na(target_symbol)) %>%
  filter(tolower(database) == "mirtarbase") %>%
  filter(grepl("reporter assay|western blot|qrt-pcr|qpcr", experiment, ignore.case = TRUE)) %>%
  pull(target_symbol) %>%
  unique() %>%
  sort()

cat("Strict validated targets:", length(validated_targets_strict), "\n")

write.table(
  validated_targets_strict,
  file.path(enrich_dir, "mir584_validated_targets_strict.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 3. Predicted targets only for SH3TC2 check
res_pred <- get_multimir(
  mirna = mir_name,
  table = "predicted",
  summary = FALSE
)

df_pred <- as_tibble(res_pred@data)

predicted_targets <- df_pred %>%
  filter(!is.na(target_symbol)) %>%
  pull(target_symbol) %>%
  unique() %>%
  sort()

cat("Total predicted unique targets:", length(predicted_targets), "\n")

sh3tc2_predicted_hits <- df_pred %>%
  filter(target_symbol == target_gene_check) %>%
  distinct()

sh3tc2_validated_hits <- df_valid %>%
  filter(target_symbol == target_gene_check) %>%
  distinct()

cat("Is", target_gene_check, "a predicted target?", nrow(sh3tc2_predicted_hits) > 0, "\n")
cat("Is", target_gene_check, "a validated target?", nrow(sh3tc2_validated_hits) > 0, "\n")

# write.csv(
#   sh3tc2_predicted_hits,
#   file.path(enrich_dir, "mir584_SH3TC2_predicted_hits.csv"),
#   row.names = FALSE
# )

# write.csv(
#   sh3tc2_validated_hits,
#   file.path(enrich_dir, "mir584_SH3TC2_validated_hits.csv"),
#   row.names = FALSE
# )

# 4. SYMBOL to ENTREZ mapping
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
  file.path(enrich_dir, "mir584_validated_targets_symbol_to_entrez.csv"),
  row.names = FALSE
)

unmapped_validated_targets <- setdiff(validated_targets, gene_df$SYMBOL)

write.table(
  unmapped_validated_targets,
  file.path(enrich_dir, "mir584_unmapped_validated_targets.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 5. Reactome enrichment for all validated targets
reactome_result <- enrichPathway(
  gene = unique(gene_df$ENTREZID),
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

reactome_df <- as.data.frame(reactome_result) %>%
  filter(Count >= 3) %>%
  arrange(p.adjust)

write.csv(
  reactome_df,
  file.path(enrich_dir, "miR-584-5p_reactome_enrichment_validated.csv"),
  row.names = FALSE
)

cat("Reactome terms with Count >= 3:", nrow(reactome_df), "\n")

if (nrow(reactome_df) > 0) {
  reactome_result_filtered <- reactome_result
  reactome_result_filtered@result <- reactome_result_filtered@result %>%
    filter(Count >= 3) %>%
    arrange(p.adjust)

  p_reactome <- dotplot(
    reactome_result_filtered,
    showCategory = 20,
    x = "GeneRatio",
    color = "p.adjust",
    size = "Count"
  ) +
    ggtitle("Reactome enrichment for validated miR-584-5p targets") +
    theme_minimal()

  ggsave(
    file.path(enrich_dir, "PathwayEnrichment_miR584_validated.png"),
    plot = p_reactome,
    width = 8,
    height = 10,
    dpi = 600
  )
}

# 6. Targeted GO glucocorticoid-response lookup
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

go_all_df <- as.data.frame(go_result_all)

gc_hits <- go_all_df %>%
  filter(ID == "GO:0051384")

write.csv(
  gc_hits,
  file.path(enrich_dir, "mir584_glucocorticoid_GO0051384_overlap.csv"),
  row.names = FALSE
)

if (nrow(gc_hits) > 0) {
  gc_overlap_genes <- strsplit(gc_hits$geneID[1], "/")[[1]]

  write.table(
    gc_overlap_genes,
    file.path(enrich_dir, "mir584_glucocorticoid_overlap_genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  cat("GO:0051384 overlap genes:", length(gc_overlap_genes), "\n")
  print(gc_overlap_genes)
} else {
  gc_overlap_genes <- character(0)
  cat("No overlap detected for GO:0051384.\n")
}

# 7. Reactome pathway context for glucocorticoid-overlap genes
if (length(gc_overlap_genes) > 0) {
  gc_gene_df <- bitr(
    gc_overlap_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, ENTREZID)

  write.csv(
    gc_gene_df,
    file.path(enrich_dir, "mir584_glucocorticoid_overlap_symbol_to_entrez.csv"),
    row.names = FALSE
  )

  gc_reactome_result <- enrichPathway(
    gene = unique(gc_gene_df$ENTREZID),
    organism = "human",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = TRUE
  )

  gc_reactome_df <- as.data.frame(gc_reactome_result) %>%
    filter(Count >= 2) %>%
    arrange(p.adjust)

  write.csv(
    gc_reactome_df,
    file.path(enrich_dir, "mir584_glucocorticoid_overlap_genes_reactome.csv"),
    row.names = FALSE
  )

  cat("Reactome terms for GC-overlap genes with Count >= 2:", nrow(gc_reactome_df), "\n")

  # Production-ready plot: top non-redundant GC-overlap Reactome pathways

  pathways_to_show <- gc_reactome_df %>%
    filter(p.adjust < 0.05) %>%
    filter(Count >= 2) %>%
    filter(Description %in% c(
      "Signaling by Interleukins",
      "Interleukin-4 and Interleukin-13 signaling",
      "MAPK family signaling cascades",
      "MAPK3 (ERK1) activation",
      "Extra-nuclear estrogen signaling",
      "ESR-mediated signaling",
      "Erythropoietin activates Phosphoinositide-3-kinase (PI3K)",
      "Interleukin receptor SHC signaling",
      "RAF activation",
      "PI3K/AKT Signaling in Cancer"
    )) %>%
    arrange(p.adjust)

  write.csv(
    pathways_to_show,
    file.path(enrich_dir, "mir584_GC_overlap_reactome_selected_pathways.csv"),
    row.names = FALSE
  )

  gc_reactome_selected <- gc_reactome_result
  gc_reactome_selected@result <- pathways_to_show

  p_gc_reactome_selected <- dotplot(
    gc_reactome_selected,
    showCategory = nrow(pathways_to_show),
    x = "GeneRatio",
    color = "p.adjust",
    size = "Count"
  ) +
    ggtitle("Reactome pathways for glucocorticoid-overlap miR-584-5p targets") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )

  ggsave(
    file.path(enrich_dir, "PathwayEnrichment_miR584_GC_overlap_selected.png"),
    plot = p_gc_reactome_selected,
    width = 9,
    height = 6,
    dpi = 600
  )

  ggsave(
    file.path(enrich_dir, "PathwayEnrichment_miR584_GC_overlap_selected.pdf"),
    plot = p_gc_reactome_selected,
    width = 9,
    height = 6
  )
}

cat("miR target enrichment analysis complete. Outputs saved to:", enrich_dir, "\n")