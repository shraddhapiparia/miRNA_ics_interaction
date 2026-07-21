set.seed(2025)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble)
  library(preprocessCore); library(ordinal); library(pROC)
  library(ggplot2); library(ggrepel)
})

# -----------------------------
# Inputs / outputs
# -----------------------------
camp_counts_path <- "data/CAMP/CAMP.492.COUNT.csv"
camp_pheno_path  <- "data/CAMP/camp_pheno_2022-08-04.tsv"
camp_link_path   <- "data/CAMP/shortened_phenotype_files_ver5.csv"

cra_counts_path     <- "data/CRA/2022-08-16_CRA_1159_exceRpt_miRNA_ReadCountscopy.csv"
cra_link_path       <- "data/CRA/cra_shortened_phenotype_files10-11.csv"
cra_hosp_pheno_path <- "data/CRA/COS_TRIO_pheno_1165.csv"
cra_er_pheno_path   <- "data/CRA/PhenotypeInfo.csv"

camp_new_counts <- "data/CAMP_NEW/New_CAMP_data_raw_read_count.csv"
camp_new_pheno  <- "data/CAMP_NEW/final_meta_all_exacerbation_phenotypes.csv"

out_dir <- if (exists("out_dir")) out_dir else "results"
dir.create(file.path(out_dir, "models"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "meta"),    recursive = TRUE, showWarnings = FALSE)

min_count        <- 5
min_prop_samples <- 0.5

comparison_plan <- tibble(
  comparison    = c("EDHOS_CAMP_vs_HospAsthma_CRA", "EDHOS_CAMP_vs_ERVisits_CRA"),
  candidate_mir = c("hsa-miR-584-5p", "hsa-miR-205-5p"),
  cra_raw_col   = c("Hospitalized_Asthma_Last_Yr", "ER_Visits_Asthma_Last_Yr"),
  cra_outcome   = c("CRA_Hospitalized_Asthma_Last_Yr_binary", "CRA_ER_Visits_Asthma_Last_Yr_binary")
)

candidate_mirs   <- c("hsa-miR-584-5p", "hsa-miR-205-5p")
camp_main_covars <- c("AGE.x", "SEX.x", "RACE.x")
cra_main_covars  <- c("age.x", "gender.x")
camp_auc_covars  <- c("PREFEVPP_RZ", "TOT.EOS_S3", "LOG10IGE_S3", "AGE.x", "SEX.x")
cra_auc_covars   <- c("pctpred_fev1_pre_BD", "log10eos", "log10Ige", "age.x", "gender.x")

# -----------------------------
# Helpers
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

or_ci <- function(beta, se) {
  ok <- !is.na(beta) & !is.na(se)
  tibble(
    OR       = ifelse(ok, exp(beta), NA_real_),
    lower_CI = ifelse(ok, exp(beta - 1.96 * se), NA_real_),
    upper_CI = ifelse(ok, exp(beta + 1.96 * se), NA_real_),
    OR_CI    = ifelse(ok,
      paste0(round(exp(beta), 2), " (",
             round(exp(beta - 1.96 * se), 2), "-",
             round(exp(beta + 1.96 * se), 2), ")"),
      NA_character_)
  )
}

empty_coef <- function(term) tibble(term = term, beta = NA_real_, se = NA_real_, statistic = NA_real_, p_value = NA_real_)

coef_term <- function(fit, term) {
  z <- tryCatch(as.data.frame(coef(summary(fit))) |> rownames_to_column("term"), error = function(e) NULL)
  if (is.null(z) || !term %in% z$term) return(empty_coef(term))
  p_col    <- grep("^Pr\\(", names(z), value = TRUE)[1]
  stat_col <- grep("z value|t value", names(z), value = TRUE)[1]
  z |> filter(.data$term == !!term) |>
    transmute(term, beta = Estimate, se = `Std. Error`, statistic = .data[[stat_col]], p_value = .data[[p_col]])
}

coef_clm <- function(fit, term) {
  z <- tryCatch(as.data.frame(coef(summary(fit))) |> rownames_to_column("term"), error = function(e) NULL)
  if (is.null(z) || !term %in% z$term) return(empty_coef(term))
  z |> filter(!grepl("\\|", .data$term), .data$term == !!term) |>
    transmute(term, beta = Estimate, se = `Std. Error`, statistic = `z value`, p_value = `Pr(>|z|)`)
}

clean_covars <- function(dat, covars) {
  covars <- intersect(covars, names(dat))
  covars[map_lgl(covars, function(v) {
    x <- dat[[v]]
    if (is.numeric(x)) length(unique(na.omit(x))) > 1 && sd(x, na.rm = TRUE) > 0
    else length(unique(na.omit(as.character(x)))) > 1
  })]
}

prep_model_dat <- function(pheno_df, log_mat, mir, cols) {
  if (!mir %in% rownames(log_mat)) stop(mir, " not found in miRNA matrix after filtering")
  pheno_df$mir_expr <- as.numeric(log_mat[mir, rownames(pheno_df)])
  pheno_df |> select(all_of(c(cols, "mir_expr"))) |> drop_na() |> droplevels()
}

rhs <- function(parts) paste(parts[nzchar(parts)], collapse = " + ")

safe_glm_coef <- function(dat, outcome, rhs_txt, term) {
  if (nrow(dat) < 10 || length(unique(dat[[outcome]])) < 2) return(empty_coef(term))
  fit <- tryCatch(glm(as.formula(paste(outcome, "~", rhs_txt)), data = dat, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) empty_coef(term) else coef_term(fit, term)
}
safe_lm_coef <- function(dat, outcome, rhs_txt, term) {
  if (nrow(dat) < 10 || length(unique(dat[[outcome]])) < 2) return(empty_coef(term))
  fit <- tryCatch(lm(as.formula(paste(outcome, "~", rhs_txt)), data = dat), error = function(e) NULL)
  if (is.null(fit)) empty_coef(term) else coef_term(fit, term)
}
safe_clm_coef <- function(dat, outcome, rhs_txt, term) {
  if (nrow(dat) < 10 || length(unique(dat[[outcome]])) < 2) return(empty_coef(term))
  fit <- tryCatch(ordinal::clm(as.formula(paste(outcome, "~", rhs_txt)), data = dat, link = "logit"), error = function(e) NULL)
  if (is.null(fit)) empty_coef(term) else coef_clm(fit, term)
}

fdr_by_model <- function(x) {
  x |> group_by(comparison, cohort, outcome, analysis, model) |>
    mutate(FDR = p.adjust(p_value, "fdr")) |> ungroup()
}

find_existing_file <- function(paths) paths[file.exists(paths)][1]

label_cohort_outcome <- function(cohort, outcome) {
  case_when(
    cohort == "CAMP"                                            ~ "CAMP EDHOS",
    cohort == "CRA" & grepl("Hosp", outcome, ignore.case = TRUE) ~ "CRA hospitalization",
    cohort == "CRA" & grepl("ER",   outcome, ignore.case = TRUE) ~ "CRA ER visits",
    TRUE ~ paste(cohort, outcome)
  )
}

fmt_p <- function(p) ifelse(p < 0.001, "p<0.001", paste0("p=", round(p, 3)))

# -----------------------------
# Load data
# -----------------------------
load_camp <- function() {
  counts <- read.csv(camp_counts_path, check.names = FALSE)
  pheno  <- read.table(camp_pheno_path, sep = "\t", header = TRUE, check.names = FALSE)
  link   <- read.csv(camp_link_path, check.names = FALSE)
  mat <- filter_mirnas(counts, min_count, min_prop_samples)
  pheno <- pheno[pheno$TG != 2, , drop = FALSE] |> merge(link, by.x = "camp", by.y = "CAMP_ID", all = FALSE)
  rownames(pheno) <- pheno$`row.names(CAMP_492_COUNT.Trans)`
  pheno <- pheno |> filter(!is.na(EDHOS_cum_Y1)) |>
    mutate(
      EDHOS_cum_Y1_binary  = factor(ifelse(EDHOS_cum_Y1 <= 1, 1, 2), levels = c(1, 2)),
      EDHOS_cum_Y1_ordered = factor(pmin(EDHOS_cum_Y1, 3), levels = 0:3, ordered = TRUE),
      treatment = relevel(factor(TG), ref = "3"),
      SEX.x = factor(SEX.x), RACE.x = factor(RACE.x)
    )
  common <- intersect(colnames(mat), rownames(pheno))
  list(log = norm_quantile_log(mat[, common, drop = FALSE]), pheno = pheno[common, , drop = FALSE])
}

load_cra <- function(raw_col, preferred_paths) {
  pheno_path <- find_existing_file(preferred_paths)
  if (is.na(pheno_path)) stop("No CRA phenotype file found among: ", paste(preferred_paths, collapse = ", "))
  counts <- read.csv(cra_counts_path, check.names = FALSE)
  pheno  <- read.csv(pheno_path, check.names = FALSE)
  link   <- read.csv(cra_link_path, check.names = FALSE)
  if (!raw_col %in% names(pheno)) stop(raw_col, " missing from ", pheno_path)
  mat <- filter_mirnas(counts, min_count, min_prop_samples)
  pheno <- merge(pheno, link, by.x = "S_SUBJECTID", by.y = "ST.ID")
  pheno$sample.id <- gsub("-", ".", pheno$sample.id)
  rownames(pheno) <- pheno$sample.id
  pheno <- pheno |> filter(!is.na(.data[[raw_col]])) |>
    mutate(
      outcome_binary = factor(.data[[raw_col]], levels = c(1, 2)),
      treatment      = relevel(factor(Inhaled_Steroids), ref = "1"),
      gender.x       = factor(gender.x)
    )
  if (any(is.na(pheno$outcome_binary))) stop(raw_col, " contains values outside expected binary levels 1/2")
  common <- intersect(colnames(mat), rownames(pheno))
  list(log = norm_quantile_log(mat[, common, drop = FALSE]), pheno = pheno[common, , drop = FALSE], pheno_path = pheno_path)
}

load_camp_new <- function() {
  counts  <- read.csv(camp_new_counts, check.names = FALSE)
  pheno   <- read.csv(camp_new_pheno,  check.names = FALSE)
  tg_vals <- unique(pheno$TG[!is.na(pheno$TG)])
  if (!all(c(1L, 3L) %in% tg_vals))
    warning("CAMP_NEW: expected TG==1 (ICS) and TG==3 (non-ICS); found: ", paste(sort(tg_vals), collapse = ", "))
  cat("\nCAMP_NEW: phenotype samples:       ", nrow(pheno), "\n")
  cat(  "CAMP_NEW: miRNAs before filtering: ", nrow(counts), "\n")
  pheno <- pheno[!is.na(pheno$TG) & pheno$TG != 2, , drop = FALSE]
  rownames(pheno) <- pheno$`LL ID`
  mat_filtered <- filter_mirnas(counts, min_count, min_prop_samples) # remove filteration to check 205-5p replication
  # Full matrix
  mat_all <- as.matrix(counts[, -1, drop = FALSE])
  rownames(mat_all) <- counts[[1]]
  storage.mode(mat_all) <- "numeric"
  candidate_present <- intersect(candidate_mirs, rownames(mat_all))
  mat <- mat_all[union(rownames(mat_filtered), candidate_present),,drop = FALSE]
  cat(  "CAMP_NEW: miRNAs after filtering:  ", nrow(mat), "\n")
  common <- intersect(colnames(mat), pheno$`LL ID`)
  cat(  "CAMP_NEW: matched samples:         ", length(common), "\n")
  if (length(common) < 20)
    warning("CAMP_NEW: unexpectedly low sample overlap (", length(common), " samples matched)")
  pheno <- pheno[common, , drop = FALSE] |>
    filter(!is.na(EDHOS_cum_Y1)) |>
    mutate(
      EDHOS_cum_Y1_binary  = factor(ifelse(EDHOS_cum_Y1 < 1, 1, 2), levels = c(1, 2)),
      EDHOS_cum_Y1_ordered = factor(pmin(EDHOS_cum_Y1, 3), levels = 0:3, ordered = TRUE),
      treatment            = relevel(factor(TG), ref = "3"),
      SEX                  = factor(SEX),
      RACE                 = factor(RACE)
    )

  
x <- as.numeric(counts[counts[[1]] == "hsa-miR-205-5p", -1])

cat("Samples with count >5:", sum(x > 5), "/", length(x), "\n")
cat("Median count:", median(x), "\n")
cat("Max count:", max(x), "\n")

  if (length(unique(na.omit(as.character(pheno$EDHOS_cum_Y1_binary)))) < 1)
    warning("CAMP_NEW: binary outcome has only one class after filtering")
  cat("CAMP_NEW: treatment counts:\n");      print(table(pheno$TG,                  useNA = "ifany"))
  cat("CAMP_NEW: EDHOS_cum_Y1 counts:\n");   print(table(pheno$EDHOS_cum_Y1,        useNA = "ifany"))
  cat("CAMP_NEW: binary outcome counts:\n"); print(table(pheno$EDHOS_cum_Y1_binary, useNA = "ifany"))
  list(log = norm_quantile_log(mat[, common, drop = FALSE]), pheno = pheno)
}

camp     <- load_camp()
cra_hosp <- load_cra("Hospitalized_Asthma_Last_Yr", c(cra_hosp_pheno_path, cra_er_pheno_path))
cra_er   <- load_cra("ER_Visits_Asthma_Last_Yr",    c(cra_er_pheno_path, cra_hosp_pheno_path))
camp_new <- load_camp_new()

cat("\nCAMP EDHOS counts:\n");                    print(table(camp$pheno$EDHOS_cum_Y1, useNA = "ifany"))
cat("\nCAMP EDHOS collapsed 0/1/2/3+ counts:\n"); print(table(camp$pheno$EDHOS_cum_Y1_ordered, useNA = "ifany"))
cat("\nCRA hospitalization phenotype file: ", cra_hosp$pheno_path, "\n", sep = ""); print(table(cra_hosp$pheno$outcome_binary, useNA = "ifany"))
cat("\nCRA ER visits phenotype file: ",       cra_er$pheno_path,   "\n", sep = ""); print(table(cra_er$pheno$outcome_binary,   useNA = "ifany"))

# -----------------------------
# Binary interaction + stratified logistic scans
# -----------------------------
run_binary_scan <- function(log_mat, pheno_df, outcome_col, treatment_col, inter_term,
                            covars, comparison, cohort, outcome_label, trt_nonics, trt_ics) {
  out <- map_dfr(rownames(log_mat), function(mir) {
    dat <- prep_model_dat(pheno_df, log_mat, mir, c(outcome_col, treatment_col, covars))
    cv  <- clean_covars(dat, covars)
    int <- safe_glm_coef(dat, outcome_col, rhs(c(paste0(treatment_col, " * mir_expr"), cv)), inter_term) |>
      mutate(model = "Interaction", n = nrow(dat), treatment_level = NA_character_)
    str <- map_dfr(c(trt_nonics, trt_ics), function(g) {
      sub <- dat |> filter(as.character(.data[[treatment_col]]) == g) |> droplevels()
      safe_glm_coef(sub, outcome_col, rhs(c("mir_expr", clean_covars(sub, covars))), "mir_expr") |>
        mutate(model = ifelse(g == trt_ics, "Stratified: ICS", "Stratified: non-ICS"),
               n = nrow(sub), treatment_level = g)
    })
    bind_rows(int, str) |> mutate(miRNA = mir)
  }) |>
    mutate(comparison = comparison, cohort = cohort, outcome = outcome_label, analysis = "Binary logistic") |>
    fdr_by_model()
  bind_cols(out, or_ci(out$beta, out$se))
}

camp_binary     <- run_binary_scan(
  camp$log, camp$pheno, "EDHOS_cum_Y1_binary", "treatment", "treatment1:mir_expr", camp_main_covars,
  "CAMP_EDHOS_all_binary", "CAMP", "CAMP_EDHOS_cum_Y1_binary", "3", "1")
camp_for_meta <- camp_binary |> mutate(comparison = "CAMP_EDHOS")
cra_hosp_binary <- run_binary_scan(
  cra_hosp$log, cra_hosp$pheno, "outcome_binary", "treatment", "treatment2:mir_expr", cra_main_covars,
  "EDHOS_CAMP_vs_HospAsthma_CRA", "CRA", "CRA_Hospitalized_Asthma_Last_Yr_binary", "1", "2")
cra_er_binary   <- run_binary_scan(
  cra_er$log, cra_er$pheno, "outcome_binary", "treatment", "treatment2:mir_expr", cra_main_covars,
  "EDHOS_CAMP_vs_ERVisits_CRA", "CRA", "CRA_ER_Visits_Asthma_Last_Yr_binary", "1", "2")

# -----------------------------
# Candidate main results table
# -----------------------------
camp_main <- comparison_plan |> pmap_dfr(function(comparison, candidate_mir, cra_raw_col, cra_outcome) {
  camp_binary |> filter(miRNA == .env$candidate_mir) |> mutate(comparison = .env$comparison)
})
main_results <- bind_rows(
  camp_main,
  cra_hosp_binary |> filter(miRNA == "hsa-miR-584-5p"),
  cra_er_binary   |> filter(miRNA == "hsa-miR-205-5p")
) |>
  select(comparison, miRNA, cohort, outcome, analysis, model, n, beta, se, OR, lower_CI, upper_CI, OR_CI, p_value, FDR)
write_csv(main_results, file.path(out_dir, "models", "candidate_main_binary_interaction_stratified_results.csv"))

# -----------------------------
# CAMP EDHOS ordinal + linear sensitivity; no CRA linear/ordinal
# -----------------------------
run_camp_sensitivity <- function(log_mat, pheno_df) {
  map_dfr(rownames(log_mat), function(mir) {
    dat_lm  <- prep_model_dat(pheno_df, log_mat, mir, c("EDHOS_cum_Y1", "treatment", camp_main_covars))
    dat_ord <- prep_model_dat(pheno_df, log_mat, mir, c("EDHOS_cum_Y1_ordered", "treatment", "AGE.x", "SEX.x"))
    lin <- bind_rows(
      safe_lm_coef(dat_lm, "EDHOS_cum_Y1",
                   rhs(c("treatment * mir_expr", clean_covars(dat_lm, camp_main_covars))), "treatment1:mir_expr") |>
        mutate(analysis = "Linear EDHOS", model = "Interaction", n = nrow(dat_lm)),
      map_dfr(c("3", "1"), function(g) {
        sub <- dat_lm |> filter(as.character(treatment) == g) |> droplevels()
        safe_lm_coef(sub, "EDHOS_cum_Y1", rhs(c("mir_expr", clean_covars(sub, camp_main_covars))), "mir_expr") |>
          mutate(analysis = "Linear EDHOS",
                 model = ifelse(g == "1", "Stratified: ICS", "Stratified: non-ICS"), n = nrow(sub))
      })
    )
    ord <- bind_rows(
      safe_clm_coef(dat_ord, "EDHOS_cum_Y1_ordered",
                    rhs(c("treatment * mir_expr", clean_covars(dat_ord, c("AGE.x", "SEX.x")))), "treatment1:mir_expr") |>
        mutate(analysis = "Ordinal EDHOS", model = "Interaction", n = nrow(dat_ord)),
      map_dfr(c("3", "1"), function(g) {
        sub <- dat_ord |> filter(as.character(treatment) == g) |> droplevels()
        safe_clm_coef(sub, "EDHOS_cum_Y1_ordered",
                      rhs(c("mir_expr", clean_covars(sub, c("AGE.x", "SEX.x")))), "mir_expr") |>
          mutate(analysis = "Ordinal EDHOS",
                 model = ifelse(g == "1", "Stratified: ICS", "Stratified: non-ICS"), n = nrow(sub))
      })
    )
    bind_rows(lin, ord) |> mutate(miRNA = mir)
  }) |>
    mutate(comparison = "CAMP_EDHOS_sensitivity", cohort = "CAMP", outcome = "EDHOS_cum_Y1") |>
    fdr_by_model() |>
    mutate(
      OR       = ifelse(analysis == "Ordinal EDHOS", exp(beta), NA_real_),
      lower_CI = ifelse(analysis == "Ordinal EDHOS", exp(beta - 1.96 * se), NA_real_),
      upper_CI = ifelse(analysis == "Ordinal EDHOS", exp(beta + 1.96 * se), NA_real_),
      OR_CI    = ifelse(analysis == "Ordinal EDHOS", or_ci(beta, se)$OR_CI, NA_character_)
    )
}

camp_sensitivity <- run_camp_sensitivity(camp$log, camp$pheno)
camp_sensitivity_candidates <- comparison_plan |> pmap_dfr(function(comparison, candidate_mir, cra_raw_col, cra_outcome) {
  camp_sensitivity |> filter(miRNA == .env$candidate_mir) |> mutate(comparison = .env$comparison)
}) |>
  select(comparison, miRNA, cohort, outcome, analysis, model, n, beta, se, OR, lower_CI, upper_CI, OR_CI, p_value, FDR)
write_csv(camp_sensitivity_candidates, file.path(out_dir, "models", "candidate_camp_edhos_ordinal_linear_results.csv"))

# -----------------------------
# Volcano plot: CAMP-wide interaction scan
# -----------------------------
volcano_df <- camp_binary |>
  filter(model == "Interaction", !is.na(p_value)) |>
  mutate(
    neg_log10_p  = -log10(p_value),
    is_candidate = miRNA %in% candidate_mirs,
    point_group  = case_when(
      is_candidate    ~ "candidate",
      p_value < 0.05  ~ "significant",
      TRUE            ~ "ns"
    ),
    mir_label = ifelse(is_candidate, gsub("hsa-", "", miRNA), NA_character_)
  )

p_volcano <- ggplot(
  volcano_df,
  aes(x = beta, y = neg_log10_p, color = point_group,
      shape = point_group, size = point_group)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_point(alpha = 0.85) +
  # candidate labels - cleaner style, white fill, thin red border
  geom_label_repel(
    data        = volcano_df |> filter(is_candidate),
    aes(label   = mir_label),
    color       = "#E31A1C", fill = "white",
    size        = 3.3, fontface = "bold",
    label.size  = 0.3,
    box.padding = 0.6, point.padding = 0.5,
    segment.color = "gray60", segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = 20, show.legend = FALSE
  ) +
  # p = 0.05 inline annotation
  annotate("text",
           x = max(volcano_df$beta, na.rm = TRUE) * 0.98,
           y = -log10(0.05) + 0.08,
           label = "p = 0.05", hjust = 1, size = 3,
           color = "gray30", fontface = "italic") +
  scale_color_manual(
    values = c(candidate = "#E31A1C", significant = "#1F78B4", ns = "gray75"),
    labels = c(candidate = "Candidate miRNA",
               significant = "p < 0.05", ns = "p ≥ 0.05")
  ) +
  scale_shape_manual(
    values = c(candidate = 18, significant = 16, ns = 16),
    labels = c(candidate = "Candidate miRNA",
               significant = "p < 0.05", ns = "p ≥ 0.05")
  ) +
  scale_size_manual(
    values = c(candidate = 4.5, significant = 2.2, ns = 1.5),
    labels = c(candidate = "Candidate miRNA",
               significant = "p < 0.05", ns = "p ≥ 0.05")
  ) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  labs(
    x = expression("Interaction "*beta*" (ICS vs. non-ICS)"),
    y = expression(-log[10]*"(p)"),
    color = NULL, shape = NULL, size = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "gray93"),
    axis.title         = element_text(face = "bold"),
    legend.text        = element_text(size = 10)
  )

ggsave(file.path(out_dir, "figures", "camp_interaction_volcano.png"), p_volcano, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "figures", "camp_interaction_volcano.pdf"), p_volcano, width = 7, height = 5.5)

# -----------------------------
# Full IVW fixed-effect meta-analysis: all miRNAs common to CAMP + CRA
# -----------------------------
run_ivw_meta <- function(combined_df) {
  combined_df |>
    filter(!is.na(beta), !is.na(se), se > 0) |>
    group_by(comparison, miRNA, model) |>
    filter(n() == 2) |>
    summarise(
      cohorts   = paste(sort(unique(cohort)), collapse = "+"),
      n_cohorts = n(),
      beta_meta = sum(beta / se^2) / sum(1 / se^2),
      se_meta   = sqrt(1 / sum(1 / se^2)),
      .groups   = "drop"
    ) |>
    mutate(
      z_meta        = beta_meta / se_meta,
      p_meta        = 2 * pnorm(-abs(z_meta)),
      OR_meta       = exp(beta_meta),
      lower_CI_meta = exp(beta_meta - 1.96 * se_meta),
      upper_CI_meta = exp(beta_meta + 1.96 * se_meta),
      OR_CI_meta    = paste0(round(OR_meta, 2), " (",
                             round(lower_CI_meta, 2), "-",
                             round(upper_CI_meta, 2), ")")
    )
}

# Relabel CAMP comparison to match each CRA comparison before stacking
camp_for_meta_hosp <- camp_binary |> mutate(comparison = "EDHOS_CAMP_vs_HospAsthma_CRA")
camp_for_meta_er   <- camp_binary |> mutate(comparison = "EDHOS_CAMP_vs_ERVisits_CRA")

meta_all <- run_ivw_meta(bind_rows(
  camp_for_meta_hosp, cra_hosp_binary,
  camp_for_meta_er,   cra_er_binary
)) |>
  group_by(comparison, model) |>
  mutate(
    FDR_meta        = p.adjust(p_meta, "BH"),
    bonferroni_meta = p.adjust(p_meta, "bonferroni")
  ) |>
  ungroup()

write_csv(meta_all, file.path(out_dir, "meta", "all_common_mirna_binary_ivw_meta_results.csv"))

meta_candidates <- meta_all |>
  filter(
    (comparison == "EDHOS_CAMP_vs_HospAsthma_CRA" & miRNA == "hsa-miR-584-5p") |
    (comparison == "EDHOS_CAMP_vs_ERVisits_CRA"   & miRNA == "hsa-miR-205-5p")
  )
write_csv(meta_candidates, file.path(out_dir, "meta", "candidate_binary_ivw_meta_results.csv"))

# -----------------------------
# Forest plots: candidate interaction and ICS-stratified effects
# -----------------------------
build_forest_df <- function(model_filter) {
  camp_rows <- camp_binary |>
    filter(miRNA %in% candidate_mirs, model == model_filter) |>
    mutate(
      comparison = case_when(
        miRNA == "hsa-miR-584-5p" ~ "EDHOS_CAMP_vs_HospAsthma_CRA",
        miRNA == "hsa-miR-205-5p" ~ "EDHOS_CAMP_vs_ERVisits_CRA"
      ),
      row_label = "CAMP (discovery)",
      row_type  = "cohort"
    )
  cra_rows <- bind_rows(
    cra_hosp_binary |> filter(miRNA == "hsa-miR-584-5p", model == model_filter) |>
      mutate(row_label = "GACRS hospitalization", row_type = "cohort"),
    cra_er_binary   |> filter(miRNA == "hsa-miR-205-5p", model == model_filter) |>
      mutate(row_label = "GACRS ER visits",       row_type = "cohort")
  )
  meta_rows <- meta_candidates |>
    filter(model == model_filter) |>
    transmute(
      comparison, miRNA,
      beta = beta_meta, se = se_meta,
      OR = OR_meta, lower_CI = lower_CI_meta, upper_CI = upper_CI_meta,
      OR_CI = OR_CI_meta, p_value = p_meta,
      n = NA_integer_, cohort = "Meta", analysis = "Binary logistic",
      FDR = FDR_meta, model = model_filter,
      row_label = "IVW meta-analysis", row_type = "meta"
    )
  bind_rows(camp_rows, cra_rows, meta_rows) |>
    mutate(
      candidate_label = case_when(
        miRNA == "hsa-miR-584-5p" ~ "miR-584-5p",
        miRNA == "hsa-miR-205-5p" ~ "miR-205-5p"
      ),
      # facet order: headline result (miR-584-5p) on top
      candidate_label = factor(candidate_label,
                               levels = c("miR-584-5p", "miR-205-5p")),
      # row order within each panel: meta at bottom, CAMP at top
      # ggplot plots y-axis bottom-to-top, so reverse the visual order you want
      row_label = factor(
        row_label,
        levels = rev(c("CAMP (discovery)",
                       "GACRS hospitalization",
                       "GACRS ER visits",
                       "IVW meta-analysis"))
      )
    )
}

draw_forest <- function(df) {
  df2 <- df |>
    group_by(candidate_label) |>
    mutate(
      label_x = max(upper_CI, na.rm = TRUE) * 1.20,
      OR_CI_label = sprintf(
        "%.2f [%.2f-%.2f]",
        OR, lower_CI, upper_CI
      ),
      p_label = paste0("p=", formatC(p_value, format = "g", digits = 2)),
      text_label = paste0(OR_CI_label, ", ", p_label)
    ) |>
    ungroup()

  x_min <- min(df2$lower_CI, na.rm = TRUE) * 0.75
  x_max <- max(df2$label_x, na.rm = TRUE) * 1.35

  ggplot(
    df2,
    aes(
      y = row_label,
      x = OR,
      xmin = lower_CI,
      xmax = upper_CI,
      color = row_type,
      shape = row_type
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "gray45",
      linewidth = 0.6
    ) +
    geom_errorbarh(
      height = 0.18,
      linewidth = 0.9
    ) +
    geom_point(
      aes(size = row_type),
      stroke = 0.8
    ) +
    geom_text(
      aes(x = label_x, label = text_label),
      hjust = 0,
      size = 4.2,
      color = "gray15",
      show.legend = FALSE
    ) +
    scale_x_log10(
      name = "Odds ratio (log scale)",
      breaks = c(0.3, 0.5, 1, 2, 3, 5, 10),
      labels = c("0.3", "0.5", "1", "2", "3", "5", "10"),
      limits = c(x_min, x_max)
    ) +
    scale_y_discrete(name = NULL) +
    scale_color_manual(
      values = c(cohort = "#2166AC", meta = "#B2182B"),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(cohort = 16, meta = 18),
      guide = "none"
    ) +
    scale_size_manual(
      values = c(cohort = 3.5, meta = 5.5),
      guide = "none"
    ) +
    facet_wrap(
      ~ candidate_label,
      ncol = 1,
      scales = "free_y"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 13, hjust = 0),
      strip.background = element_rect(fill = "gray90", color = "gray55", linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray88", linewidth = 0.5),
      axis.title.x = element_text(face = "bold", size = 15, margin = margin(t = 8)),
      axis.text.x = element_text(size = 12, color = "gray25"),
      axis.text.y = element_text(size = 12, color = "gray25"),
      panel.spacing = unit(1.1, "lines"),
      plot.margin = margin(8, 60, 8, 8)
    )
}

forest_int <- build_forest_df("Interaction")
print(forest_int |>
  group_by(candidate_label) |>
  summarise(rows = paste(row_label, collapse = " | "), .groups = "drop"))
p_forest_int <- draw_forest(forest_int)
ggsave(file.path(out_dir, "figures", "candidate_interaction_forest_plot.png"), p_forest_int, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "figures", "candidate_interaction_forest_plot.pdf"), p_forest_int, width = 12, height = 5.5)

forest_ics <- build_forest_df("Stratified: ICS")
p_forest_ics <- draw_forest(forest_ics)
ggsave(file.path(out_dir, "figures", "candidate_ics_stratified_forest_plot.png"), p_forest_ics, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "figures", "candidate_ics_stratified_forest_plot.pdf"), p_forest_ics, width = 12, height = 5.5)

forest_nonics <- build_forest_df("Stratified: non-ICS")
p_forest_nonics <- draw_forest(forest_nonics)
ggsave(file.path(out_dir, "figures", "candidate_nonics_stratified_forest_plot.png"), p_forest_nonics, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "figures", "candidate_nonics_stratified_forest_plot.pdf"), p_forest_nonics, width = 12, height = 5.5)

# -----------------------------
# AUC summary
# -----------------------------
auc_one <- function(log_mat, pheno_df, mir, outcome_col, treatment_col, covars,
                    comparison, cohort, outcome_label, trt_nonics, trt_ics) {
  dat <- prep_model_dat(pheno_df, log_mat, mir, c(outcome_col, treatment_col, covars))
  cv  <- clean_covars(dat, covars)
  fit <- glm(as.formula(paste(outcome_col, "~", rhs(c(paste0(treatment_col, " * mir_expr"), cv)))),
             data = dat, family = binomial())
  all_row <- tibble(scope = "All", n = nrow(dat),
    auc = as.numeric(pROC::auc(pROC::roc(dat[[outcome_col]], predict(fit, type = "response"), quiet = TRUE))))
  str_rows <- map_dfr(c(trt_nonics, trt_ics), function(g) {
    sub   <- dat |> filter(as.character(.data[[treatment_col]]) == g) |> droplevels()
    scope <- ifelse(g == trt_ics, "ICS", "Non-ICS")
    if (nrow(sub) < 10 || length(unique(sub[[outcome_col]])) < 2)
      return(tibble(scope = scope, n = nrow(sub), auc = NA_real_))
    fit_s <- glm(as.formula(paste(outcome_col, "~", rhs(c("mir_expr", clean_covars(sub, covars))))),
                 data = sub, family = binomial())
    tibble(scope = scope, n = nrow(sub),
           auc = as.numeric(pROC::auc(pROC::roc(sub[[outcome_col]], predict(fit_s, type = "response"), quiet = TRUE))))
  })
  bind_rows(all_row, str_rows) |>
    mutate(comparison = comparison, miRNA = mir, cohort = cohort, outcome = outcome_label, .before = 1)
}

auc_summary <- bind_rows(
  auc_one(camp$log,     camp$pheno,     "hsa-miR-584-5p", "EDHOS_cum_Y1_binary", "treatment", camp_auc_covars, "EDHOS_CAMP_vs_HospAsthma_CRA", "CAMP", "CAMP_EDHOS_cum_Y1_binary",                "3", "1"),
  auc_one(cra_hosp$log, cra_hosp$pheno, "hsa-miR-584-5p", "outcome_binary",       "treatment", cra_auc_covars,  "EDHOS_CAMP_vs_HospAsthma_CRA", "CRA",  "CRA_Hospitalized_Asthma_Last_Yr_binary",  "1", "2"),
  auc_one(camp$log,     camp$pheno,     "hsa-miR-205-5p", "EDHOS_cum_Y1_binary", "treatment", camp_auc_covars, "EDHOS_CAMP_vs_ERVisits_CRA",   "CAMP", "CAMP_EDHOS_cum_Y1_binary",                "3", "1"),
  auc_one(cra_er$log,   cra_er$pheno,   "hsa-miR-205-5p", "outcome_binary",       "treatment", cra_auc_covars,  "EDHOS_CAMP_vs_ERVisits_CRA",   "CRA",  "CRA_ER_Visits_Asthma_Last_Yr_binary",     "1", "2")
)
write_csv(auc_summary, file.path(out_dir, "models", "candidate_auc_summary_all_ics_nonics.csv"))

# -----------------------------
# Adjusted predicted probability values
# -----------------------------
predprob_one <- function(log_mat, pheno_df, mir, outcome_col, treatment_col, covars,
                         comparison, cohort, outcome_label) {
  dat <- prep_model_dat(pheno_df, log_mat, mir, c(outcome_col, treatment_col, covars))
  cv  <- clean_covars(dat, covars)
  fit <- glm(as.formula(paste(outcome_col, "~", rhs(c(paste0(treatment_col, " * mir_expr"), cv)))),
             data = dat, family = binomial())
  grid <- expand.grid(
    mir_expr = seq(quantile(dat$mir_expr, 0.1), quantile(dat$mir_expr, 0.9), length.out = 100),
    tmp      = levels(dat[[treatment_col]])
  )
  names(grid)[2] <- treatment_col
  grid[[treatment_col]] <- factor(grid[[treatment_col]], levels = levels(dat[[treatment_col]]))
  for (v in cv)
    grid[[v]] <- if (is.numeric(dat[[v]])) mean(dat[[v]], na.rm = TRUE) else
      factor(levels(dat[[v]])[1], levels = levels(dat[[v]]))
  pr <- predict(fit, newdata = grid, type = "link", se.fit = TRUE)
  grid |> mutate(
    fit       = pr$fit, se = pr$se.fit,
    pred_prob = plogis(fit),
    lower     = plogis(fit - 1.96 * se),
    upper     = plogis(fit + 1.96 * se),
    comparison = comparison, miRNA = mir, cohort = cohort, outcome = outcome_label, .before = 1
  )
}

predprob_all <- bind_rows(
  predprob_one(camp$log,     camp$pheno,     "hsa-miR-584-5p", "EDHOS_cum_Y1_binary", "treatment", camp_main_covars, "EDHOS_CAMP_vs_HospAsthma_CRA", "CAMP", "CAMP_EDHOS_cum_Y1_binary"),
  predprob_one(cra_hosp$log, cra_hosp$pheno, "hsa-miR-584-5p", "outcome_binary",       "treatment", cra_main_covars,  "EDHOS_CAMP_vs_HospAsthma_CRA", "CRA",  "CRA_Hospitalized_Asthma_Last_Yr_binary"),
  predprob_one(camp$log,     camp$pheno,     "hsa-miR-205-5p", "EDHOS_cum_Y1_binary", "treatment", camp_main_covars, "EDHOS_CAMP_vs_ERVisits_CRA",   "CAMP", "CAMP_EDHOS_cum_Y1_binary"),
  predprob_one(cra_er$log,   cra_er$pheno,   "hsa-miR-205-5p", "outcome_binary",       "treatment", cra_main_covars,  "EDHOS_CAMP_vs_ERVisits_CRA",   "CRA",  "CRA_ER_Visits_Asthma_Last_Yr_binary")
)
write_csv(predprob_all, file.path(out_dir, "models", "candidate_adjusted_predprob_values.csv"))

# -----------------------------
# Predicted probability plots with interaction OR + p annotations
# -----------------------------
plot_df <- predprob_all |>
  mutate(
    cohort_outcome  = label_cohort_outcome(cohort, outcome),
    treatment_label = case_when(
      cohort == "CAMP" & as.character(treatment) == "1" ~ "ICS",
      cohort == "CAMP" & as.character(treatment) == "3" ~ "Non-ICS",
      cohort == "CRA"  & as.character(treatment) == "2" ~ "ICS",
      cohort == "CRA"  & as.character(treatment) == "1" ~ "Non-ICS",
      TRUE ~ as.character(treatment)
    ),
    treatment_label = factor(treatment_label, levels = c("Non-ICS", "ICS"))
  )

# Interaction OR + p pulled from model results, keyed to panel facet variable
predprob_annot <- main_results |>
  filter(model == "Interaction", !is.na(OR_CI)) |>
  mutate(
    cohort_outcome = label_cohort_outcome(cohort, outcome),
    OR_CI_label = OR_CI |>
      stringr::str_replace("\\(", "[") |>
      stringr::str_replace("\\)", "]"),
    annot_label = paste0("OR ", OR_CI_label, "\n", fmt_p(p_value))
  ) |>
  select(comparison, miRNA, cohort_outcome, annot_label)

draw_predprob <- function(mir_id, comparison_id, x_label) {
  pd <- plot_df |>
    filter(miRNA == mir_id, comparison == comparison_id)

  an <- predprob_annot |>
    filter(miRNA == mir_id, comparison == comparison_id)

  # Put annotation near top-right but inside each panel
  ann_pos <- pd |>
    group_by(cohort_outcome) |>
    summarise(
      x_pos = max(mir_expr, na.rm = TRUE),
      y_pos = max(upper, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(an, by = "cohort_outcome")

  ggplot(pd, aes(x = mir_expr, y = pred_prob,
                 color = treatment_label, fill = treatment_label)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.20, color = NA) +
    geom_line(linewidth = 1.35) +
    geom_text(
      data = ann_pos,
      aes(x = x_pos, y = y_pos, label = annot_label),
      inherit.aes = FALSE,
      hjust = 1.02, vjust = 1.05,
      size = 4.2,
      color = "gray15",
      lineheight = 0.95
    ) +
    facet_wrap(~ cohort_outcome, scales = "free_x") +
    scale_color_manual(values = c("Non-ICS" = "#E41A1C", "ICS" = "#377EB8")) +
    scale_fill_manual(values = c("Non-ICS" = "#E41A1C", "ICS" = "#377EB8")) +
    labs(
      x = x_label,
      y = "Adjusted predicted probability",
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "gray88", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray93", linewidth = 0.3),
      strip.background = element_rect(fill = "gray85", color = "gray35", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 13),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12, color = "gray25"),
      legend.position = "bottom",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.spacing = unit(1.0, "lines"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

p_mir205_er <- draw_predprob(
  "hsa-miR-205-5p", "EDHOS_CAMP_vs_ERVisits_CRA",
  "miR-205-5p expression (log2)")+coord_cartesian(ylim = c(0, 0.95))

p_mir584_hosp <- draw_predprob(
  "hsa-miR-584-5p", "EDHOS_CAMP_vs_HospAsthma_CRA",
  "miR-584-5p expression (log2)")+coord_cartesian(ylim = c(0, 0.40))

ggsave(file.path(out_dir, "figures", "mir205_edhos_er_predprob_adjusted.png"),              p_mir205_er,   width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "figures", "mir205_edhos_er_predprob_adjusted.pdf"),              p_mir205_er,   width = 10, height = 5)
ggsave(file.path(out_dir, "figures", "mir584_edhos_hospitalization_predprob_adjusted.png"), p_mir584_hosp, width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "figures", "mir584_edhos_hospitalization_predprob_adjusted.pdf"), p_mir584_hosp, width = 10, height = 5)

# -----------------------------
# CAMP_NEW independent replication (candidate miRNAs only)
# -----------------------------
camp_new_covars <- c("AGE", "SEX", "RACE")

camp_new_candidates <- intersect(candidate_mirs, rownames(camp_new$log))
missing_cands_new   <- setdiff(candidate_mirs, rownames(camp_new$log))
if (length(missing_cands_new) > 0)
  warning("CAMP_NEW: candidate miRNA(s) missing after filtering: ",
          paste(missing_cands_new, collapse = ", "))

camp_new_log_cands <- camp_new$log[camp_new_candidates, , drop = FALSE]

# Model A + B: binary interaction and stratified logistic
camp_new_binary <- run_binary_scan(
  camp_new_log_cands, camp_new$pheno,
  "EDHOS_cum_Y1_binary", "treatment", "treatment1:mir_expr", camp_new_covars,
  "CAMP_NEW_EDHOS_replication", "CAMP_NEW", "CAMP_NEW_EDHOS_cum_Y1_binary", "3", "1"
)
camp_new_for_meta <- camp_new_binary |> mutate(comparison = "CAMP_EDHOS")

camp_meta <- run_ivw_meta(
  bind_rows(
    camp_for_meta,
    camp_new_for_meta
  )
)
camp_meta_candidates <- camp_meta |> filter(miRNA %in% candidate_mirs)

write_csv(
  camp_meta_candidates,
  file.path(out_dir,
            "meta",
            "candidate_camp_old_new_ivw_meta_results.csv")
)
# Model C: ordinal sensitivity (treatment * mir_expr + AGE + SEX)
camp_new_sensitivity_raw <- map_dfr(camp_new_candidates, function(mir) {
  dat_ord <- prep_model_dat(camp_new$pheno, camp_new_log_cands, mir,
                            c("EDHOS_cum_Y1_ordered", "treatment", "AGE", "SEX"))
  ord_cv  <- clean_covars(dat_ord, c("AGE", "SEX"))
  bind_rows(
    safe_clm_coef(dat_ord, "EDHOS_cum_Y1_ordered",
                  rhs(c("treatment * mir_expr", ord_cv)), "treatment1:mir_expr") |>
      mutate(model = "Interaction", n = nrow(dat_ord)),
    map_dfr(c("3", "1"), function(g) {
      sub <- dat_ord |> filter(as.character(treatment) == g) |> droplevels()
      safe_clm_coef(sub, "EDHOS_cum_Y1_ordered",
                    rhs(c("mir_expr", clean_covars(sub, c("AGE", "SEX")))), "mir_expr") |>
        mutate(model = ifelse(g == "1", "Stratified: ICS", "Stratified: non-ICS"), n = nrow(sub))
    })
  ) |> mutate(miRNA = mir)
}) |>
  mutate(comparison = "CAMP_NEW_EDHOS_replication", cohort = "CAMP_NEW",
         outcome = "EDHOS_cum_Y1_ordered", analysis = "Ordinal sensitivity") |>
  group_by(comparison, cohort, outcome, analysis, model) |>
  mutate(FDR = p.adjust(p_value, "fdr")) |>
  ungroup()
camp_new_sensitivity <- bind_cols(camp_new_sensitivity_raw,
                                  or_ci(camp_new_sensitivity_raw$beta,
                                        camp_new_sensitivity_raw$se))

camp_new_results <- bind_rows(
  camp_new_binary      |> select(cohort, miRNA, analysis, model, n, beta, se, OR, lower_CI, upper_CI, OR_CI, p_value, FDR),
  camp_new_sensitivity |> select(cohort, miRNA, analysis, model, n, beta, se, OR, lower_CI, upper_CI, OR_CI, p_value, FDR)
)
write_csv(camp_new_results,
          file.path(out_dir, "models", "camp_new_candidate_replication_results.csv"))
cat("\nCAMP_NEW replication results:\n"); print(camp_new_results, n = Inf)

# -----------------------------
# Summary
# -----------------------------
cat("\nSaved outputs:\n")
invisible(lapply(c(
  file.path(out_dir, "models",  "candidate_main_binary_interaction_stratified_results.csv"),
  file.path(out_dir, "models",  "candidate_camp_edhos_ordinal_linear_results.csv"),
  file.path(out_dir, "models",  "candidate_auc_summary_all_ics_nonics.csv"),
  file.path(out_dir, "models",  "candidate_adjusted_predprob_values.csv"),
  file.path(out_dir, "models",  "all_common_mirna_binary_ivw_meta_results.csv"),
  file.path(out_dir, "models",  "candidate_binary_ivw_meta_results.csv"),
  file.path(out_dir, "figures", "camp_interaction_volcano.png"),
  file.path(out_dir, "figures", "camp_interaction_volcano.pdf"),
  file.path(out_dir, "figures", "candidate_interaction_forest_plot.png"),
  file.path(out_dir, "figures", "candidate_interaction_forest_plot.pdf"),
  file.path(out_dir, "figures", "candidate_ics_stratified_forest_plot.png"),
  file.path(out_dir, "figures", "candidate_ics_stratified_forest_plot.pdf"),
  file.path(out_dir, "figures", "mir205_edhos_er_predprob_adjusted.png"),
  file.path(out_dir, "figures", "mir205_edhos_er_predprob_adjusted.pdf"),
  file.path(out_dir, "figures", "mir584_edhos_hospitalization_predprob_adjusted.png"),
  file.path(out_dir, "figures", "mir584_edhos_hospitalization_predprob_adjusted.pdf"),
  file.path(out_dir, "models",  "camp_new_candidate_replication_results.csv")
), function(f) cat(" -", f, "\n")))

cat("\nMain binary results:\n");                  print(main_results, n = Inf)
cat("\nCAMP ordinal/linear candidate results:\n"); print(camp_sensitivity_candidates, n = Inf)
cat("\nAUC summary:\n");                           print(auc_summary, n = Inf)
cat("\nCandidate IVW meta results:\n");             print(meta_candidates, n = Inf)
