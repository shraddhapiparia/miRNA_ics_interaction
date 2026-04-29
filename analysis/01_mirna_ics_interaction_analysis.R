set.seed(2025)

library(edgeR)
library(preprocessCore)
library(DESeq2)
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pROC)
library(tibble)
library(purrr)
library(readr)
library(stringr)
library(forcats)
library(patchwork)
library(matrixStats)
library(scales)
library(meta)
library(ordinal)

# Paths (relative to repo root)
cra_counts_path    <- "data/CRA/2022-08-16_CRA_1159_exceRpt_miRNA_ReadCountscopy.csv"
cra_pheno_path     <- "data/CRA/PhenotypeInfo.csv"
cra_link_path      <- "data/CRA/cra_shortened_phenotype_files10-11.csv"
camp_counts_path   <- "data/CAMP/CAMP.492.COUNT.csv"
camp_pheno_path    <- "data/CAMP/camp_pheno_2022-08-04.tsv"
camp_link_path     <- "data/CAMP/shortened_phenotype_files_ver5.csv"
if (!exists("out_dir")) out_dir <- "results"
highlight_mir      <- "hsa-miR-584-5p"
min_count          <- 5
min_prop_samples   <- 0.5

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "demographics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Helper functions

filter_mirnas <- function(count_df, min_count = 5, min_prop_samples = 0.5) {
  mat <- as.matrix(count_df[, -1, drop = FALSE])
  rownames(mat) <- count_df[[1]]
  class(mat) <- "numeric"
  keep <- rowSums(mat > min_count) >= (ncol(mat) * min_prop_samples)
  mat[keep, , drop = FALSE]
}

norm_raw_log <- function(mat) log2(mat + 1)

norm_tmm_log <- function(mat) {
  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge, method = "TMM")
  log2(cpm(dge, log = FALSE, prior.count = 0.5) + 1)
}

norm_deseq2_log <- function(mat) {
  coldata <- data.frame(row.names = colnames(mat), group = factor(rep("all", ncol(mat))))
  dds <- DESeqDataSetFromMatrix(countData = round(mat), colData = coldata, design = ~1)
  dds <- estimateSizeFactors(dds)
  normalized <- counts(dds, normalized = TRUE)
  log2(normalized + 1)
}

norm_quantile_log <- function(mat) {
  qn <- normalize.quantiles(as.matrix(mat))
  rownames(qn) <- rownames(mat)
  colnames(qn) <- colnames(mat)
  log2(qn + 1)
}

get_rle_df <- function(log_mat, method_name, cohort_name) {
  med <- rowMedians(log_mat, na.rm = TRUE)
  rle_mat <- sweep(log_mat, 1, med, FUN = "-")
  as.data.frame(rle_mat) |>
    rownames_to_column("miRNA") |>
    pivot_longer(-miRNA, names_to = "sample_id", values_to = "rle") |>
    mutate(method = method_name, cohort = cohort_name)
}

get_mean_sd_df <- function(log_mat, method_name, cohort_name) {
  tibble(
    miRNA = rownames(log_mat),
    mean_expr = rowMeans(log_mat, na.rm = TRUE),
    sd_expr = matrixStats::rowSds(log_mat, na.rm = TRUE),
    method = method_name,
    cohort = cohort_name
  )
}

get_density_df <- function(log_mat, method_name, cohort_name) {
  as.data.frame(log_mat) |>
    rownames_to_column("miRNA") |>
    pivot_longer(-miRNA, names_to = "sample_id", values_to = "expression") |>
    mutate(method = method_name, cohort = cohort_name)
}

run_pca_df <- function(log_mat, method_name, cohort_name, pheno_df, treatment_col) {
  pca <- prcomp(t(log_mat), center = TRUE, scale. = TRUE)
  out <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
    rownames_to_column("sample_id") |>
    left_join(pheno_df |> rownames_to_column("sample_id"), by = "sample_id") |>
    mutate(method = method_name, cohort = cohort_name, treatment = .data[[treatment_col]])
  attr(out, "variance_explained") <- summary(pca)$importance[2, 1:2]
  out
}

plot_density <- function(df, title_text) {
  ggplot(df, aes(x = expression, group = sample_id, color = sample_id)) +
    geom_density(alpha = 0.35) +
    facet_wrap(~ method, scales = "free") +
    labs(title = title_text, x = "log2 expression", y = "Density") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
}

plot_rle <- function(df, title_text) {
  ggplot(df, aes(x = sample_id, y = rle)) +
    geom_boxplot(outlier.size = 0.15) +
    facet_wrap(~ method, scales = "free_x") +
    labs(title = title_text, x = "Sample", y = "Relative log expression") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

plot_mean_sd <- function(df, title_text) {
  ggplot(df, aes(x = mean_expr, y = sd_expr)) +
    geom_point(alpha = 0.45, size = 0.8) +
    geom_smooth(se = FALSE, linewidth = 0.8) +
    facet_wrap(~ method, scales = "free") +
    labs(title = title_text, x = "Mean log2 expression", y = "SD") +
    theme_bw(base_size = 12)
}

plot_pca <- function(df, title_text) {
  ggplot(df, aes(x = PC1, y = PC2, color = treatment)) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~ method, scales = "free") +
    labs(title = title_text, x = "PC1", y = "PC2", color = "Treatment") +
    theme_bw(base_size = 12)
}

make_corr_plot <- function(mat1, mat2, label1, label2, cohort_name) {
  common_rows <- intersect(rownames(mat1), rownames(mat2))
  common_cols <- intersect(colnames(mat1), colnames(mat2))
  x <- as.vector(mat1[common_rows, common_cols, drop = FALSE])
  y <- as.vector(mat2[common_rows, common_cols, drop = FALSE])
  corr <- cor(x, y, use = "complete.obs", method = "pearson")

  df <- tibble(x = x, y = y)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.12, size = 0.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    annotate("text", x = quantile(x, 0.1), y = quantile(y, 0.95),
             label = paste0("r = ", round(corr, 3)), hjust = 0, size = 4) +
    labs(
      title = paste0(cohort_name, ": ", label1, " vs ", label2),
      x = label1,
      y = label2
    ) +
    theme_bw(base_size = 12)
}

summarize_continuous_by_group <- function(df, vars, group_col) {
  map_dfr(vars, function(v) {
    sub <- df |> dplyr::select(all_of(group_col), all_of(v)) |> filter(!is.na(.data[[v]]))
    pval <- tryCatch(t.test(sub[[v]] ~ sub[[group_col]])$p.value, error = function(e) NA_real_)
    sub |>
      group_by(.data[[group_col]]) |>
      summarise(n = n(), mean = mean(.data[[v]], na.rm = TRUE), sd = sd(.data[[v]], na.rm = TRUE), .groups = "drop") |>
      mutate(variable = v, p_value = pval) |>
      relocate(variable)
  })
}

summarize_categorical_by_group <- function(df, vars, group_col) {
  map_dfr(vars, function(v) {
    tab <- table(df[[group_col]], df[[v]], useNA = "no")
    pval <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
    as.data.frame(tab) |>
      rename(group = Var1, level = Var2, n = Freq) |>
      mutate(variable = v, p_value = pval) |>
      relocate(variable)
  })
}

plot_demo_violin <- function(df, vars, labels, group_col, title_text) {
  long_df <- df |>
    dplyr::select(all_of(group_col), all_of(vars)) |>
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value")

  ggplot(long_df, aes(x = .data[[group_col]], y = value, fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = labels)) +
    labs(title = title_text, x = NULL, y = "Value") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
}

extract_glm_term <- function(model, term_name) {
  coef_tab <- summary(model)$coefficients
  if (!term_name %in% rownames(coef_tab)) {
    return(tibble(term = term_name, estimate = NA_real_, se = NA_real_, p_value = NA_real_, or = NA_real_))
  }
  tibble(
    term = term_name,
    estimate = coef_tab[term_name, "Estimate"],
    se = coef_tab[term_name, "Std. Error"],
    p_value = coef_tab[term_name, "Pr(>|z|)"],
    or = exp(coef_tab[term_name, "Estimate"])
  )
}

run_interaction_scan <- function(
  log_mat,
  pheno_df,
  outcome_col,
  treatment_col,
  treatment_term,
  covars = NULL
) {
  mirs <- rownames(log_mat)

  map_dfr(mirs, function(mir) {
    dat <- pheno_df
    dat$mir_expr <- as.numeric(log_mat[mir, rownames(pheno_df)])

    needed_cols <- unique(c(outcome_col, treatment_col, "mir_expr", covars))
    dat <- dat |>
      dplyr::select(all_of(needed_cols)) |>
      tidyr::drop_na()

    full_rhs <- paste(c(
      paste0(treatment_col, " * mir_expr"),
      covars
    ), collapse = " + ")

    fit_all <- glm(
      as.formula(paste(outcome_col, "~", full_rhs)),
      data = dat,
      family = binomial()
    )

    treated_level <- levels(dat[[treatment_col]])[2]
    untreated_level <- levels(dat[[treatment_col]])[1]

    strat_rhs <- paste(c("mir_expr", covars), collapse = " + ")

    fit_treated <- glm(
      as.formula(paste(outcome_col, "~", strat_rhs)),
      data = dat |> filter(.data[[treatment_col]] == treated_level),
      family = binomial()
    )

    fit_untreated <- glm(
      as.formula(paste(outcome_col, "~", strat_rhs)),
      data = dat |> filter(.data[[treatment_col]] == untreated_level),
      family = binomial()
    )

    inter <- extract_glm_term(fit_all, treatment_term)
    treated <- extract_glm_term(fit_treated, "mir_expr")
    untreated <- extract_glm_term(fit_untreated, "mir_expr")

    tibble(
      miRNA = mir,
      n_all = nobs(fit_all),

      beta_interaction = log(inter$or),
      OR_interaction = inter$or,
      p_val_interaction = inter$p_value,
      se_interaction = inter$se,

      beta_treated = log(treated$or),
      OR_treated = treated$or,
      p_val_treated = treated$p_value,
      se_treated = treated$se,

      beta_untreated = log(untreated$or),
      OR_untreated = untreated$or,
      p_val_untreated = untreated$p_value,
      se_untreated = untreated$se
    )
  }) |>
    mutate(
      padj_interaction = p.adjust(p_val_interaction, method = "fdr"),
      padj_treated = p.adjust(p_val_treated, method = "fdr"),
      padj_untreated = p.adjust(p_val_untreated, method = "fdr")
    )
}

make_volcano <- function(results_df, title_text, highlight_mir) {
  plot_df <- results_df |>
    mutate(
      neglog10p = -log10(p_val_interaction),
      category = case_when(
        miRNA == highlight_mir ~ "Highlighted",
        p_val_interaction < 0.05 ~ "Significant",
        TRUE ~ "Non-significant"
      )
    )

  ggplot(plot_df, aes(x = OR_interaction, y = neglog10p)) +
    geom_point(aes(color = category, shape = category), alpha = 0.85, size = 2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = subset(plot_df, category != "Non-significant"),
      aes(label = miRNA),
      size = 3,
      max.overlaps = 30
    ) +
    scale_color_manual(values = c("Highlighted" = "blue", "Significant" = "firebrick", "Non-significant" = "gray70")) +
    scale_shape_manual(values = c("Highlighted" = 17, "Significant" = 16, "Non-significant" = 16)) +
    labs(title = title_text, x = "Odds ratio (interaction)", y = expression(-log[10](p))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

make_predprob_df <- function(log_mat, pheno_df, mir_name, outcome_col, treatment_col, treated_level) {
  dat <- pheno_df
  dat$mir_expr <- as.numeric(log_mat[mir_name, rownames(pheno_df)])
  fit <- glm(as.formula(paste(outcome_col, "~", treatment_col, "* mir_expr")), data = dat, family = binomial())

  mir_grid <- seq(quantile(dat$mir_expr, 0.1, na.rm = TRUE), quantile(dat$mir_expr, 0.9, na.rm = TRUE), length.out = 100)
  newdat <- expand.grid(mir_expr = mir_grid, tmp_group = levels(dat[[treatment_col]]))
  names(newdat)[2] <- treatment_col
  pred <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
  newdat$fit <- pred$fit
  newdat$se <- pred$se.fit
  newdat$pred_prob <- plogis(newdat$fit)
  newdat$lower <- plogis(newdat$fit - 1.96 * newdat$se)
  newdat$upper <- plogis(newdat$fit + 1.96 * newdat$se)
  list(model = fit, pred_df = newdat)
}

plot_predprob <- function(pred_df, treatment_col, title_text) {
  ggplot(pred_df, aes(x = mir_expr, y = pred_prob, color = .data[[treatment_col]], fill = .data[[treatment_col]])) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, linewidth = 0) +
    geom_line(linewidth = 1) +
    labs(title = title_text, x = "miRNA expression", y = "Predicted probability", color = "Treatment", fill = "Treatment") +
    theme_bw(base_size = 12)
}

make_predprob_df_adjusted <- function(log_mat, pheno_df, mir_name, outcome_col,
                                       treatment_col, treated_level, covars) {
  dat <- pheno_df
  dat$mir_expr <- as.numeric(log_mat[mir_name, rownames(pheno_df)])

  needed_cols <- unique(c(outcome_col, treatment_col, "mir_expr", covars))
  dat <- dat |>
    dplyr::select(all_of(needed_cols)) |>
    tidyr::drop_na()

  rhs <- paste(c(paste0(treatment_col, " * mir_expr"), covars), collapse = " + ")
  fit <- glm(as.formula(paste(outcome_col, "~", rhs)), data = dat, family = binomial())

  mir_grid <- seq(quantile(dat$mir_expr, 0.1, na.rm = TRUE),
                  quantile(dat$mir_expr, 0.9, na.rm = TRUE), length.out = 100)
  newdat <- expand.grid(mir_expr = mir_grid, tmp_group = levels(dat[[treatment_col]]))
  names(newdat)[2] <- treatment_col

  for (cv in covars) {
    v <- dat[[cv]]
    if (is.numeric(v)) {
      newdat[[cv]] <- mean(v, na.rm = TRUE)
    } else {
      newdat[[cv]] <- factor(levels(v)[1], levels = levels(v))
    }
  }

  pred <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
  newdat$fit       <- pred$fit
  newdat$se        <- pred$se.fit
  newdat$pred_prob <- plogis(newdat$fit)
  newdat$lower     <- plogis(newdat$fit - 1.96 * newdat$se)
  newdat$upper     <- plogis(newdat$fit + 1.96 * newdat$se)

  list(model = fit, pred_df_adjusted = newdat)
}

plot_predprob_adjusted <- function(pred_df, treatment_col, title_text) {
  ggplot(pred_df, aes(x = mir_expr, y = pred_prob, color = .data[[treatment_col]], fill = .data[[treatment_col]])) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, linewidth = 0) +
    geom_line(linewidth = 1) +
    labs(title = title_text, x = "miRNA expression", y = "Adjusted predicted probability",
         color = "Treatment", fill = "Treatment") +
    theme_bw(base_size = 12)
}

compute_auc_summary <- function(log_mat, pheno_df, mir_name, outcome_col, treatment_col, covars = NULL) {
  dat <- pheno_df
  dat$Expression <- as.numeric(log_mat[mir_name, rownames(pheno_df)])

  if (!is.null(covars)) {
    for (cv in covars) {
      if (cv %in% names(dat) && anyNA(dat[[cv]])) {
        dat[[cv]][is.na(dat[[cv]])] <- mean(dat[[cv]], na.rm = TRUE)
      }
    }
  }

  rhs_full <- paste(c(paste0(treatment_col, " * Expression"), covars), collapse = " + ")
  model_full <- glm(as.formula(paste(outcome_col, "~", rhs_full)), data = dat, family = binomial())
  pred_full <- predict(model_full, type = "response")
  roc_full <- roc(dat[[outcome_col]], pred_full, quiet = TRUE)

  stratified <- map_dfr(levels(dat[[treatment_col]]), function(g) {
    sub <- dat |> filter(.data[[treatment_col]] == g)
    rhs_sub <- paste(c("Expression", covars), collapse = " + ")
    fit <- glm(as.formula(paste(outcome_col, "~", rhs_sub)), data = sub, family = binomial())
    pred <- predict(fit, type = "response")
    roc_obj <- roc(sub[[outcome_col]], pred, quiet = TRUE)
    tibble(group = as.character(g), auc = as.numeric(auc(roc_obj)))
  })

  list(
    interaction_auc = as.numeric(auc(roc_full)),
    stratified_auc = stratified,
    model = model_full,
    roc_full = roc_full
  )
}

make_volcano_pub <- function(results_df, title_text, highlight_mir) {
  plot_df <- results_df |>
    mutate(
      neglog10p = -log10(p_val_interaction),
      category  = case_when(
        miRNA == highlight_mir   ~ "Highlighted",
        p_val_interaction < 0.05 ~ "Significant",
        TRUE                     ~ "Non-significant"
      )
    )
  x_range  <- range(plot_df$OR_interaction, na.rm = TRUE)
  y_max    <- max(plot_df$neglog10p, na.rm = TRUE)
  sig_line <- -log10(0.05)

  ggplot(plot_df, aes(x = OR_interaction, y = neglog10p)) +
    geom_point(aes(color = category, shape = category), alpha = 0.8, size = 1.8) +
    geom_vline(xintercept = 1,        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    geom_hline(yintercept = sig_line, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    annotate("text", x = x_range[2], y = y_max,
             label = "Higher risk\nin ICS", hjust = 1, vjust = 1,
             size = 2.8, color = "gray30", fontface = "italic") +
    annotate("text", x = x_range[1], y = y_max,
             label = "Lower risk\nin ICS", hjust = 0, vjust = 1,
             size = 2.8, color = "gray30", fontface = "italic") +
    annotate("text", x = 1, y = sig_line * 0.45,
             label = "lower evidence", hjust = 0.5,
             size = 2.5, color = "gray50") +
    geom_text_repel(
      data         = subset(plot_df, category != "Non-significant"),
      aes(label    = miRNA),
      size         = 2.8, max.overlaps = 30,
      segment.size = 0.3, box.padding  = 0.4
    ) +
    scale_color_manual(
      values = c("Highlighted" = "#2166ac", "Significant" = "#d6604d",
                 "Non-significant" = "gray75"),
      name = NULL
    ) +
    scale_shape_manual(
      values = c("Highlighted" = 17, "Significant" = 16, "Non-significant" = 16),
      name = NULL
    ) +
    labs(title = title_text,
         x = "Odds ratio (interaction)",
         y = expression(-log[10](italic(p)))) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "top",
          panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 11, face = "bold"))
}

plot_predprob_combined <- function(camp_pred_df, cra_pred_df,
                                    camp_ics_level, cra_ics_level) {
  pal <- c("ICS" = "#d6604d", "Non-ICS" = "#4393c3")

  add_label <- function(df, ics_level) {
    df |> mutate(
      group_label = factor(
        ifelse(as.character(treatment) == as.character(ics_level), "ICS", "Non-ICS"),
        levels = c("ICS", "Non-ICS")
      )
    )
  }

  cp <- add_label(camp_pred_df, camp_ics_level)
  cr <- add_label(cra_pred_df,  cra_ics_level)

  make_panel <- function(df, tag_label, y_lab) {
    ggplot(df, aes(x = mir_expr, y = pred_prob,
                   color = group_label, fill = group_label)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = pal, name = NULL) +
      scale_fill_manual(values  = pal, name = NULL) +
      labs(tag = tag_label, x = "miR-584-5p expression", y = y_lab) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
  }

  p_a <- make_panel(cp, "A", "Predicted probability\nof exacerbation (CAMP)")
  p_b <- make_panel(cr, "B", "Predicted probability\nof hospitalization (CRA)")

  (p_a | p_b) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

make_demographics_table <- function(cohort_name, pheno_df, cont_vars, cat_vars,
                                     group_col, ics_level) {
  pheno_df$group_label <- ifelse(
    as.character(pheno_df[[group_col]]) == as.character(ics_level), "ICS", "Non-ICS"
  )
  label_map <- pheno_df |>
    dplyr::select(group = all_of(group_col), group_label) |>
    distinct() |>
    mutate(group = as.character(group))

  cont_rows <- map_dfr(cont_vars, function(v) {
    sub <- pheno_df |>
      dplyr::select(group = all_of(group_col), group_label, val = all_of(v)) |>
      filter(!is.na(val))
    pval <- tryCatch(t.test(sub$val ~ sub$group)$p.value, error = function(e) NA_real_)
    sub |>
      group_by(group, group_label) |>
      summarise(n = n(), mean = mean(val), sd = sd(val), .groups = "drop") |>
      mutate(cohort = cohort_name, variable = v, variable_type = "continuous",
             level = NA_character_, percent = NA_real_, p_value = pval,
             group = as.character(group))
  })

  cat_rows <- map_dfr(cat_vars, function(v) {
    tab  <- table(pheno_df[[group_col]], pheno_df[[v]], useNA = "no")
    pval <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
    as.data.frame(tab) |>
      rename(group = Var1, level = Var2, n = Freq) |>
      mutate(group = as.character(group)) |>
      group_by(group) |>
      mutate(percent = 100 * n / sum(n)) |>
      ungroup() |>
      left_join(label_map, by = "group") |>
      mutate(cohort = cohort_name, variable = v, variable_type = "categorical",
             mean = NA_real_, sd = NA_real_, p_value = pval)
  })

  bind_rows(cont_rows, cat_rows) |>
    dplyr::select(cohort, variable, variable_type, group, group_label,
                  level, n, percent, mean, sd, p_value)
}

detect_sensitivity_outcomes <- function(pheno_df, exclude_cols) {
  candidates  <- setdiff(names(pheno_df), exclude_cols)
  cont_col    <- NULL
  ordinal_col <- NULL
  rows        <- list()

  for (col in candidates) {
    v      <- pheno_df[[col]]
    pct_na <- mean(is.na(v))
    if (pct_na >= 0.5) next

    if (is.numeric(v)) {
      n_uniq <- length(unique(na.omit(v)))
      sdv    <- sd(v, na.rm = TRUE)
      if (!is.na(sdv) && sdv > 0 && n_uniq >= 5) {
        if (is.null(cont_col)) cont_col <- col
        rows <- c(rows, list(data.frame(
          analysis_type    = "continuous",
          candidate_column = col,
          selected_column  = cont_col,
          status           = if (col == cont_col) "selected" else "candidate",
          reason           = sprintf("sd=%.3f, n_unique=%d", sdv, n_uniq),
          stringsAsFactors = FALSE
        )))
      }
    }

    if (is.factor(v) && nlevels(v) >= 3) {
      if (is.null(ordinal_col)) ordinal_col <- col
      rows <- c(rows, list(data.frame(
        analysis_type    = "ordinal",
        candidate_column = col,
        selected_column  = ordinal_col,
        status           = if (col == ordinal_col) "selected" else "candidate",
        reason           = sprintf("levels=%d", nlevels(v)),
        stringsAsFactors = FALSE
      )))
    }
  }

  list(cont = cont_col, ordinal = ordinal_col,
       check = if (length(rows)) bind_rows(rows) else data.frame())
}

# CAMP cohort

camp_counts <- read.csv(camp_counts_path, check.names = FALSE)
camp_pheno <- read.table(camp_pheno_path, sep = "\t", header = TRUE, check.names = FALSE)
camp_link <- read.csv(camp_link_path, check.names = FALSE)

camp_pheno <- camp_pheno[camp_pheno$TG != 2, , drop = FALSE]
camp_mat <- filter_mirnas(camp_counts, min_count, min_prop_samples)

camp_pheno_merged <- merge(camp_pheno, camp_link, by.x = "camp", by.y = "CAMP_ID", all=FALSE)
rownames(camp_pheno_merged) <- camp_pheno_merged$`row.names(CAMP_492_COUNT.Trans)`

camp_pheno_merged <- camp_pheno_merged |>
  filter(!is.na(EDHOS_cum_Y1)) |>
  mutate(
    outcome = factor(ifelse(EDHOS_cum_Y1 < 1, "low", "high")),
    outcome = relevel(outcome, ref = "low"),
    treatment = factor(TG),
    treatment = relevel(treatment, ref = "3"),
    SEX.x = factor(SEX.x)
  )

camp_pheno_merged$RACE.x = factor(camp_pheno_merged$RACE.x)
common_camp <- intersect(colnames(camp_mat), rownames(camp_pheno_merged))
camp_mat <- camp_mat[, common_camp, drop = FALSE]
camp_pheno_merged <- camp_pheno_merged[common_camp, , drop = FALSE]

# CAMP demographics
camp_cont_vars <- c("AGE.x", "PREFEVPP_RZ", "BMI_RZ")
camp_cat_vars <- c("SEX.x", "outcome")

camp_demo <- make_demographics_table(
  "CAMP", camp_pheno_merged, c(camp_cont_vars, "EDHOS_cum_Y1"), camp_cat_vars, "treatment", ics_level = "1"
)
write_csv(camp_demo, file.path(out_dir, "demographics", "CAMP_demographics_summary.csv"))

camp_demo_plot <- plot_demo_violin(
  camp_pheno_merged,
  vars = camp_cont_vars,
  labels = c(AGE.x = "Age", PREFEVPP_RZ = "FEV1 (% predicted)", BMI_RZ = "BMI"),
  group_col = "treatment",
  title_text = "CAMP demographics by treatment group"
)

# CAMP normalization comparison
camp_norms <- list(
  Raw = norm_raw_log(camp_mat),
  TMM = norm_tmm_log(camp_mat),
  DESeq2 = norm_deseq2_log(camp_mat),
  Quantile = norm_quantile_log(camp_mat)
)

camp_density_df <- imap_dfr(camp_norms, ~ get_density_df(.x, .y, "CAMP"))
camp_rle_df <- imap_dfr(camp_norms, ~ get_rle_df(.x, .y, "CAMP"))
camp_meansd_df <- imap_dfr(camp_norms, ~ get_mean_sd_df(.x, .y, "CAMP"))
camp_pca_df <- imap_dfr(camp_norms, ~ run_pca_df(.x, .y, "CAMP", camp_pheno_merged, "treatment"))

p_camp_density <- plot_density(camp_density_df, "CAMP: expression distributions across normalization methods")
p_camp_rle <- plot_rle(camp_rle_df, "CAMP: RLE across normalization methods")
p_camp_meansd <- plot_mean_sd(camp_meansd_df, "CAMP: mean-SD trend across normalization methods")
p_camp_pca <- plot_pca(camp_pca_df, "CAMP: PCA across normalization methods")

ggsave(file.path(out_dir, "qc", "camp_density.png"), p_camp_density, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "qc", "camp_rle.png"), p_camp_rle, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "qc", "camp_mean_sd.png"), p_camp_meansd, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "qc", "camp_pca.png"), p_camp_pca, width = 12, height = 8, dpi = 300)

p_camp_corr <- make_corr_plot(camp_norms$DESeq2, camp_norms$Quantile, "DESeq2 log-normalized", "Quantile log-normalized", "CAMP")
ggsave(file.path(out_dir, "qc", "camp_deseq2_vs_quantile_corr.png"), p_camp_corr, width = 6, height = 5, dpi = 300)

# CAMP interaction analysis
log_quantile_counts_camp <- camp_norms$Quantile

camp_results <- run_interaction_scan(
  log_mat = log_quantile_counts_camp,
  pheno_df = camp_pheno_merged,
  outcome_col = "outcome",
  treatment_col = "treatment",
  treatment_term = "treatment1:mir_expr",
  covars = c("AGE.x", "SEX.x", "RACE.x")
)

write_csv(camp_results, file.path(out_dir, "models", "camp_interaction_results.csv"))

p_camp_volcano <- make_volcano(camp_results, "CAMP: miRNA × ICS interaction effects", highlight_mir)
ggsave(file.path(out_dir, "figures", "camp_interaction_volcano.png"), p_camp_volcano, width = 7, height = 6, dpi = 300)

p_camp_volcano_pub <- make_volcano_pub(camp_results, "CAMP: miRNA × ICS interaction effects", highlight_mir)
ggsave(file.path(out_dir, "figures", "camp_interaction_volcano_pub.png"), p_camp_volcano_pub, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "figures", "camp_interaction_volcano_pub.pdf"), p_camp_volcano_pub, width = 7, height = 6)

# CAMP miR-584-5p follow-up
camp_pred <- make_predprob_df(
  log_mat = log_quantile_counts_camp,
  pheno_df = camp_pheno_merged,
  mir_name = highlight_mir,
  outcome_col = "outcome",
  treatment_col = "treatment",
  treated_level = "1"
)

p_camp_pred <- plot_predprob(camp_pred$pred_df, "treatment", "CAMP: predicted exacerbation risk by miR-584-5p and ICS")
ggsave(file.path(out_dir, "figures", "camp_mir584_predprob.png"), p_camp_pred, width = 7, height = 5, dpi = 300)

camp_pred_adj <- make_predprob_df_adjusted(
  log_mat       = log_quantile_counts_camp,
  pheno_df      = camp_pheno_merged,
  mir_name      = highlight_mir,
  outcome_col   = "outcome",
  treatment_col = "treatment",
  treated_level = "1",
  covars        = c("AGE.x", "SEX.x", "RACE.x")
)
p_camp_pred_adj <- plot_predprob_adjusted(
  camp_pred_adj$pred_df_adjusted, "treatment",
  "CAMP: adjusted predicted exacerbation risk by miR-584-5p and ICS"
)
ggsave(file.path(out_dir, "figures", "camp_mir584_predprob_adjusted.png"),
       p_camp_pred_adj, width = 7, height = 5, dpi = 300)

camp_auc <- compute_auc_summary(
  log_mat = log_quantile_counts_camp,
  pheno_df = camp_pheno_merged,
  mir_name = highlight_mir,
  outcome_col = "outcome",
  treatment_col = "treatment",
  covars = c("PREFEVPP_RZ", "AGEON", "BMI_RZ", "AGE.x", "TOT.EOS_S3", "LOG10IGE_S3")
)

write_csv(camp_auc$stratified_auc, file.path(out_dir, "models", "camp_mir584_auc_stratified.csv"))
print(camp_auc$interaction_auc)
print(camp_auc$stratified_auc)

# CAMP sensitivity analysis — explicit outcome models for miR-584-5p
camp_sens_dat <- camp_pheno_merged
camp_sens_dat$mir_expr <- as.numeric(
  log_quantile_counts_camp[highlight_mir, rownames(camp_pheno_merged)]
)
camp_sens_dat <- camp_sens_dat[
  !is.na(camp_sens_dat$EDHOS_cum_Y1) & !is.na(camp_sens_dat$mir_expr), ]

fit_camp_lm <- tryCatch(
  lm(EDHOS_cum_Y1 ~ treatment * mir_expr + AGE.x + SEX.x + RACE.x,
     data = camp_sens_dat),
  error = function(e) { message("CAMP lm sensitivity: ", e$message); NULL }
)
if (!is.null(fit_camp_lm)) {
  write_csv(
    broom::tidy(fit_camp_lm) |> mutate(outcome = "EDHOS_cum_Y1"),
    file.path(out_dir, "models", "camp_linear_interaction_results.csv")
  )
}

camp_ord_n_levels <- length(unique(na.omit(camp_sens_dat$EDHOS_cum_Y1)))
fit_camp_clm <- NULL

if (camp_ord_n_levels > 2) {
  camp_sens_dat$outcome_ord <- factor(camp_sens_dat$EDHOS_cum_Y1, ordered = TRUE)

  fit_camp_clm <- tryCatch(
    ordinal::clm(
      outcome_ord ~ treatment * mir_expr + AGE.x + SEX.x + RACE.x,
      data = camp_sens_dat,
      link = "logit"
    ),
    error = function(e) { message("CAMP clm sensitivity: ", e$message); NULL }
  )

  if (!is.null(fit_camp_clm)) {
    coef_table <- coef(summary(fit_camp_clm))

    clm_out <- as.data.frame(coef_table) |>
      tibble::rownames_to_column("term") |>
      dplyr::filter(!grepl("\\|", term)) |>
      dplyr::rename(
        estimate = Estimate,
        std.error = `Std. Error`,
        statistic = `z value`,
        p_value = `Pr(>|z|)`
      ) |>
      dplyr::mutate(
        OR = exp(estimate),
        outcome = "EDHOS_cum_Y1_ordered"
      )

    write_csv(
      clm_out,
      file.path(out_dir, "models", "camp_ordinal_interaction_results.csv")
    )
  }
}

camp_check_df <- data.frame(
  cohort = "CAMP",
  analysis_type = c("linear", "ordinal"),
  selected_column = c("EDHOS_cum_Y1",
                      if (camp_ord_n_levels > 2) "EDHOS_cum_Y1" else NA_character_),
  status = c(if (!is.null(fit_camp_lm)) "fitted" else "error",
             if (!is.null(fit_camp_clm)) "fitted" else "no_predefined_ordinal_outcome"),
  reason = c("predefined: EDHOS_cum_Y1",
             if (camp_ord_n_levels > 2) "predefined: EDHOS_cum_Y1 as ordered factor"
             else "no predefined ordinal exacerbation column available"),
  stringsAsFactors = FALSE
)

# CRA cohort

cra_counts <- read.csv(cra_counts_path, check.names = FALSE)
cra_pheno <- read.csv(cra_pheno_path, check.names = FALSE)
cra_link <- read.csv(cra_link_path, check.names = FALSE)

cra_mat <- filter_mirnas(cra_counts, min_count, min_prop_samples)

cra_pheno_merged <- merge(cra_pheno, cra_link, by.x = "S_SUBJECTID", by.y = "ST.ID")
cra_pheno_merged$sample.id <- gsub("\\-", ".", cra_pheno_merged$sample.id)
rownames(cra_pheno_merged) <- cra_pheno_merged$sample.id

cra_pheno_merged <- cra_pheno_merged |>
  filter(!is.na(Hospitalized_Asthma_Last_Yr)) |>
  mutate(
    outcome = factor(Hospitalized_Asthma_Last_Yr),
    treatment = factor(Inhaled_Steroids),
    treatment = relevel(treatment, ref = "1"),
    gender.x = factor(gender.x)
  )

common_cra <- intersect(colnames(cra_mat), rownames(cra_pheno_merged))
cra_mat <- cra_mat[, common_cra, drop = FALSE]
cra_pheno_merged <- cra_pheno_merged[common_cra, , drop = FALSE]

# CRA demographics
cra_cont_vars <- c("age.x", "pctpred_fev1_pre_BD", "bmi_pct")
cra_cat_vars <- c("gender.x", "Hospitalized_Asthma_Last_Yr")

cra_demo <- make_demographics_table(
  "CRA", cra_pheno_merged, c(cra_cont_vars, "Hospitalized_Asthma_Last_Yr"), cra_cat_vars, "treatment", ics_level = "2"
)
write_csv(cra_demo, file.path(out_dir, "demographics", "CRA_demographics_summary.csv"))

cra_demo_plot <- plot_demo_violin(
  cra_pheno_merged,
  vars = cra_cont_vars,
  labels = c(age.x = "Age", pctpred_fev1_pre_BD = "FEV1 (% predicted)", bmi_pct = "BMI percentile"),
  group_col = "treatment",
  title_text = "CRA demographics by treatment group"
)

# CRA normalization comparison
cra_norms <- list(
  Raw = norm_raw_log(cra_mat),
  TMM = norm_tmm_log(cra_mat),
  DESeq2 = norm_deseq2_log(cra_mat),
  Quantile = norm_quantile_log(cra_mat)
)

cra_density_df <- imap_dfr(cra_norms, ~ get_density_df(.x, .y, "CRA"))
cra_rle_df <- imap_dfr(cra_norms, ~ get_rle_df(.x, .y, "CRA"))
cra_meansd_df <- imap_dfr(cra_norms, ~ get_mean_sd_df(.x, .y, "CRA"))
cra_pca_df <- imap_dfr(cra_norms, ~ run_pca_df(.x, .y, "CRA", cra_pheno_merged, "treatment"))

p_cra_density <- plot_density(cra_density_df, "CRA: expression distributions across normalization methods")
p_cra_rle <- plot_rle(cra_rle_df, "CRA: RLE across normalization methods")
p_cra_meansd <- plot_mean_sd(cra_meansd_df, "CRA: mean-SD trend across normalization methods")
p_cra_pca <- plot_pca(cra_pca_df, "CRA: PCA across normalization methods")

ggsave(file.path(out_dir, "qc", "cra_density.png"), p_cra_density, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "qc", "cra_rle.png"), p_cra_rle, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "qc", "cra_mean_sd.png"), p_cra_meansd, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "qc", "cra_pca.png"), p_cra_pca, width = 12, height = 8, dpi = 300)

p_cra_corr <- make_corr_plot(cra_norms$DESeq2, cra_norms$Quantile, "DESeq2 log-normalized", "Quantile log-normalized", "CRA")
ggsave(file.path(out_dir, "qc", "cra_deseq2_vs_quantile_corr.png"), p_cra_corr, width = 6, height = 5, dpi = 300)

# CRA interaction analysis
log_quantile_counts_cra <- cra_norms$Quantile

cra_results <- run_interaction_scan(
  log_mat = log_quantile_counts_cra,
  pheno_df = cra_pheno_merged,
  outcome_col = "outcome",
  treatment_col = "treatment",
  treatment_term = "treatment2:mir_expr",
  covars = c("age.x", "gender.x")
)

write_csv(cra_results, file.path(out_dir, "models", "cra_interaction_results.csv"))

p_cra_volcano <- make_volcano(cra_results, "CRA: miRNA × ICS interaction effects", highlight_mir)
ggsave(file.path(out_dir, "figures", "cra_interaction_volcano.png"), p_cra_volcano, width = 7, height = 6, dpi = 300)

p_cra_volcano_pub <- make_volcano_pub(cra_results, "CRA: miRNA × ICS interaction effects", highlight_mir)
ggsave(file.path(out_dir, "figures", "cra_interaction_volcano_pub.png"), p_cra_volcano_pub, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "figures", "cra_interaction_volcano_pub.pdf"), p_cra_volcano_pub, width = 7, height = 6)

# CRA miR-584-5p follow-up
cra_pred <- make_predprob_df(
  log_mat = log_quantile_counts_cra,
  pheno_df = cra_pheno_merged,
  mir_name = highlight_mir,
  outcome_col = "outcome",
  treatment_col = "treatment",
  treated_level = "2"
)

p_cra_pred <- plot_predprob(cra_pred$pred_df, "treatment", "CRA: predicted exacerbation risk by miR-584-5p and ICS")
ggsave(file.path(out_dir, "figures", "cra_mir584_predprob.png"), p_cra_pred, width = 7, height = 5, dpi = 300)

cra_pred_adj <- make_predprob_df_adjusted(
  log_mat       = log_quantile_counts_cra,
  pheno_df      = cra_pheno_merged,
  mir_name      = highlight_mir,
  outcome_col   = "outcome",
  treatment_col = "treatment",
  treated_level = "2",
  covars        = c("age.x", "gender.x")
)
p_cra_pred_adj <- plot_predprob_adjusted(
  cra_pred_adj$pred_df_adjusted, "treatment",
  "CRA: adjusted predicted hospitalization risk by miR-584-5p and ICS"
)
ggsave(file.path(out_dir, "figures", "cra_mir584_predprob_adjusted.png"),
       p_cra_pred_adj, width = 7, height = 5, dpi = 300)

cra_auc <- compute_auc_summary(
  log_mat = log_quantile_counts_cra,
  pheno_df = cra_pheno_merged,
  mir_name = highlight_mir,
  outcome_col = "outcome",
  treatment_col = "treatment",
  covars = c("pctpred_fev1_pre_BD", "`G19B_(age_onset)`", "bmi_pct", "age.x")
)

write_csv(cra_auc$stratified_auc, file.path(out_dir, "models", "cra_mir584_auc_stratified.csv"))
print(cra_auc$interaction_auc)
print(cra_auc$stratified_auc)

# CRA sensitivity analysis — explicit outcome models for miR-584-5p
cra_sens_dat <- cra_pheno_merged
cra_sens_dat$mir_expr <- as.numeric(
  log_quantile_counts_cra[highlight_mir, rownames(cra_pheno_merged)]
)
cra_sens_dat <- cra_sens_dat[
  !is.na(cra_sens_dat$Hospitalized_Asthma_Last_Yr) & !is.na(cra_sens_dat$mir_expr), ]

fit_cra_lm <- tryCatch(
  lm(Hospitalized_Asthma_Last_Yr ~ treatment * mir_expr + age.x + gender.x,
     data = cra_sens_dat),
  error = function(e) { message("CRA lm sensitivity: ", e$message); NULL }
)
if (!is.null(fit_cra_lm)) {
  write_csv(
    broom::tidy(fit_cra_lm) |> mutate(outcome = "Hospitalized_Asthma_Last_Yr"),
    file.path(out_dir, "models", "cra_linear_interaction_results.csv")
  )
}

cra_ord_n_levels <- length(unique(na.omit(cra_sens_dat$Hospitalized_Asthma_Last_Yr)))
fit_cra_polr <- NULL
if (cra_ord_n_levels > 2) {
  cra_sens_dat$outcome_ord <- factor(cra_sens_dat$Hospitalized_Asthma_Last_Yr, ordered = TRUE)
  fit_cra_polr <- tryCatch(
    MASS::polr(outcome_ord ~ treatment * mir_expr + age.x + gender.x,
               data = cra_sens_dat, Hess = TRUE),
    error = function(e) { message("CRA polr sensitivity: ", e$message); NULL }
  )
  if (!is.null(fit_cra_polr)) {
    coef_table <- coef(summary(fit_cra_polr))
    is_coef    <- !grepl("|", rownames(coef_table), fixed = TRUE)
    p_value    <- 2 * pnorm(-abs(coef_table[, "t value"]))
    polr_out   <- as.data.frame(coef_table[is_coef, , drop = FALSE]) |>
      tibble::rownames_to_column("term") |>
      mutate(p_value = p_value[is_coef], outcome = "Hospitalized_Asthma_Last_Yr_ordered")
    write_csv(polr_out, file.path(out_dir, "models", "cra_ordinal_interaction_results.csv"))
  }
}

cra_check_df <- data.frame(
  cohort = "CRA",
  analysis_type = c("linear", "ordinal"),
  selected_column = c("Hospitalized_Asthma_Last_Yr",
                      if (cra_ord_n_levels > 2) "Hospitalized_Asthma_Last_Yr" else NA_character_),
  status = c(if (!is.null(fit_cra_lm)) "fitted" else "error",
             if (!is.null(fit_cra_polr)) "fitted" else "no_predefined_ordinal_outcome"),
  reason = c("predefined: Hospitalized_Asthma_Last_Yr",
             if (cra_ord_n_levels > 2) "predefined: Hospitalized_Asthma_Last_Yr as ordered factor"
             else "no predefined ordinal exacerbation column available"),
  stringsAsFactors = FALSE
)

# Combined sensitivity column check (both cohorts)
write_csv(
  bind_rows(camp_check_df, cra_check_df),
  file.path(out_dir, "models", "sensitivity_analysis_column_check.csv")
)

# Combined predicted probability figure (publication)
p_predprob_combined <- plot_predprob_combined(
  camp_pred_adj$pred_df_adjusted, cra_pred_adj$pred_df_adjusted,
  camp_ics_level = "1", cra_ics_level = "2"
)
ggsave(file.path(out_dir, "figures", "mir584_predprob_combined_pub.png"),
       p_predprob_combined, width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "figures", "mir584_predprob_combined_pub.pdf"),
       p_predprob_combined, width = 10, height = 5)

# Shared signal across cohorts
shared_hits <- inner_join(
  cra_results |>
    dplyr::select(
      miRNA,
      beta_interaction_cra = beta_interaction,
      se_interaction_cra = se_interaction,
      OR_interaction_cra = OR_interaction,
      p_interaction_cra = p_val_interaction
    ) |> filter(p_interaction_cra < 0.05),
  camp_results |>
    dplyr::select(
      miRNA,
      beta_interaction_camp = beta_interaction,
      se_interaction_camp = se_interaction,
      OR_interaction_camp = OR_interaction,
      p_interaction_camp = p_val_interaction
    ) |> filter(p_interaction_camp < 0.05),
  by = "miRNA"
) |>
  arrange(p_interaction_cra, p_interaction_camp)

write_csv(shared_hits, file.path(out_dir, "models", "shared_miRNA_interaction_summary.csv"))
print(head(shared_hits))

# IVW meta-analysis
meta_ivw <- shared_hits |>
  mutate(
    w_camp = 1 / (se_interaction_camp^2),
    w_cra  = 1 / (se_interaction_cra^2),
    beta_meta = (
      beta_interaction_camp * w_camp +
      beta_interaction_cra  * w_cra
    ) / (w_camp + w_cra),
    se_meta = sqrt(1 / (w_camp + w_cra)),
    z_meta = beta_meta / se_meta,
    p_meta = 2 * pnorm(-abs(z_meta)),
    OR_meta = exp(beta_meta),
    lower_CI_meta = exp(beta_meta - 1.96 * se_meta),
    upper_CI_meta = exp(beta_meta + 1.96 * se_meta)
  ) |>
  arrange(p_meta)

# Sensitivity checks: p-value combination methods
meta_pval <- shared_hits |>
  mutate(
    fisher_stat = -2 * (log(p_interaction_camp) + log(p_interaction_cra)),
    fisher_p = pchisq(fisher_stat, df = 4, lower.tail = FALSE),

    z_camp_unweighted = qnorm(1 - p_interaction_camp / 2),
    z_cra_unweighted  = qnorm(1 - p_interaction_cra / 2),
    z_stouffer = (z_camp_unweighted + z_cra_unweighted) / sqrt(2),
    stouffer_p = 2 * pnorm(-abs(z_stouffer)),

    z_camp_weighted = sign(beta_interaction_camp) * qnorm(1 - p_interaction_camp / 2),
    z_cra_weighted  = sign(beta_interaction_cra)  * qnorm(1 - p_interaction_cra / 2),
    w_camp_st = 1 / se_interaction_camp,
    w_cra_st  = 1 / se_interaction_cra,
    z_stouffer_weighted = (
      w_camp_st * z_camp_weighted +
      w_cra_st  * z_cra_weighted
    ) / sqrt(w_camp_st^2 + w_cra_st^2),
    stouffer_weighted_p = 2 * pnorm(-abs(z_stouffer_weighted))
  ) |>
  dplyr::select(
    miRNA,
    fisher_p,
    stouffer_p,
    stouffer_weighted_p
  )
meta_results_all <- meta_ivw |>
  left_join(meta_pval, by = "miRNA")

print(meta_results_all)
write_csv(meta_results_all, file.path(out_dir, "models", "meta", "meta_interaction_results.csv"))