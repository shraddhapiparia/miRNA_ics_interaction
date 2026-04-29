options(stringsAsFactors = FALSE)

library(dplyr)
library(broom)
library(preprocessCore)
library(ggplot2)
library(stringr)
library(tidyr)

# Paths (relative to repo root)
if (!exists("out_dir")) out_dir <- "results"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- list(
  qtl_tables_dir         = file.path(out_dir, "qtl", "tables"),
  qtl_plots_dir          = file.path(out_dir, "qtl", "plots"),
  qtl_followup_dir       = file.path(out_dir, "qtl", "followup"),
  qtl_followup_plots_dir = file.path(out_dir, "qtl", "followup_plots"),
  meta_tables_dir        = file.path(out_dir, "qtl", "meta")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# CAMP input files
SNP     <- read.table("data/CAMP/camp_mir584_cis.raw", header = TRUE)
covar   <- read.csv("data/CAMP/shortened_phenotype_files_ver5.csv")
miRs    <- read.csv("data/CAMP/CAMP.492.COUNT.csv")
mapping <- read.table("data/CAMP/camp_pheno_2022-08-04.tsv", sep = "\t", header = TRUE)
pcs     <- read.table("data/CAMP/camp_pcs.txt", header = TRUE)

# CAMP miRNA normalization
miRNAs <- miRs[rowSums(miRs > 5) >= (ncol(miRs) * 0.5), ]

mir_mat <- as.matrix(miRNAs[, -1])
rownames(mir_mat) <- miRNAs$MIRNA

miR_norm <- normalize.quantiles(mir_mat)
rownames(miR_norm) <- rownames(mir_mat)
colnames(miR_norm) <- colnames(mir_mat)

miR584 <- as.data.frame(miR_norm["hsa-miR-584-5p", ])
miR584$svid <- rownames(miR584)
rownames(miR584) <- NULL
colnames(miR584)[1] <- "miRNA584"

# CAMP prepare phenotype and covariate data
mapping <- mapping[mapping$TG != 2, ]
mapping <- mapping %>% mutate(across(c(camp, U_PEDIGREEID, S_SUBJECTID, TG), as.character))
covar   <- covar %>% mutate(across(c(CAMP_ID, row.names.CAMP_492_COUNT.Trans.), as.character))
miR584  <- miR584 %>% mutate(svid = as.character(svid))

# CAMP standardize IDs
mapping_std <- mapping %>%
  rename(
    campid = camp,
    pedid  = U_PEDIGREEID,
    subid  = S_SUBJECTID,
    TG     = TG
  ) %>%
  dplyr::select(campid, pedid, subid, TG, AGE, SEX)

covar_std <- covar %>%
  rename(
    campid = CAMP_ID,
    svid   = row.names.CAMP_492_COUNT.Trans.
  )

snp2 <- SNP %>%
  rename(pedid = FID, subid = IID) %>%
  mutate(across(c(pedid, subid), as.character))

# CAMP merge data
covar_map <- covar_std %>%
  left_join(mapping_std, by = "campid")

pcs$subid <- pcs$S_SUBJECTID

df <- covar_map %>%
  inner_join(pcs, by = "subid") %>%
  inner_join(snp2, by = "subid") %>%
  inner_join(miR584, by = "svid")

df$treatment_grp <- ifelse(df$TG == 1, "ICS", "Placebo")
df$treatment_grp <- relevel(factor(df$treatment_grp), ref = "Placebo")
df$miRNA584_norm <- qnorm((rank(df$miRNA584) - 0.5) / length(df$miRNA584))

# CAMP quick checks
dim(df)
table(df$treatment_grp)
summary(df$miRNA584_norm)

# CAMP identify SNP columns
snp_cols <- setdiff(
  names(snp2),
  c("pedid", "subid", "PAT", "MAT", "SEX", "PHENOTYPE")
)
length(snp_cols)
head(snp_cols)

# CAMP base model: miRNA584 ~ SNP + covariates
results_base <- data.frame()

for (snp in snp_cols) {

  snp_bt <- paste0("`", snp, "`")

  model_formula <- as.formula(
    paste0(
      "miRNA584_norm ~ ", snp_bt,
      " + AGE.x + SEX.x",
      " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  )

  fit <- try(lm(model_formula, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) next

  tmp <- broom::tidy(fit)

  snp_row <- tmp[tmp$term == snp, ]
  if (nrow(snp_row) == 0) next

  results_base <- rbind(
    results_base,
    data.frame(
      SNP = snp,
      beta_snp = snp_row$estimate,
      p_snp = snp_row$p.value
    )
  )
}

if (nrow(results_base) > 0) {
  results_base$p_snp_fdr <- p.adjust(results_base$p_snp, method = "fdr")
  results_base <- results_base %>% arrange(p_snp)
}

head(results_base, 20)

# CAMP interaction model: miRNA584 ~ SNP + ICS + SNP*ICS + covariates
results_int <- data.frame()

for (snp in snp_cols) {

  snp_bt <- paste0("`", snp, "`")

  model_formula <- as.formula(
    paste0(
      "miRNA584_norm ~ ", snp_bt,
      " + treatment_grp + ",
      snp, ":treatment_grp",
      " + AGE.x + SEX.x",
      " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  )

  fit <- try(lm(model_formula, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) next

  tmp <- broom::tidy(fit)

  snp_row <- tmp[tmp$term == snp, ]
  int_row <- tmp[tmp$term == paste0(snp, ":treatment_grpICS"), ]

  if (nrow(snp_row) == 0) next

  if (nrow(int_row) == 0) {
    int_beta <- NA
    int_p <- NA
  } else {
    int_beta <- int_row$estimate
    int_p <- int_row$p.value
  }

  results_int <- rbind(
    results_int,
    data.frame(
      SNP = snp,
      beta_snp = snp_row$estimate,
      p_snp = snp_row$p.value,
      beta_int = int_beta,
      p_int = int_p
    )
  )
}

results_int$p_snp_fdr <- p.adjust(results_int$p_snp, method = "fdr")
results_int$p_int_fdr <- p.adjust(results_int$p_int, method = "fdr")

results_int <- results_int %>% arrange(p_int, p_snp)

head(results_int, 20)

# CAMP simple regional plot
results_base$pos <- as.numeric(sub("X5\\.(\\d+)\\..*", "\\1", results_base$SNP))
p_camp_base <- ggplot(results_base, aes(pos, -log10(p_snp))) +
  geom_point() +
  theme_bw() +
  labs(x = "chr5 position", y = "-log10(p)", title = "CAMP cis-miR-QTL base model")

ggsave(file.path(dirs[["qtl_plots_dir"]], "CAMP_base_plot.png"),
       p_camp_base, width = 7, height = 5, dpi = 300)

camp_plot_df <- results_int
camp_plot_df$index <- 1:nrow(camp_plot_df)
camp_plot_df$logp_int <- -log10(camp_plot_df$p_int)

p_camp_int <- ggplot(camp_plot_df, aes(x = index, y = logp_int)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "CAMP: MIR584 cis-mirQTL interaction scan",
    x = "SNP index across extracted region",
    y = "-log10(p interaction)"
  ) +
  theme_bw()

ggsave(file.path(dirs[["qtl_plots_dir"]], "CAMP_interaction_regional_plot.png"),
       p_camp_int, width = 7, height = 5, dpi = 300)

# CAMP save results
# NOTE: original Rmd used bare `qtl_tables_dir` (undefined); corrected to dirs[["qtl_tables_dir"]]
write.csv(results_base,
          file.path(dirs[["qtl_tables_dir"]], "CAMP_miR584_cis_mirQTL_base_model.csv"),
          row.names = FALSE)

write.csv(results_int,
          file.path(dirs[["qtl_tables_dir"]], "CAMP_miR584_cis_mirQTL_interaction_model.csv"),
          row.names = FALSE)

# CRA input files
SNP_cra     <- read.table("data/CRA/gacrs_mir584_cis.raw", header = TRUE)
covar_cra   <- read.csv("data/CRA/cra_shortened_phenotype_files10-11.csv")
miRs_cra    <- read.csv("data/CRA/2022-08-16_CRA_1159_exceRpt_miRNA_ReadCountscopy.csv", sep = ",", header = TRUE)
mapping_cra <- read.table("data/CRA/PhenotypeInfo.csv", sep = ",", header = TRUE)
cra_pcs     <- read.table("data/CRA/gacrs_pcs.txt", header = TRUE)

# CRA miRNA normalization
miRNAs_cra <- miRs_cra[rowSums(miRs_cra > 5) >= (ncol(miRs_cra) * 0.5), ]

mir_mat_cra <- as.matrix(miRNAs_cra[, -1])
rownames(mir_mat_cra) <- miRNAs_cra$X

miR_norm_cra <- normalize.quantiles(mir_mat_cra)
rownames(miR_norm_cra) <- rownames(mir_mat_cra)
colnames(miR_norm_cra) <- colnames(mir_mat_cra)

miR584_cra <- as.data.frame(miR_norm_cra["hsa-miR-584-5p", ])
miR584_cra$svid <- rownames(miR584_cra)
rownames(miR584_cra) <- NULL
colnames(miR584_cra)[1] <- "miRNA584"

# CRA prepare phenotype and covariate data
mapping_cra <- mapping_cra %>% mutate(across(c(studyid, S_SUBJECTID, Inhaled_Steroids), as.character))

covar_cra <- covar_cra %>% mutate(across(c(ST.ID, sample.id), as.character))

covar_cra$sample.id <- gsub("\\-", ".", covar_cra$sample.id)

miR584_cra <- miR584_cra %>% mutate(svid = as.character(svid))

# CRA standardize IDs
mapping_std_cra <- mapping_cra %>%
  rename(
    subid = S_SUBJECTID,
    TG    = Inhaled_Steroids
  ) %>%
  dplyr::select(studyid, subid, TG)

covar_std_cra <- covar_cra %>%
  rename(
    subid = ST.ID,
    svid  = sample.id
  )

snp2_cra <- SNP_cra %>%
  rename(pedid = FID, subid = IID) %>%
  mutate(across(c(pedid, subid), as.character)) %>%
  distinct(subid, .keep_all = TRUE)

mapping_std_cra <- mapping_std_cra %>% distinct(subid, .keep_all = TRUE)
covar_std_cra   <- covar_std_cra %>% distinct(subid, svid, .keep_all = TRUE)
miR584_cra      <- miR584_cra %>% distinct(svid, .keep_all = TRUE)

# CRA merge data
covar_map_cra <- covar_std_cra %>%
  left_join(mapping_std_cra, by = "subid")

cra_pcs$subid <- as.character(cra_pcs$IID)

df_cra <- covar_map_cra %>%
  inner_join(cra_pcs, by = "subid") %>%
  inner_join(snp2_cra, by = "subid") %>%
  inner_join(miR584_cra, by = "svid")

df_cra$treatment_grp <- ifelse(df_cra$TG == 2, "ICS", "Placebo")
df_cra$treatment_grp <- relevel(factor(df_cra$treatment_grp), ref = "Placebo")
df_cra$miRNA584_norm <- qnorm((rank(df_cra$miRNA584) - 0.5) / length(df_cra$miRNA584))

# CRA quick checks
dim(df_cra)
table(df_cra$treatment_grp)
summary(df_cra$miRNA584_norm)

# CRA identify SNP columns
snp_cols_cra <- setdiff(names(snp2_cra),
  c("pedid", "subid", "PAT", "MAT", "SEX", "PHENOTYPE"))
length(snp_cols_cra)
head(snp_cols_cra)

# CRA base model: miRNA584 ~ SNP + covariates
results_base_cra <- data.frame()

for (snp in snp_cols_cra) {

  model_formula <- as.formula(
    paste0(
      "miRNA584_norm ~ ", snp,
      " + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  )

  fit <- try(lm(model_formula, data = df_cra), silent = TRUE)
  if (inherits(fit, "try-error")) next

  tmp <- broom::tidy(fit)
  snp_row <- tmp[tmp$term == snp, ]
  if (nrow(snp_row) == 0) next

  results_base_cra <- rbind(
    results_base_cra,
    data.frame(
      SNP = snp,
      beta_snp = snp_row$estimate,
      p_snp = snp_row$p.value
    )
  )
}

results_base_cra$p_snp_fdr <- p.adjust(results_base_cra$p_snp, method = "fdr")
results_base_cra <- results_base_cra %>% arrange(p_snp)

head(results_base_cra, 20)

# CRA interaction model: miRNA584 ~ SNP + ICS + SNP*ICS + covariates
results_int_cra <- data.frame()

for (snp in snp_cols_cra) {

  model_formula <- as.formula(
    paste0(
      "miRNA584_norm ~ ", snp,
      " + treatment_grp + ",
      snp, ":treatment_grp",
      " + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  )

  fit <- try(lm(model_formula, data = df_cra), silent = TRUE)
  if (inherits(fit, "try-error")) next

  tmp <- broom::tidy(fit)

  snp_row <- tmp[tmp$term == snp, ]
  int_row <- tmp[tmp$term == paste0(snp, ":treatment_grpICS"), ]

  if (nrow(snp_row) == 0) next

  if (nrow(int_row) == 0) {
    int_beta <- NA
    int_p <- NA
  } else {
    int_beta <- int_row$estimate
    int_p <- int_row$p.value
  }

  results_int_cra <- rbind(
    results_int_cra,
    data.frame(
      SNP = snp,
      beta_snp = snp_row$estimate,
      p_snp = snp_row$p.value,
      beta_int = int_beta,
      p_int = int_p
    )
  )
}

results_int_cra$p_snp_fdr <- p.adjust(results_int_cra$p_snp, method = "fdr")
results_int_cra$p_int_fdr <- p.adjust(results_int_cra$p_int, method = "fdr")

results_int_cra <- results_int_cra %>% arrange(p_int, p_snp)

head(results_int_cra, 20)

# CRA simple regional plot
results_base_cra$pos <- as.numeric(sub("X5\\.(\\d+)\\..*", "\\1", results_base_cra$SNP))
p_cra_base <- ggplot(results_base_cra, aes(pos, -log10(p_snp))) +
  geom_point() +
  theme_bw() +
  labs(x = "chr5 position", y = "-log10(p)", title = "CRA cis-miR-QTL base model")

ggsave(file.path(dirs[["qtl_plots_dir"]], "CRA_base_plot.png"),
       p_cra_base, width = 7, height = 5, dpi = 300)

cra_plot_df <- results_int_cra
cra_plot_df$index <- 1:nrow(cra_plot_df)
cra_plot_df$logp_int <- -log10(cra_plot_df$p_int)

p_cra_int <- ggplot(cra_plot_df, aes(x = index, y = logp_int)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "CRA: MIR584 cis-mirQTL interaction scan",
    x = "SNP index across extracted region",
    y = "-log10(p interaction)"
  ) +
  theme_bw()

ggsave(file.path(dirs[["qtl_plots_dir"]], "CRA_interaction_regional_plot.png"),
       p_cra_int, width = 7, height = 5, dpi = 300)

# CRA save results
# NOTE: original Rmd used bare `qtl_tables_dir` (undefined); corrected to dirs[["qtl_tables_dir"]]
write.csv(results_base_cra,
          file.path(dirs[["qtl_tables_dir"]], "CRA_miR584_cis_mirQTL_base_model.csv"),
          row.names = FALSE)

write.csv(results_int_cra,
          file.path(dirs[["qtl_tables_dir"]], "CRA_miR584_cis_mirQTL_interaction_model.csv"),
          row.names = FALSE)

# Notes
cat("CAMP base model:\n")
cat("miRNA584 ~ SNP + AGE + SEX + PCs\n\n")

cat("CAMP interaction model:\n")
cat("miRNA584 ~ SNP + treatment_grp + SNP*treatment_grp + AGE + SEX + PCs\n\n")

cat("CRA base model:\n")
cat("miRNA584 ~ SNP + age + gender\n\n")

cat("CRA interaction model:\n")
cat("miRNA584 ~ SNP + treatment_grp + SNP*treatment_grp + age + gender\n\n")

cat("Interpretation:\n")
cat("- base model tests overall SNP association with miR-584-5p\n")
cat("- interaction model tests whether SNP effect differs by ICS status\n")
cat("- focus on p_snp_fdr in base model and p_int_fdr in interaction model\n")
cat("- simple plots show signal pattern across the extracted regional SNP set\n")

top_camp <- results_base %>% slice_min(p_snp, n = 50)
top_cra  <- results_base_cra %>% slice_min(p_snp, n = 50)

intersect(top_camp$SNP, top_cra$SNP)

top_camp_pos <- results_base %>%
  slice_min(p_snp, n = 50) %>%
  mutate(
    chr = str_extract(SNP, "(?<=X)\\d+(?=\\.)"),
    pos = str_extract(SNP, "(?<=X\\d\\.)\\d+")
  ) %>%
  mutate(
    chr = as.integer(chr),
    pos = as.integer(pos)
  )

top_cra_pos <- results_base_cra %>%
  slice_min(p_snp, n = 50) %>%
  mutate(
    chr = str_extract(SNP, "(?<=X)\\d+(?=\\.)"),
    pos = str_extract(SNP, "(?<=X\\d\\.)\\d+")
  ) %>%
  mutate(
    chr = as.integer(chr),
    pos = as.integer(pos)
  )

overlap_snps <- inner_join(
  top_camp_pos %>% dplyr::select(SNP_camp = SNP, chr, pos, p_camp = p_snp),
  top_cra_pos  %>% dplyr::select(SNP_cra  = SNP, chr, pos, p_cra  = p_snp),
  by = c("chr", "pos")
)
print(overlap_snps)

camp_hits <- results_base %>%
  filter(pos %in% overlap_snps$pos) %>%
  dplyr::select(SNP, pos, beta_snp, p_snp, p_snp_fdr)

cra_hits <- results_base_cra %>%
  filter(pos %in% overlap_snps$pos) %>%
  dplyr::select(SNP, pos, beta_snp, p_snp, p_snp_fdr)

combined_df <- inner_join(
  camp_hits,
  cra_hits,
  by = "pos",
  suffix = c("_camp", "_cra")
)
print(combined_df)

results_int <- results_int %>%
  mutate(
    chr = as.integer(str_extract(SNP, "(?<=X)\\d+(?=\\.)")),
    pos = as.integer(str_extract(SNP, "(?<=X\\d\\.)\\d+"))
  )

results_int_cra <- results_int_cra %>%
  mutate(
    chr = as.integer(str_extract(SNP, "(?<=X)\\d+(?=\\.)")),
    pos = as.integer(str_extract(SNP, "(?<=X\\d\\.)\\d+"))
  )

camp_int_leads <- results_int %>%
  filter(pos %in% overlap_snps$pos) %>%
  arrange(match(pos, overlap_snps$pos))

cra_int_leads <- results_int_cra %>%
  filter(pos %in% overlap_snps$pos) %>%
  arrange(match(pos, overlap_snps$pos))

print(camp_int_leads)
print(cra_int_leads)

# Define lead SNPs: top 3 base-model hits by p_snp, one per cohort
if (!"SNP"   %in% names(results_base)) stop("results_base has no SNP column")
if (!"p_snp" %in% names(results_base)) stop("results_base has no p_snp column")

camp_leads <- results_base %>%
  arrange(p_snp) %>%
  slice_head(n = 3) %>%
  dplyr::select(SNP, beta_snp, p_snp, p_snp_fdr, pos)

if (nrow(camp_leads) == 0) stop("camp_leads has zero rows — check results_base")
missing_camp <- setdiff(camp_leads$SNP, names(df))
if (length(missing_camp) > 0) stop("camp_leads SNPs not found in df: ", paste(missing_camp, collapse = ", "))
message("camp_leads SNPs: ", paste(camp_leads$SNP, collapse = ", "))

if (!"SNP"   %in% names(results_base_cra)) stop("results_base_cra has no SNP column")
if (!"p_snp" %in% names(results_base_cra)) stop("results_base_cra has no p_snp column")

cra_leads <- results_base_cra %>%
  arrange(p_snp) %>%
  slice_head(n = 3) %>%
  dplyr::select(SNP, beta_snp, p_snp, p_snp_fdr, pos)

if (nrow(cra_leads) == 0) stop("cra_leads has zero rows — check results_base_cra")
missing_cra <- setdiff(cra_leads$SNP, names(df_cra))
if (length(missing_cra) > 0) stop("cra_leads SNPs not found in df_cra: ", paste(missing_cra, collapse = ", "))
message("cra_leads SNPs: ", paste(cra_leads$SNP, collapse = ", "))

camp_long <- df %>%
  dplyr::select(miRNA584_norm, treatment_grp, all_of(camp_leads$SNP)) %>%
  pivot_longer(
    cols = all_of(camp_leads$SNP),
    names_to = "SNP",
    values_to = "Genotype"
  )

ggplot(camp_long, aes(x = factor(Genotype), y = miRNA584_norm)) +
  geom_boxplot() +
  facet_wrap(~SNP, scales = "free_x") +
  theme_bw() +
  labs(
    title = "CAMP: miR-584-5p vs genotype",
    x = "Genotype",
    y = "miR-584-5p expression"
  )

ggplot(camp_long, aes(x = factor(Genotype), y = miRNA584_norm, color = treatment_grp)) +
  geom_boxplot() +
  facet_wrap(~SNP, scales = "free_x") +
  theme_bw() +
  labs(
    title = "CAMP: miR-584-5p vs genotype by ICS",
    x = "Genotype",
    y = "miR-584-5p expression",
    color = "ICS group"
  )

cra_long <- df_cra %>%
  dplyr::select(miRNA584_norm, treatment_grp, all_of(cra_leads$SNP)) %>%
  pivot_longer(
    cols = all_of(cra_leads$SNP),
    names_to = "SNP",
    values_to = "Genotype"
  )

ggplot(cra_long, aes(x = factor(Genotype), y = miRNA584_norm)) +
  geom_boxplot() +
  facet_wrap(~SNP, scales = "free_x") +
  theme_bw() +
  labs(
    title = "CRA: miR-584-5p vs genotype",
    x = "Genotype",
    y = "miR-584-5p expression"
  )

ggplot(cra_long, aes(x = factor(Genotype), y = miRNA584_norm, color = treatment_grp)) +
  geom_boxplot() +
  facet_wrap(~SNP, scales = "free_x") +
  theme_bw() +
  labs(
    title = "CAMP: miR-584-5p vs genotype by ICS",
    x = "Genotype",
    y = "miR-584-5p expression",
    color = "ICS group"
  )

camp_snps <- camp_leads$SNP

cond_results_camp <- data.frame()

pairs_camp <- combn(camp_snps, 2, simplify = FALSE)

for (pair in pairs_camp) {

  snp1 <- pair[1]
  snp2 <- pair[2]

  snp1_bt <- paste0("`", snp1, "`")
  snp2_bt <- paste0("`", snp2, "`")

  model_formula <- as.formula(
    paste0(
      "miRNA584_norm ~ ", snp1_bt, " + ", snp2_bt,
      " + AGE.x + SEX.x",
      " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  )

  fit <- lm(model_formula, data = df)
  tmp <- broom::tidy(fit)

  row1 <- tmp[tmp$term == snp1, ]
  row2 <- tmp[tmp$term == snp2, ]

  cond_results_camp <- rbind(
    cond_results_camp,
    data.frame(
      snp1 = snp1,
      snp2 = snp2,
      beta_snp1 = row1$estimate,
      p_snp1 = row1$p.value,
      beta_snp2 = row2$estimate,
      p_snp2 = row2$p.value
    )
  )
}

print(cond_results_camp)

s1 <- paste0("`", camp_snps[1], "`")
s2 <- paste0("`", camp_snps[2], "`")
s3 <- paste0("`", camp_snps[3], "`")

fit3_camp <- lm(
  as.formula(
    paste0(
      "miRNA584_norm ~ ", s1, " + ", s2, " + ", s3,
      " + AGE.x + SEX.x",
      " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  ),
  data = df
)

print(broom::tidy(fit3_camp))

cra_snps <- cra_leads$SNP

cond_results_cra <- data.frame()

pairs_cra <- combn(cra_snps, 2, simplify = FALSE)

for (pair in pairs_cra) {

  snp1 <- pair[1]
  snp2 <- pair[2]

  snp1_bt <- paste0("`", snp1, "`")
  snp2_bt <- paste0("`", snp2, "`")

  model_formula <- as.formula(
    paste0(
      "miRNA584_norm ~ ", snp1_bt, " + ", snp2_bt,
      " + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  )

  fit <- lm(model_formula, data = df_cra)
  tmp <- broom::tidy(fit)

  row1 <- tmp[tmp$term == snp1, ]
  row2 <- tmp[tmp$term == snp2, ]

  cond_results_cra <- rbind(
    cond_results_cra,
    data.frame(
      snp1 = snp1,
      snp2 = snp2,
      beta_snp1 = row1$estimate,
      p_snp1 = row1$p.value,
      beta_snp2 = row2$estimate,
      p_snp2 = row2$p.value
    )
  )
}

print(cond_results_cra)

s1 <- paste0("`", cra_snps[1], "`")
s2 <- paste0("`", cra_snps[2], "`")
s3 <- paste0("`", cra_snps[3], "`")

fit3_cra <- lm(
  as.formula(
    paste0(
      "miRNA584_norm ~ ", s1, " + ", s2, " + ", s3,
      " + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
    )
  ),
  data = df_cra
)

print(broom::tidy(fit3_cra))

# Meta-analysis: merge cohorts by genomic position
replicated_hits <- inner_join(
  camp_int_leads %>%
    dplyr::select(SNP_camp = SNP, chr, pos, beta_int_camp = beta_int, p_int_camp = p_int),

  cra_int_leads %>%
    dplyr::select(SNP_cra = SNP, chr, pos, beta_int_cra = beta_int, p_int_cra = p_int),

  by = c("chr", "pos")
)

replicated_hits <- replicated_hits %>%
  filter(sign(beta_int_camp) == sign(beta_int_cra))

replicated_hits <- replicated_hits %>%
  arrange(p_int_camp)

replicated_hits <- replicated_hits %>%
  filter(p_int_camp < 0.05 | p_int_cra < 0.05)

meta_input <- replicated_hits %>%
  mutate(
    z_camp = qnorm(p_int_camp / 2, lower.tail = FALSE),
    z_cra  = qnorm(p_int_cra  / 2, lower.tail = FALSE),

    se_camp = abs(beta_int_camp) / z_camp,
    se_cra  = abs(beta_int_cra)  / z_cra,

    var_camp = se_camp^2,
    var_cra  = se_cra^2,

    w_camp = 1 / var_camp,
    w_cra  = 1 / var_cra,

    beta_meta = (beta_int_camp * w_camp + beta_int_cra * w_cra) / (w_camp + w_cra),
    se_meta   = sqrt(1 / (w_camp + w_cra)),

    z_meta = beta_meta / se_meta,
    p_meta = 2 * pnorm(abs(z_meta), lower.tail = FALSE),

    q_stat = w_camp * (beta_int_camp - beta_meta)^2 +
             w_cra  * (beta_int_cra  - beta_meta)^2,

    q_df = 1,
    p_het = pchisq(q_stat, df = q_df, lower.tail = FALSE),

    i2 = pmax(0, (q_stat - q_df) / q_stat) * 100
  ) %>%
  arrange(p_meta)

write.csv(meta_input,
          file.path(dirs[["meta_tables_dir"]], "CAMP_CRA_leadSNP_meta_analysis.csv"),
          row.names = FALSE)

print(meta_input)

# NOTE: original Rmd used bare `qtl_followup_tables_dir` (undefined); corrected to dirs[["qtl_followup_dir"]]
write.csv(top_camp_pos,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_top50_base_hits_with_position.csv"),
          row.names = FALSE)

write.csv(top_cra_pos,
          file.path(dirs[["qtl_followup_dir"]], "CRA_top50_base_hits_with_position.csv"),
          row.names = FALSE)

write.csv(overlap_snps,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_CRA_overlap_top50_by_position.csv"),
          row.names = FALSE)

write.csv(camp_hits,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_overlap_hits.csv"),
          row.names = FALSE)

write.csv(cra_hits,
          file.path(dirs[["qtl_followup_dir"]], "CRA_overlap_hits.csv"),
          row.names = FALSE)

write.csv(replicated_hits,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_CRA_replicated_hits.csv"),
          row.names = FALSE)

write.csv(camp_leads,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_lead3_base_hits.csv"),
          row.names = FALSE)

write.csv(cra_leads,
          file.path(dirs[["qtl_followup_dir"]], "CRA_lead3_base_hits.csv"),
          row.names = FALSE)

write.csv(camp_int_leads,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_lead3_interaction_hits.csv"),
          row.names = FALSE)

write.csv(cra_int_leads,
          file.path(dirs[["qtl_followup_dir"]], "CRA_lead3_interaction_hits.csv"),
          row.names = FALSE)

write.csv(cond_results_camp,
          file.path(dirs[["qtl_followup_dir"]], "CAMP_pairwise_conditional_models.csv"),
          row.names = FALSE)

write.csv(cond_results_cra,
          file.path(dirs[["qtl_followup_dir"]], "CRA_pairwise_conditional_models.csv"),
          row.names = FALSE)

write.csv(broom::tidy(fit3_camp),
          file.path(dirs[["qtl_followup_dir"]], "CAMP_three_snp_model.csv"),
          row.names = FALSE)

write.csv(broom::tidy(fit3_cra),
          file.path(dirs[["qtl_followup_dir"]], "CRA_three_snp_model.csv"),
          row.names = FALSE)
