options(stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(preprocessCore)
library(broom)
library(mediation)
library(stringr)

# Paths (relative to repo root)
cra_counts_path <- "data/CRA/2022-08-16_CRA_1159_exceRpt_miRNA_ReadCountscopy.csv"
cra_pheno_path  <- "data/CRA/PhenotypeInfo.csv"
cra_link_path   <- "data/CRA/cra_shortened_phenotype_files10-11.csv"
cra_snp_path    <- "data/CRA/gacrs_mir584_cis.raw"
cra_pcs_path    <- "data/CRA/gacrs_pcs.txt"
camp_counts_path <- "data/CAMP/CAMP.492.COUNT.csv"
camp_pheno_path  <- "data/CAMP/camp_pheno_2022-08-04.tsv"
camp_link_path   <- "data/CAMP/shortened_phenotype_files_ver5.csv"
camp_snp_path    <- "data/CAMP/camp_mir584_cis.raw"
camp_pcs_path    <- "data/CAMP/camp_pcs.txt"
if (!exists("out_dir")) out_dir <- "results"
highlight_mir    <- "hsa-miR-584-5p"
min_count        <- 5
min_prop_samples <- 0.5

med_tables_dir <- file.path(out_dir, "mediation", "tables")
med_plots_dir  <- file.path(out_dir, "mediation", "plots")

dir.create(med_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(med_plots_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: quantile-normalize one miRNA
extract_miR584_quantile <- function(count_df, mir_col, mir_name = "hsa-miR-584-5p",
                                    min_count = 5, min_prop = 0.5) {
  keep <- rowSums(count_df[, -1] > min_count) >= (ncol(count_df[, -1]) * min_prop)
  miRNAs <- count_df[keep, ]

  mir_mat <- as.matrix(miRNAs[, -1])
  rownames(mir_mat) <- miRNAs[[mir_col]]

  mir_norm <- normalize.quantiles(mir_mat)
  rownames(mir_norm) <- rownames(mir_mat)
  colnames(mir_norm) <- colnames(mir_mat)

  out <- as.data.frame(mir_norm[mir_name, ])
  out$svid <- rownames(out)
  rownames(out) <- NULL
  colnames(out)[1] <- "miRNA584"

  out
}

# CAMP input files
camp_snp    <- read.table(camp_snp_path, header = TRUE)
camp_link   <- read.csv(camp_link_path)
camp_counts <- read.csv(camp_counts_path)
camp_pheno  <- read.table(camp_pheno_path, sep = "\t", header = TRUE)
camp_pcs    <- read.table(camp_pcs_path, header = TRUE)

# CAMP miRNA normalization
camp_miR584 <- extract_miR584_quantile(
  count_df  = camp_counts,
  mir_col   = "MIRNA",
  mir_name  = highlight_mir,
  min_count = min_count,
  min_prop  = min_prop_samples
)

# CAMP prepare and merge
camp_pheno2 <- camp_pheno %>%
  filter(TG != 2) %>%
  mutate(across(c(camp, U_PEDIGREEID, S_SUBJECTID, TG), as.character))

camp_link2 <- camp_link %>%
  mutate(across(c(CAMP_ID, row.names.CAMP_492_COUNT.Trans.), as.character))

camp_miR584 <- camp_miR584 %>%
  mutate(svid = as.character(svid))
colnames(camp_snp) <- trimws(colnames(camp_snp))
camp_snp2 <- camp_snp %>%
  dplyr::rename(pedid = FID, subid = IID) %>%
  mutate(across(c(pedid, subid), as.character))

# Verify all required columns exist in camp_pheno2 before rename+select
.required_camp_pheno <- c("camp", "U_PEDIGREEID", "S_SUBJECTID", "TG", "AGE", "SEX", "EDHOS_cum_Y1")
.missing_camp_pheno  <- setdiff(.required_camp_pheno, names(camp_pheno2))
if (length(.missing_camp_pheno) > 0) {
  message("Available columns in camp_pheno2: ", paste(names(camp_pheno2), collapse = ", "))
  stop("camp_map select: required column(s) missing from camp_pheno2: ",
       paste(.missing_camp_pheno, collapse = ", "))
}
rm(.required_camp_pheno, .missing_camp_pheno)

camp_map <- camp_pheno2 %>%
  dplyr::rename(
    campid = camp,
    pedid = U_PEDIGREEID,
    subid = S_SUBJECTID
  ) %>%
  dplyr::select(campid, pedid, subid, TG, AGE, SEX, EDHOS_cum_Y1)

camp_link_std <- camp_link2 %>%
  dplyr::rename(
    campid = CAMP_ID,
    svid = row.names.CAMP_492_COUNT.Trans.
  )

camp_pcs$subid <- as.character(camp_pcs$S_SUBJECTID)

camp_df <- camp_link_std %>%
  left_join(camp_map, by = "campid") %>%
  inner_join(camp_pcs, by = "subid") %>%
  inner_join(camp_snp2, by = "subid") %>%
  inner_join(camp_miR584, by = "svid")

camp_df$treatment_grp <- ifelse(camp_df$TG == 1, "ICS", "Placebo")
camp_df$treatment_grp <- relevel(factor(camp_df$treatment_grp), ref = "Placebo")

camp_df$miRNA584_norm <- qnorm((rank(camp_df$miRNA584) - 0.5) / length(camp_df$miRNA584))

# CAMP define outcome and SNP
lead_snp_camp <- c("X5.149061758.C.T_C", "X5.149031357.A.C_A", "X5.149121003.G.A_A")

camp_df <- camp_df |>
  filter(!is.na(EDHOS_cum_Y1)) |>
  mutate(
    outcome = factor(ifelse(EDHOS_cum_Y1 < 1, "low", "high")),
    outcome = relevel(outcome, ref = "low"),
    treatment = factor(TG),
    treatment = relevel(treatment, ref = "3"),
    SEX.x = factor(SEX.x)
  )

if ("AGE.x" %in% names(camp_df) && !"AGE" %in% names(camp_df)) camp_df$AGE <- camp_df$AGE.x
if ("SEX.x" %in% names(camp_df) && !"SEX" %in% names(camp_df)) camp_df$SEX <- camp_df$SEX.x
camp_df$SEX <- factor(camp_df$SEX)

# CAMP SNP -> exacerbation check
results_list <- list()

for (snp in lead_snp_camp) {

  formula_str <- paste0(
    "outcome ~ `", snp, "` + AGE + SEX + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  )

  model <- glm(
    as.formula(formula_str),
    data = camp_df,
    family = binomial()
  )

  res <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == snp) %>%
    mutate(SNP = snp)

  results_list[[snp]] <- res
}

final_results <- bind_rows(results_list)
print(final_results)

# Named model for primary lead SNP (rs36047 = X5.149061758.C.T_C), used in save block
camp_snp_outcome <- glm(
  outcome ~ `X5.149061758.C.T_C` + AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df,
  family = binomial()
)
camp_out_fit_primary <- camp_snp_outcome

# CAMP SNP -> exacerbation check for ICS subgroup
results_list <- list()
ics_df <- camp_df[camp_df$TG == 1, ]
for (snp in lead_snp_camp) {

  formula_str <- paste0(
    "outcome ~ `", snp, "` + AGE + SEX + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  )

  model <- glm(
    as.formula(formula_str),
    data = ics_df,
    family = binomial()
  )

  res <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == snp) %>%
    mutate(SNP = snp)

  results_list[[snp]] <- res
}

ics_results <- bind_rows(results_list)
print(ics_results)

# CAMP SNP -> miRNA mediator model
results_list <- list()
covars <- c("AGE", "SEX", paste0("PC", 1:10))

for (snp in lead_snp_camp) {

  formula_str <- paste0(
    "miRNA584_norm ~ `", snp, "` + AGE + SEX + treatment_grp + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  )

  camp_med_fit <- lm(
    as.formula(formula_str),
    data = camp_df
  )

  res <- broom::tidy(camp_med_fit, conf.int = TRUE) %>%
    filter(!term %in% covars) %>%
    mutate(SNP_tested = snp)

  results_list[[snp]] <- res
}

final_results <- bind_rows(results_list)
print(final_results)

# Explicit mediator model for primary SNP (rs36047 = X5.149061758.C.T_C)
camp_med_fit_primary <- lm(
  miRNA584_norm ~ `X5.149061758.C.T_C` + AGE + SEX + treatment_grp +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df
)

# CAMP miRNA -> exacerbation outcome model with interaction
covars <- c("AGE", "SEX", paste0("PC", 1:10))

results_list <- list()
model_list <- list()

for (snp in lead_snp_camp) {

  formula_str <- paste0(
    "outcome ~ `", snp, "` + miRNA584_norm * treatment_grp + ",
    "AGE + SEX + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  )

  fit <- glm(
    as.formula(formula_str),
    data = camp_df,
    family = binomial()
  )

  model_list[[snp]] <- fit

  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(!term %in% covars) %>%
    mutate(SNP_tested = snp)

  results_list[[snp]] <- res
}

final_results_int <- bind_rows(results_list)
print(final_results_int)

# Named model for primary lead SNP (rs36047), used in save block
camp_out_fit_int <- model_list[["X5.149061758.C.T_C"]]
if (is.null(camp_out_fit_int)) stop("camp_out_fit_int: model for X5.149061758.C.T_C not found in model_list")
camp_out_fit_int_primary <- camp_out_fit_int

# CAMP pooled mediation without interaction
camp_out_fit_pooled <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm +
    AGE + SEX + treatment_grp +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df,
  family = binomial()
)

camp_med_pooled <- mediate(
  model.m = camp_med_fit_primary,
  model.y = camp_out_fit_pooled,
  treat = "X5.149061758.C.T_C",
  mediator = "miRNA584_norm",
  boot = TRUE,
  sims = 1000
)

summary(camp_med_pooled)

# CAMP stratified mediation: ICS only
camp_df_ics <- camp_df %>% filter(treatment_grp == "ICS")

camp_med_fit_ics <- lm(
  miRNA584_norm ~ `X5.149061758.C.T_C` + AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df_ics
)

camp_out_fit_ics <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm +
    AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df_ics,
  family = binomial()
)

camp_med_ics <- mediate(
  model.m = camp_med_fit_ics,
  model.y = camp_out_fit_ics,
  treat = "X5.149061758.C.T_C",
  mediator = "miRNA584_norm",
  boot = TRUE,
  sims = 1000
)

summary(camp_med_ics)

# CAMP stratified mediation: Placebo only
camp_df_placebo <- camp_df %>% filter(treatment_grp == "Placebo")

camp_med_fit_placebo <- lm(
  miRNA584_norm ~ `X5.149061758.C.T_C` + AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df_placebo
)

camp_out_fit_placebo <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm +
    AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df_placebo,
  family = binomial()
)

camp_med_placebo <- mediate(
  model.m = camp_med_fit_placebo,
  model.y = camp_out_fit_placebo,
  treat = "X5.149061758.C.T_C",
  mediator = "miRNA584_norm",
  boot = TRUE,
  sims = 1000
)

summary(camp_med_placebo)

# CRA / GACRS

# CRA input files
cra_snp    <- read.table(cra_snp_path, header = TRUE)
cra_link   <- read.csv(cra_link_path)
cra_counts <- read.csv(cra_counts_path)
cra_pheno  <- read.table(cra_pheno_path, sep = ",", header = TRUE)
cra_pcs    <- read.table(cra_pcs_path, header = TRUE)

# CRA miRNA normalization
cra_miR584 <- extract_miR584_quantile(
  count_df  = cra_counts,
  mir_col   = "X",
  mir_name  = highlight_mir,
  min_count = min_count,
  min_prop  = min_prop_samples
)

# CRA prepare and merge
cra_pheno2 <- cra_pheno %>%
  mutate(across(c(studyid, S_SUBJECTID, Inhaled_Steroids), as.character))

cra_link2 <- cra_link %>%
  mutate(across(c(ST.ID, sample.id), as.character))

cra_link2$sample.id <- gsub("\\-", ".", cra_link2$sample.id)

cra_miR584 <- cra_miR584 %>%
  mutate(svid = as.character(svid))

cra_snp2 <- cra_snp %>%
  rename(pedid = FID, subid = IID) %>%
  mutate(across(c(pedid, subid), as.character)) %>%
  distinct(subid, .keep_all = TRUE)

# Verify all required columns exist in cra_pheno2 before rename+select
# TG is Inhaled_Steroids after rename; check pre-rename names here
.required_cra_pheno <- c("studyid", "S_SUBJECTID", "Inhaled_Steroids", "Hospitalized_Asthma_Last_Yr")
.missing_cra_pheno  <- setdiff(.required_cra_pheno, names(cra_pheno2))
if (length(.missing_cra_pheno) > 0) {
  message("Available columns in cra_pheno2: ", paste(names(cra_pheno2), collapse = ", "))
  stop("cra_map select: required column(s) missing from cra_pheno2: ",
       paste(.missing_cra_pheno, collapse = ", "))
}
rm(.required_cra_pheno, .missing_cra_pheno)

cra_map <- cra_pheno2 %>%
  dplyr::rename(
    subid = S_SUBJECTID,
    TG = Inhaled_Steroids
  ) %>%
  dplyr::select(studyid, subid, TG, Hospitalized_Asthma_Last_Yr)

cra_link_std <- cra_link2 %>%
  rename(
    subid = ST.ID,
    svid = sample.id
  ) %>%
  distinct(subid, svid, .keep_all = TRUE)

cra_pcs$subid <- as.character(cra_pcs$IID)

cra_df <- cra_link_std %>%
  left_join(cra_map, by = "subid") %>%
  inner_join(cra_snp2, by = "subid") %>%
  inner_join(cra_pcs, by = "subid") %>%
  inner_join(cra_miR584, by = "svid")

cra_df$treatment_grp <- ifelse(cra_df$TG == 2, "ICS", "Placebo")
cra_df$treatment_grp <- relevel(factor(cra_df$treatment_grp), ref = "Placebo")

cra_df$miRNA584_norm <- qnorm((rank(cra_df$miRNA584) - 0.5) / length(cra_df$miRNA584))

# CRA define outcome and SNP
lead_snp_cra <- "X5.149061758.C.T_C"

cra_df <- cra_df |>
  filter(!is.na(Hospitalized_Asthma_Last_Yr)) |>
  mutate(
    outcome = factor(Hospitalized_Asthma_Last_Yr),
    treatment = factor(TG),
    treatment = relevel(treatment, ref = "1"),
    gender.x = factor(gender)
  )

if ("gender.x" %in% names(cra_df) && !"gender" %in% names(cra_df)) cra_df$gender <- cra_df$gender.x
cra_df$gender <- factor(cra_df$gender)

stopifnot(lead_snp_cra %in% names(cra_df))

# CRA SNP -> exacerbation check
cra_snp_outcome <- glm(
  outcome ~ `X5.149061758.C.T_C` + age + gender +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df,
  family = binomial()
)

summary(cra_snp_outcome)
print(broom::tidy(cra_snp_outcome, exponentiate = TRUE, conf.int = TRUE))

# CRA SNP -> miRNA mediator model
cra_med_fit <- lm(
  miRNA584_norm ~ `X5.149061758.C.T_C` + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + treatment_grp,
  data = cra_df
)

summary(cra_med_fit)
print(broom::tidy(cra_med_fit))

# CRA miRNA -> exacerbation outcome model with interaction
cra_out_fit_int <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm * treatment_grp + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df,
  family = binomial()
)

summary(cra_out_fit_int)
print(broom::tidy(cra_out_fit_int, exponentiate = TRUE, conf.int = TRUE))

# CRA pooled mediation without interaction
cra_out_fit_pooled <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm + treatment_grp + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df,
  family = binomial()
)

cra_med_pooled <- mediate(
  model.m = cra_med_fit,
  model.y = cra_out_fit_pooled,
  treat = "X5.149061758.C.T_C",
  mediator = "miRNA584_norm",
  boot = TRUE,
  sims = 1000
)

summary(cra_med_pooled)

# CRA stratified mediation: ICS only
cra_df_ics <- cra_df %>% filter(treatment_grp == "ICS")

cra_med_fit_ics <- lm(
  miRNA584_norm ~ `X5.149061758.C.T_C` + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df_ics
)

cra_out_fit_ics <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df_ics,
  family = binomial()
)

cra_med_ics <- mediate(
  model.m = cra_med_fit_ics,
  model.y = cra_out_fit_ics,
  treat = "X5.149061758.C.T_C",
  mediator = "miRNA584_norm",
  boot = TRUE,
  sims = 1000
)

summary(cra_med_ics)

# CRA stratified mediation: Placebo only
cra_df_placebo <- cra_df %>% filter(treatment_grp == "Placebo")

cra_med_fit_placebo <- lm(
  miRNA584_norm ~ `X5.149061758.C.T_C` + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df_placebo
)

cra_out_fit_placebo <- glm(
  outcome ~ `X5.149061758.C.T_C` + miRNA584_norm + age + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df_placebo,
  family = binomial()
)

cra_med_placebo <- mediate(
  model.m = cra_med_fit_placebo,
  model.y = cra_out_fit_placebo,
  treat = "X5.149061758.C.T_C",
  mediator = "miRNA584_norm",
  boot = TRUE,
  sims = 1000
)

summary(cra_med_placebo)

# Stratified total effect check: SNP -> outcome without mediator
camp_interaction_check <- glm(
  outcome ~ `X5.149061758.C.T_C` * treatment_grp + AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df,
  family = binomial()
)

camp_total_ics <- glm(
  outcome ~ `X5.149061758.C.T_C` + AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df[camp_df$treatment_grp == "ICS", ],
  family = binomial()
)

camp_total_placebo <- glm(
  outcome ~ `X5.149061758.C.T_C` + AGE + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = camp_df[camp_df$treatment_grp == "Placebo", ],
  family = binomial()
)

cra_interaction_check <- glm(
  outcome ~ `X5.149061758.C.T_C` * treatment_grp + age + gender +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df,
  family = binomial()
)

cra_total_ics <- glm(
  outcome ~ `X5.149061758.C.T_C` + age + gender +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df[cra_df$treatment_grp == "ICS", ],
  family = binomial()
)

cra_total_placebo <- glm(
  outcome ~ `X5.149061758.C.T_C` + age + gender +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = cra_df[cra_df$treatment_grp == "Placebo", ],
  family = binomial()
)

print(list(
  interaction  = broom::tidy(camp_interaction_check, exponentiate = TRUE),
  pooled  = broom::tidy(camp_out_fit_primary, exponentiate = TRUE),
  ics     = broom::tidy(camp_total_ics,        exponentiate = TRUE),
  placebo = broom::tidy(camp_total_placebo,    exponentiate = TRUE)
))

print(list(
  interaction  = broom::tidy(cra_interaction_check, exponentiate = TRUE),
  pooled  = broom::tidy(cra_snp_outcome, exponentiate = TRUE),
  cra_ics     = broom::tidy(cra_total_ics,     exponentiate = TRUE),
  cra_placebo = broom::tidy(cra_total_placebo, exponentiate = TRUE)
))

# Helper: tidy a model and tag with cohort + model name
tidy_with_model <- function(model, model_name, cohort, exponentiate = FALSE) {
  broom::tidy(model, exponentiate = exponentiate, conf.int = TRUE) %>%
    mutate(cohort = cohort, model = model_name, .before = 1)
}

# Save consolidated outputs
message("CAMP rs36047 outputs use primary SNP: X5.149061758.C.T_C ",
        "(camp_out_fit_primary, camp_med_fit_primary, camp_out_fit_int_primary)")

camp_mediation_models <- bind_rows(
  tidy_with_model(camp_out_fit_primary,    "snp_to_exacerbation_pooled",           "CAMP", TRUE),
  tidy_with_model(camp_med_fit_primary,    "snp_to_miRNA_mediator",                "CAMP", FALSE),
  tidy_with_model(camp_out_fit_int_primary,"miRNA_treatment_interaction_outcome",   "CAMP", TRUE),
  tidy_with_model(camp_interaction_check,  "snp_treatment_interaction_total_effect","CAMP", TRUE),
  tidy_with_model(camp_total_ics,          "total_effect_ICS",                     "CAMP", TRUE),
  tidy_with_model(camp_total_placebo,      "total_effect_placebo",                 "CAMP", TRUE)
)

write.csv(camp_mediation_models,
          file.path(med_tables_dir, "CAMP_rs36047_mediation_models.csv"),
          row.names = FALSE)

cra_mediation_models <- bind_rows(
  tidy_with_model(cra_snp_outcome, "snp_to_exacerbation_pooled",          "CRA", TRUE),
  tidy_with_model(cra_med_fit,     "snp_to_miRNA_mediator",               "CRA", FALSE),
  tidy_with_model(cra_out_fit_int,       "miRNA_treatment_interaction_outcome",   "CRA", TRUE),
  tidy_with_model(cra_interaction_check, "snp_treatment_interaction_total_effect", "CRA", TRUE),
  tidy_with_model(cra_total_ics,         "total_effect_ICS",                       "CRA", TRUE),
  tidy_with_model(cra_total_placebo,"total_effect_placebo",               "CRA", TRUE)
)

write.csv(cra_mediation_models,
          file.path(med_tables_dir, "CRA_rs36047_mediation_models.csv"),
          row.names = FALSE)

write.csv(
  bind_rows(camp_mediation_models, cra_mediation_models),
  file.path(med_tables_dir, "combined_rs36047_mediation_models.csv"),
  row.names = FALSE
)

# Summary objects from mediation output
camp_med_ics_sum     <- summary(camp_med_ics)
camp_med_placebo_sum <- summary(camp_med_placebo)
cra_med_ics_sum      <- summary(cra_med_ics)
cra_med_placebo_sum  <- summary(cra_med_placebo)
