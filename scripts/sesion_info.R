# scripts/session_info.R
# Generate a reproducibility record after running the analysis pipeline

library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)
library(clusterProfiler)
library(ReactomePA)
library(multiMiR)
library(meta)
library(mediation)
library(pROC)
library(patchwork)
library(ggrepel)
library(knitr)
library(broom)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

message("Wrote sessionInfo.txt")