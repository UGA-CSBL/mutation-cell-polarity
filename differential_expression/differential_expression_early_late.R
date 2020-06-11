library(tidyverse)
library(edgeR)

rm(list = ls())
set.seed(1005)

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-READ", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

# Perform DESeq2 on input projects with missing RData files ------------------
annot <- data.table::fread("~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv") %>%
  filter(project %in% projects) %>%
  mutate(tumor_stage = case_when(
    sample_type == "Solid Tissue Normal" ~ "Normal",
    str_detect(tumor_stage, "(^i$)|(\\si[abc]?$)|(1)") ~ "Stage I",
    str_detect(tumor_stage, "(^ii$)|(\\si{2}[abc]?$)|(2)") ~ "Stage II",
    str_detect(tumor_stage, "(^iii$)|(\\si{3}[abc]?$)|(3)") ~ "Stage III",
    str_detect(tumor_stage, "(^iv$)|(\\siv[abc]?$)|(4)") ~ "Stage IV",
    T ~ "Unknown"
  )) %>%
  filter(tumor_stage != "Unknown") %>%
  mutate(tumor_stage = case_when(
    tumor_stage %in% c("Stage I", "Stage II") ~ "Early",
    tumor_stage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage)

if (!dir.exists("~/storage/data/TCGA/DEA/EarlyLateVsNormal")) {
  dir.create("~/storage/data/TCGA/DEA/EarlyLateVsNormal")
}

for (proj in projects) {
  print(str_glue("Working on {proj}:"))
  proj.annot <- annot %>%
    filter(Project == proj) %>%
    select(Barcode, Stage = TumorStage) %>%
    arrange(Stage) %>%
    mutate(Stage = factor(Stage),
           Stage = relevel(Stage, ref = "Normal"))

  design <- model.matrix(~ 0 + Stage, data = proj.annot)
  rownames(design) <- rownames(proj.annot)

  counts <- data.table::fread(str_glue(
    "~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
    column_to_rownames("Ensembl") %>%
    select(one_of(proj.annot$Barcode)) %>%
    head(-5)  # remove  __no_feature, __ambiguous,  __too_low_aQual, __not_aligned, and __alignment_not_unique
  print(str_glue("  Data loaded."))

  dge <- DGEList(counts = counts)
  # Filter out rows (genes) with very low or all-zero counts
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = F]
  print(str_glue("  Filtering finished."))

  # Use the TMM method to scale normalize RNA-Seq counts
  dge <- calcNormFactors(dge)
  v <- voom(dge, design = design, plot = F)
  print(str_glue("  TMM scale normalization finished."))

  # limma pipeline for differential expression
  fit <- lmFit(v, design)

  print(str_glue("  Comparing Early vs. Normal..."))
  contr <- makeContrasts(StageEarly - StageNormal, levels = colnames(coef(fit)))
  fit.early <- contrasts.fit(fit, contr)
  fit.early <- eBayes(fit.early)
  res.early <- topTable(fit.early, number = Inf, sort.by = "P")

  print(str_glue("  Comparing Late vs. Normal..."))
  contr <- makeContrasts(StageLate - StageNormal, levels = colnames(coef(fit)))
  fit.late <- contrasts.fit(fit, contr)
  fit.late <- eBayes(fit.late)
  res.late <- topTable(fit.late, number = Inf, sort.by = "P")

  data.table::fwrite(res.early,
                     file=str_glue("~/storage/data/TCGA/DEA/EarlyLateVsNormal/{proj}_early_vs_normal.csv"),
                     row.names = T)
  data.table::fwrite(res.late,
                     file=str_glue("~/storage/data/TCGA/DEA/EarlyLateVsNormal/{proj}_late_vs_normal.csv"),
                     row.names = T)
}

sessionInfo()
# R version 4.0.0 (2020-04-24)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
#
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
#
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# [1] edgeR_3.30.0    limma_3.44.1    forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0     tibble_3.0.1    ggplot2_3.3.0
# [11] tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.4.6      cellranger_1.1.0  pillar_1.4.4      compiler_4.0.0    dbplyr_1.4.3      tools_4.0.0       jsonlite_1.6.1    lubridate_1.7.8   lifecycle_0.2.0
# [10] nlme_3.1-147      gtable_0.3.0      lattice_0.20-41   pkgconfig_2.0.3   rlang_0.4.6       reprex_0.3.0      cli_2.0.2         DBI_1.1.0         rstudioapi_0.11
# [19] haven_2.2.0       withr_2.2.0       xml2_1.3.2        httr_1.4.1        fs_1.4.1          generics_0.0.2    vctrs_0.3.0       hms_0.5.3         locfit_1.5-9.4
# [28] grid_4.0.0        tidyselect_1.1.0  glue_1.4.1        data.table_1.12.8 R6_2.4.1          fansi_0.4.1       readxl_1.3.1      modelr_0.1.8      magrittr_1.5
# [37] backports_1.1.7   scales_1.1.1      ellipsis_0.3.1    rvest_0.3.5       assertthat_0.2.1  colorspace_1.4-1  stringi_1.4.6     munsell_0.5.0     broom_0.5.6
# [46] crayon_1.3.4
