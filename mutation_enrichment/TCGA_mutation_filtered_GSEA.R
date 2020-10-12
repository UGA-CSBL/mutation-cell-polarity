library(tidyverse)
library(fgsea)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

clinical <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/GDC-PANCAN.basic_phenotype.tsv"
) %>%
  filter(str_starts(project_id, "TCGA-")) %>%
  filter(
    sample_type %in% c("Solid Tissue Normal", "Primary Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood")
  )

fpkm.annot <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/fpkm_annot.csv"
) %>%
  mutate(tumor_stage = case_when(
    sample_type == "Solid Tissue Normal" ~ "Normal",
    str_detect(tumor_stage, "(^i$)|(\\si[abc]?$)|(1)") ~ "Stage I",
    str_detect(tumor_stage, "(^ii$)|(\\si{2}[abc]?$)|(2)") ~ "Stage II",
    str_detect(tumor_stage, "(^iii$)|(\\si{3}[abc]?$)|(3)") ~ "Stage III",
    str_detect(tumor_stage, "(^iv$)|(\\siv[abc]?$)|(4)") ~ "Stage IV",
    T ~ "Unknown"
  )) %>%
  filter(tumor_stage != "Unknown") %>%
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage) %>%
  mutate(SampleID = str_sub(Barcode, 1, 16))

annot <- clinical %>%
  select(Project = project_id, SampleID = sample) %>%
  inner_join(fpkm.annot, by = c("Project", "SampleID")) %>%
  distinct(SampleID, .keep_all = T)

## Stages I-IV to early & late
annot <- annot %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  filter(TumorStage != "Normal") %>%
  filter(Project %in% projects)

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
# mutation.dir <- "~/storage/data/archive/muscle/mutation_filtered_by_frequency/combined-stage"
# mutations <- NULL
# for (proj in projects) {
#   for (stage in stages) {
#     dat <- data.table::fread(str_glue("{mutation.dir}/{proj}_{stage}.csv")) %>%
#       select(-V1) %>%
#       mutate(TumorStage = stage)
#     mutations <- bind_rows(mutations, dat)
#   }
# }
data.dir <- "~/storage/data/archive/muscle/selected_mutations"
mutation.dir <- str_glue("{data.dir}/mutation_by_selection/adjpval/1e-3")
mutations <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- data.table::fread(str_glue("{mutation.dir}/{proj}_{stage}.csv")) %>%
      select(-V1) %>%
      mutate(TumorStage = stage)
    mutations <- bind_rows(mutations, dat)
  }
}

mutations <- mutations %>%
  select(Project, TumorStage, Symbol, MutationNum)

# GSEA of mutations -----
gsea.analysis <- function(dat, pws) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pws, stats = stats, eps = 0, scoreType = "pos") %>%
    arrange(padj, pval)
}

if (!dir.exists(str_glue("~/{data.dir}/GSEA_early_late"))) {
  dir.create(str_glue("{data.dir}/GSEA_early_late"), recursive = T)
}

for (project in projects) {
  print(str_glue("Working on {project}"))
  res.mut <- vector("list", length = length(stages))
  names(res.mut) <- stages

  for (stage in stages) {
    set.seed(1005)
    proj.annot <- annot %>%
      filter(Project == project & TumorStage == stage)

    proj.mutations <- mutations %>%
      filter(Project == project & TumorStage == stage) %>%
      select(Symbol, MutationNum) %>%
      arrange(desc(MutationNum))

    res.mut[[stage]] <- gsea.analysis(proj.mutations, pws = pathways)
  }
  openxlsx::write.xlsx(res.mut, file = str_glue(
    "{data.dir}/GSEA_early_late/{project}.xlsx"))
}

# Highly informative pathways -----
# pw.annot <- data.table::fread(
#   "http://supfam.org/SUPERFAMILY/Domain2GO/SDFO.all.txt",
#   skip = 1)
# go.pws <- pw.annot %>%
#   filter(GO_subontology == "biological_process" &
#            str_starts(SDFO_level, "Highly")) %>%
#   pull(GO_name) %>%
#   str_to_upper() %>%
#   str_replace_all("[\\s-]+", "_")
# go.pws <- paste0("GO_", go.pws)
# reactome.pws <- names(pathways)[str_starts(names(pathways), "REACTOME_")]
# pathways.filtered <- pathways[names(pathways) %in% c(go.pws, reactome.pws)]
#
# if (!dir.exists(str_glue("{data.dir}/GSEA_early_late_filtered_pws"))) {
#   dir.create(str_glue("{data.dir}/GSEA_early_late_filtered_pws"),
#              recursive = T)
# }
#
# for (project in projects) {
#   print(str_glue("Working on {project}"))
#   res.mut <- vector("list", length = length(stages))
#   names(res.mut) <- stages
#
#   for (stage in stages) {
#     set.seed(1005)
#     proj.annot <- annot %>%
#       filter(Project == project & TumorStage == stage)
#
#     proj.mutations <- mutations %>%
#       filter(Project == project & TumorStage == stage) %>%
#       select(Symbol, MutationNum) %>%
#       arrange(desc(MutationNum))
#
#     res.mut[[stage]] <- gsea.analysis(proj.mutations, pws = pathways)
#   }
#   openxlsx::write.xlsx(res.mut, file = str_glue(
#     "{data.dir}/GSEA_early_late_filtered_pws/{project}.xlsx"))
# }

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
# [1] fgsea_1.15.1    forcats_0.5.0   stringr_1.4.0   dplyr_1.0.0     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0     tibble_3.0.1    ggplot2_3.3.2   tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.0    haven_2.3.1         lattice_0.20-41     colorspace_1.4-1    vctrs_0.3.1         generics_0.0.2      blob_1.2.1          rlang_0.4.6
# [9] pillar_1.4.4        glue_1.4.1          withr_2.2.0         DBI_1.1.0           BiocParallel_1.22.0 dbplyr_1.4.4        modelr_0.1.8        readxl_1.3.1
# [17] lifecycle_0.2.0     munsell_0.5.0       gtable_0.3.0        cellranger_1.1.0    rvest_0.3.5         zip_2.0.4           parallel_4.0.0      fansi_0.4.1
# [25] broom_0.5.6         Rcpp_1.0.4.6        scales_1.1.1        backports_1.1.8     jsonlite_1.6.1      fs_1.4.1            gridExtra_2.3       fastmatch_1.1-0
# [33] hms_0.5.3           stringi_1.4.7       openxlsx_4.1.5      grid_4.0.0          cli_2.0.2           tools_4.0.0         magrittr_1.5        crayon_1.3.4
# [41] pkgconfig_2.0.3     Matrix_1.2-18       ellipsis_0.3.1      data.table_1.12.8   xml2_1.3.2          reprex_0.3.0        lubridate_1.7.9     assertthat_0.2.1
# [49] httr_1.4.1          rstudioapi_0.11     R6_2.4.1            nlme_3.1-148        compiler_4.0.0
