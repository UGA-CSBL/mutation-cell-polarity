library(tidyverse)
library(fgsea)

rm(list = ls())
set.seed(1005)

# Load annotation data ------------------------------
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

projects <- annot %>%
  count(Project, TumorStage) %>%
  filter(n >= 15 | Project == "TCGA-LIHC") %>%
  spread(TumorStage, n) %>%
  select(-Normal) %>%
  filter(rowSums(is.na(.)) == 0) %>%
  pull(Project)

annot <- annot %>%
  filter(Project %in% projects) %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  filter(TumorStage != "Normal")
stages <- sort(unique(annot$TumorStage))

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
mutation.dir <- "~/storage/data/archive/muscle/mutation_filtered_by_frequency/combined-stage"
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
  select(Project, TumorStage, Symbol, n = MutationNum)

# GSEA of mutations -----
gsea.analysis <- function(dat) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pathways,
        stats = stats,
        eps = 0,
        # minSize = 15, maxSize = 500,
        scoreType = "pos"
  ) %>%
    arrange(padj, pval)
}

res.mut <- NULL
for (proj in projects) {
  print(str_glue("Working on {proj}"))

  proj.mutations <- mutations %>%
    filter(Project == proj) %>%
    select(Symbol, MutationNum = n) %>%
    arrange(desc(MutationNum))

  proj.res.mut <- gsea.analysis(proj.mutations) %>%
    mutate(Project = proj) %>%
    select(Project, everything())

  res.mut <- bind_rows(res.mut, proj.res.mut)
}
data.table::fwrite(res.mut, file = "~/storage/data/archive/muscle/mutation_freq/filtered_mutations_all_samples_GSEA.csv")

# All mutations -----
mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
  group_by(Impact) %>%
  summarise(Term = paste(Term, collapse = "|"))
mutations <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/DNASeq/GDC-PANCAN.mutect2_snv.tsv"
) %>%
  filter(Sample_ID %in% annot$SampleID) %>%
  mutate(Impact = case_when(
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "HIGH"]) ~ "HIGH",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODERATE"]) ~ "MODERATE",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODIFIER"]) ~ "MODIFIER",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "LOW"]) ~ "LOW",
    T ~ "OTHER"  # not present
  )) %>%
  filter(Impact != "LOW") %>%
  select(SampleID = Sample_ID, Symbol = gene, Impact)

transform.alias <- function(genes) {
  bind_cols(
    Alias = genes,
    Symbol = limma::alias2SymbolTable(genes)
  )
}

mutated.genes <- unique(mutations$Symbol)
mutated.genes <- transform.alias(mutated.genes) %>%
  filter(!is.na(Symbol))

mutations <- mutations %>%
  inner_join(mutated.genes, by = c("Symbol" = "Alias")) %>%
  select(SampleID, Symbol = Symbol.y, Impact)

res.mut <- NULL
for (proj in projects) {
  print(str_glue("Working on {proj}"))

  proj.mutations <- mutations %>%
    filter(SampleID %in% annot$SampleID[annot$Project == proj]) %>%
    count(Symbol, Impact) %>%
    spread(Impact, n) %>%
    replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
    mutate(MutationNum = HIGH + MODERATE + MODIFIER) %>%
    arrange(desc(MutationNum), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
    select(Symbol, MutationNum)

  proj.res.mut <- gsea.analysis(proj.mutations) %>%
    mutate(Project = proj) %>%
    select(Project, everything())

  res.mut <- bind_rows(res.mut, proj.res.mut)
}
data.table::fwrite(res.mut, file = "~/storage/data/archive/muscle/mutation_freq/all_mutations_all_samples_GSEA.csv")

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
# [1] fgsea_1.14.0    forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0     tibble_3.0.1    ggplot2_3.3.0   tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.4.6         lubridate_1.7.8      lattice_0.20-41      assertthat_0.2.1     digest_0.6.25        R6_2.4.1             cellranger_1.1.0
# [8] backports_1.1.7      reprex_0.3.0         stats4_4.0.0         RSQLite_2.2.0        httr_1.4.1           pillar_1.4.4         rlang_0.4.6
# [15] readxl_1.3.1         rstudioapi_0.11      data.table_1.12.8    blob_1.2.1           S4Vectors_0.26.1     Matrix_1.2-18        BiocParallel_1.22.0
# [22] bit_1.1-15.2         munsell_0.5.0        broom_0.5.6          compiler_4.0.0       modelr_0.1.8         BiocGenerics_0.34.0  pkgconfig_2.0.3
# [29] tidyselect_1.1.0     gridExtra_2.3        IRanges_2.22.2       fansi_0.4.1          crayon_1.3.4         dbplyr_1.4.3         withr_2.2.0
# [36] grid_4.0.0           nlme_3.1-148         jsonlite_1.6.1       gtable_0.3.0         lifecycle_0.2.0      DBI_1.1.0            magrittr_1.5
# [43] scales_1.1.1         zip_2.0.4            cli_2.0.2            stringi_1.4.7        fs_1.4.1             limma_3.44.1         xml2_1.3.2
# [50] ellipsis_0.3.1       generics_0.0.2       vctrs_0.3.0          openxlsx_4.1.5       fastmatch_1.1-0      org.Hs.eg.db_3.11.4  tools_4.0.0
# [57] bit64_0.9-7          Biobase_2.48.0       glue_1.4.1           hms_0.5.3            parallel_4.0.0       AnnotationDbi_1.50.0 colorspace_1.4-1
# [64] rvest_0.3.5          memoise_1.1.0        haven_2.3.0
