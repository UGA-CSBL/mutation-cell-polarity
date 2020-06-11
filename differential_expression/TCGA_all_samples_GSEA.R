library(tidyverse)
library(fgsea)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

# GSEA ------------------------------
gsea.analysis <- function(dat, pathways, score.type = "std") {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pathways, stats = stats, eps = 0, scoreType = score.type) %>%
    arrange(padj, pval)
}

res <- NULL
res.up <- NULL
res.down <- NULL
for (proj in projects) {
  set.seed(1005)
  print(str_glue("Working on {proj}"))

  proj.dat <- data.table::fread(
    str_glue("~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal/{proj}.csv")) %>%
    rename(Ensembl = V1) %>%
    mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    filter(!is.na(pvalue) & pvalue != 0) %>%
    mutate(Score = -log10(pvalue) * sign(log2FoldChange)) %>%
    select(Symbol, Score) %>%
    distinct(Symbol, .keep_all = T) %>%
    arrange(desc(Score))

  stat.up <- proj.dat %>%
    filter(Score > 0) %>%
    arrange(desc(Score))
  stat.down <- proj.dat %>%
    filter(Score < 0) %>%
    arrange(desc(Score))

  proj.res.all <- gsea.analysis(proj.dat, pathways, score.type = "std") %>%
    mutate(Project = proj) %>%
    select(Project, everything())
  res <- bind_rows(res, proj.res.all)

  proj.res.up <- gsea.analysis(stat.up, pathways, score.type = "pos") %>%
    mutate(Project = proj) %>%
    select(Project, everything())
  res.up <- bind_rows(res.up, proj.res.up)

  proj.res.down <- gsea.analysis(stat.down, pathways, score.type = "neg") %>%
    mutate(Project = proj) %>%
    select(Project, everything())
  res.down <- bind_rows(res.down, proj.res.down)

}
data.table::fwrite(res, file = "~/storage/data/archive/muscle/DEA_GSEA_all/GSEA_all_samples.csv")
data.table::fwrite(res.up, file = "~/storage/data/archive/muscle/DEA_GSEA_all/GSEA_all_samples_up.csv")
data.table::fwrite(res.down, file = "~/storage/data/archive/muscle/DEA_GSEA_all/GSEA_all_samples_down.csv")

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
