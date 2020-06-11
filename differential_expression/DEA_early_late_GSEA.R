library(tidyverse)
library(fgsea)

rm(list = ls())

# GSEA results (by cancer type and cancer stage) -----
gsea.dir <- "~/storage/data/archive/muscle/mutation_freq/GSEA_early_late"
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

# DEA results -----
id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Ensembl = ensembl_gene_id, Symbol = external_gene_name)

dea.res <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- data.table::fread(str_glue(
      "~/storage/data/TCGA/DEA/EarlyLateVsNormal/{proj}_{tolower(stage)}_vs_normal.csv")) %>%
      rename(Ensembl = V1) %>%
      mutate(Project = proj,
             Stage = stage,
             Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
      inner_join(id.map, by = "Ensembl") %>%
      select(-Ensembl) %>%
      select(Project, Stage, Symbol, log2FC = logFC, p_adjusted = adj.P.Val)
    dea.res <- bind_rows(dea.res, dat)
  }
}

# GSEA using DEA results -----
pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

gsea.analysis <- function(dat, pathways, score.type = "std") {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pathways,
        stats = stats,
        eps = 0,
        # minSize = 15, maxSize = 500,
        scoreType = score.type
  ) %>%
    arrange(padj, pval)
}

res <- NULL
for (proj in projects) {
  print(str_glue("Working on {proj}"))
  set.seed(1005)
  res.up <- vector("list", length(stages))
  res.down <- vector("list", length(stages))
  names(res.up) <- stages
  names(res.down) <- stages
  for (stage in stages) {
    print(str_glue("---{stage}"))

    stat <- dea.res %>%
      filter(Project == proj & Stage == stage) %>%
      filter(!is.na(p_adjusted) & p_adjusted != 0) %>%
      mutate(Score = -log10(p_adjusted) * sign(log2FC)) %>%
      select(Symbol, Score) %>%
      distinct(Symbol, .keep_all = T)

    stat.up <- stat %>%
      filter(Score > 0) %>%
      arrange(desc(Score))
    stat.down <- stat %>%
      filter(Score < 0) %>%
      arrange(desc(Score))

    res.all <- gsea.analysis(stat, pathways, score.type = "std") %>%
      mutate(Project = proj, Stage = stage) %>%
      select(Project, Stage, everything())
    res <- bind_rows(res, res.all)

    res.up[[stage]] <- gsea.analysis(stat.up, pathways, score.type = "pos")
    res.down[[stage]] <- gsea.analysis(stat.down, pathways, score.type = "neg")
  }
  openxlsx::write.xlsx(res.up, file = str_glue("~/storage/data/archive/muscle/DEA_GSEA_early_late/{proj}_up.xlsx"))
  openxlsx::write.xlsx(res.down, file = str_glue("~/storage/data/archive/muscle/DEA_GSEA_early_late/{proj}_down.xlsx"))
}
data.table::fwrite(res, "~/storage/data/archive/muscle/DEA_GSEA_early_late/all.csv")

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
# [1] tidyselect_1.1.0    haven_2.2.0         lattice_0.20-41     colorspace_1.4-1    vctrs_0.3.0         generics_0.0.2      rlang_0.4.6         pillar_1.4.4
# [9] glue_1.4.1          withr_2.2.0         DBI_1.1.0           BiocParallel_1.22.0 dbplyr_1.4.3        modelr_0.1.8        readxl_1.3.1        lifecycle_0.2.0
# [17] munsell_0.5.0       gtable_0.3.0        cellranger_1.1.0    zip_2.0.4           rvest_0.3.5         parallel_4.0.0      fansi_0.4.1         broom_0.5.6
# [25] Rcpp_1.0.4.6        scales_1.1.1        backports_1.1.7     jsonlite_1.6.1      fs_1.4.1            gridExtra_2.3       fastmatch_1.1-0     hms_0.5.3
# [33] stringi_1.4.6       openxlsx_4.1.5      grid_4.0.0          cli_2.0.2           tools_4.0.0         magrittr_1.5        crayon_1.3.4        pkgconfig_2.0.3
# [41] ellipsis_0.3.1      Matrix_1.2-18       data.table_1.12.8   xml2_1.3.2          reprex_0.3.0        lubridate_1.7.8     assertthat_0.2.1    httr_1.4.1
# [49] rstudioapi_0.11     R6_2.4.1            nlme_3.1-147        compiler_4.0.0
