library(tidyverse)
library(fgsea)
library(reactable)
library(doSNOW)

rm(list = ls())
set.seed(1005)

# Load annotation data ------------------------------
id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

clinical <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/GDC-PANCAN.basic_phenotype.tsv"
) %>%
  filter(project_id %in% projects) %>%
  filter(sample_type %in% c("Solid Tissue Normal",
                            "Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood"))

fpkm.annot <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/fpkm_annot.csv"
) %>%
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
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage) %>%
  mutate(SampleID = str_sub(Barcode, 1, 16))

annot <- clinical %>%
  select(Project = project_id, SampleID = sample) %>%
  inner_join(fpkm.annot, by = c("Project", "SampleID")) %>%
  distinct(SampleID, .keep_all = T) %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  ))

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

mutation.summary <- mutations %>%
  group_by(Project, TumorStage) %>%
  summarise(nGenes = n(),
            nMutations = sum(n)) %>%
  inner_join(count(annot, Project, TumorStage),
             by = c("Project", "TumorStage")) %>%
  rename(SampleNum = n)

if (!dir.exists("~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary")) {
  dir.create("~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary")
}

if (!file.exists("~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary/mutations.xlsx")) {
  openxlsx::write.xlsx(list(Mutations = mutations, Summary = mutation.summary),
                       file = "~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary/mutations.xlsx")
}

# DEA results -----
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

# Frequency of mutations in interested pathways -----
## Note that genes with only 1 mutation was filtered out during GSEA
## So the `size` doesn't match the actual `nPathwayGenes`
pathways <- enframe(c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)) %>%
  rename(Pathway = name, Symbols = value) %>%
  mutate(nPathwayGenes = map_int(Symbols, length))

cl <- makeCluster(length(projects)+1, outfile="")
registerDoSNOW(cl)
foreach (proj = projects,
         .packages = c("tibble", "tidyr", "purrr", "dplyr", "stringr", "reactable", "openxlsx")) %dopar%
  {
    proj.res <- vector("list", length(stages))
    names(proj.res) <- stages
    print(str_glue("Working on project {proj}"))
    for (stage in stages) {
      print(str_glue("---{stage}"))
      gsea.res <- openxlsx::read.xlsx(str_glue(
        "~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/{proj}.xlsx"
      ), sheet = stage) %>%
        as_tibble() %>%
        filter(pval <= 0.2)

      if (nrow(gsea.res) == 0) {
        next()
      }
      proj.mutations <- mutations %>%
        filter(Project == proj & TumorStage == stage) %>%
        select(Symbol, n)

      proj.dea <- dea.res %>%
        filter(Project == proj & Stage == stage) %>%
        select(Symbol, log2FC, p_adjusted)

      proj.res[[stage]] <- gsea.res %>%
        inner_join(pathways, by = c("pathway" = "Pathway")) %>%
        mutate(leadingEdge = str_split(leadingEdge, ", "),
               nLeadingEdge = map_int(leadingEdge, length),
               PrMutatedGenes = size / nPathwayGenes,
               PrLeadingEdge = nLeadingEdge / nPathwayGenes) %>%
        mutate(
          nPathwayMutations = map_int(Symbols, ~ proj.mutations %>%
                                        filter(Symbol %in% .x) %>%
                                        pull(n) %>% sum()),
          nLeadingEdgeMutations = map_int(leadingEdge, ~ proj.mutations %>%
                                            filter(Symbol %in% .x) %>%
                                            pull(n) %>% sum())
        ) %>%
        select(Pathway = pathway, `p-value` = pval, `Adjusted p-value` = padj,
               nPathwayGenes, nMutatedGenes = size, PrMutatedGenes,
               nLeadingEdge, PrLeadingEdge,
               nPathwayMutations, nLeadingEdgeMutations,
               LeadingEdge = leadingEdge, PathwayGenes = Symbols
        )

      proj.stage.res <- proj.res[[stage]] %>%
        select(-PathwayGenes) %>%
        mutate(LeadingEdge = map(LeadingEdge, ~ tibble(Symbol = .x) %>%
                                   left_join(proj.dea, by = "Symbol") %>%
                                   left_join(proj.mutations, by = "Symbol") %>%
                                   select(Symbol, log2FC, p_adjusted, NumOfMutations = n) %>%
                                   arrange(desc(log2FC)))) %>%
        unnest(LeadingEdge)

      pw.data <- proj.stage.res %>%
        distinct(Pathway, `p-value`, `Adjusted p-value`, nPathwayGenes, nMutatedGenes,
                 PrMutatedGenes, nLeadingEdge, PrLeadingEdge, nPathwayMutations, nLeadingEdgeMutations)

      reactable(pw.data, searchable = T, defaultPageSize = 20,
                details = function(index) {
                  data <- proj.stage.res %>%
                    filter(Pathway == pw.data$Pathway[index]) %>%
                    select(Symbol, log2FC, p_adjusted, NumOfMutations)
                  htmltools::div(style = "padding: 16px",
                                 reactable(data, outlined = TRUE))
                }) %>%
        htmltools::save_html(file = str_glue("~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary/{proj}_{stage}.html"))
    }
    openxlsx::write.xlsx(
      proj.res,
      file = str_glue("~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary/{proj}.xlsx")
    )
  }
stopCluster(cl)

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
# [1] doSNOW_1.0.18     snow_0.4-3        iterators_1.0.12  foreach_1.5.0     reactable_0.1.0.1 fgsea_1.14.0      forcats_0.5.0     stringr_1.4.0     dplyr_0.8.5
# [10] purrr_0.3.4       readr_1.3.1       tidyr_1.1.0       tibble_3.0.1      ggplot2_3.3.0     tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.0    haven_2.2.0         lattice_0.20-41     colorspace_1.4-1    vctrs_0.3.0         generics_0.0.2      htmltools_0.4.0     rlang_0.4.6
# [9] pillar_1.4.4        glue_1.4.1          withr_2.2.0         DBI_1.1.0           BiocParallel_1.22.0 dbplyr_1.4.3        modelr_0.1.8        readxl_1.3.1
# [17] lifecycle_0.2.0     munsell_0.5.0       gtable_0.3.0        cellranger_1.1.0    zip_2.0.4           rvest_0.3.5         htmlwidgets_1.5.1   codetools_0.2-16
# [25] parallel_4.0.0      fansi_0.4.1         broom_0.5.6         Rcpp_1.0.4.6        scales_1.1.1        backports_1.1.7     jsonlite_1.6.1      fs_1.4.1
# [33] gridExtra_2.3       fastmatch_1.1-0     digest_0.6.25       hms_0.5.3           openxlsx_4.1.5      stringi_1.4.6       grid_4.0.0          cli_2.0.2
# [41] tools_4.0.0         magrittr_1.5        crayon_1.3.4        pkgconfig_2.0.3     ellipsis_0.3.1      Matrix_1.2-18       data.table_1.12.8   xml2_1.3.2
# [49] reprex_0.3.0        lubridate_1.7.8     assertthat_0.2.1    httr_1.4.1          rstudioapi_0.11     R6_2.4.1            nlme_3.1-147        compiler_4.0.0
