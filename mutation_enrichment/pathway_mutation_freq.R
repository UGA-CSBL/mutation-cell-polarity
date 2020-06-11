library(tidyverse)
library(fgsea)
library(reactable)

rm(list = ls())
set.seed(1005)

# Load annotation data ------------------------------
id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-LIHC",
              "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- str_glue("Stage {c('I', 'II', 'III', 'IV')}")

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
  distinct(SampleID, .keep_all = T)

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
# See the link above for impact rating
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
    T ~ "OTHER"
  )) %>%
  filter(Impact != "LOW") %>%
  select(SampleID = Sample_ID, Symbol = gene, Impact) %>%
  inner_join(annot, by = "SampleID") %>%
  count(Project, TumorStage, Symbol, Impact) %>%
  spread(Impact, n) %>%
  replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
  mutate(n = HIGH + MODERATE + MODIFIER) %>%
  arrange(Project, TumorStage, desc(n), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
  select(Project, TumorStage, Symbol, n) %>%
  mutate(n = as.integer(n))

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
  select(Project, TumorStage, Symbol = Symbol.y, n)

mutation.summary <- mutations %>%
  group_by(Project, TumorStage) %>%
  summarise(nGenes = n(),
            nMutations = sum(n)) %>%
  inner_join(count(annot, Project, TumorStage),
             by = c("Project", "TumorStage")) %>%
  rename(SampleNum = n)

if (!file.exists("~/storage/data/archive/muscle/mutation_freq/mutations.xlsx")) {
  openxlsx::write.xlsx(list(Mutations = mutations, Summary = mutation.summary),
                       file = "~/storage/data/archive/muscle/mutation_freq/mutations.xlsx")
}

# DEA results -----
if (file.exists("~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")) {
  dea.res <- data.table::fread("~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")
} else {
  con <- DBI::dbConnect(RSQLite::SQLite(), "~/CSBL_shared/RNASeq/TCGA/DEA/results.db")
  dea.res <- tbl(con, "stages_vs_normal") %>%
    filter(project %in% projects) %>%
    filter(comparison %in% c("I_vs_N", "II_vs_N", "III_vs_N", "IV_vs_N")) %>%
    collect() %>%
    mutate(ensembl = str_sub(ensembl, 1, 15)) %>%  # remove Ensembl version
    inner_join(id.map, by = c("ensembl" = "Ensembl")) %>%
    select(-ensembl) %>%
    select(Symbol, everything()) %>%
    filter(!is.na(p_adjusted) & p_adjusted != 0)
  DBI::dbDisconnect(con)

  for (stage in str_remove(stages, "^Stage ")) {
    lihc.dea <- data.table::fread(str_glue("~/storage/data/TCGA/DEA/csv/TCGA-LIHC_{stage}_vs_N.csv")) %>%
      rename(ensembl = V1) %>%
      mutate(project = "TCGA-LIHC", comparison = paste0(stage, "_vs_N"),
             ensembl = str_sub(ensembl, 1, 15)) %>%  # remove Ensembl version
      inner_join(id.map, by = c("ensembl" = "Ensembl")) %>%
      select(-ensembl) %>%
      select(Symbol, project, comparison, base_mean = baseMean,
             log2_fold_change = log2FoldChange, p_value = pvalue, p_adjusted = padj) %>%
      filter(!is.na(p_adjusted) & p_adjusted != 0)
    dea.res <- bind_rows(dea.res, lihc.dea)
  }

  data.table::fwrite(dea.res, file = "~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")
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

for (proj in projects) {
  proj.res <- vector("list", length(stages))
  names(proj.res) <- stages
  print(str_glue("Working on project {proj}"))
  for (stage in stages) {
    print(str_glue("---{stage}"))
    gsea.res <- openxlsx::read.xlsx(str_glue(
      "~/storage/data/archive/muscle/mutation_freq/GSEA_stages/{proj}.xlsx"
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
      filter(project == proj &
               comparison == paste0(str_remove(stage, "^Stage "), "_vs_N")) %>%
      select(Symbol, log2FC = log2_fold_change, p_adjusted)

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
      htmltools::save_html(file = str_glue("~/storage/data/archive/muscle/mutation_freq/GSEA/{proj}_{stage}.html"))
  }
  openxlsx::write.xlsx(
    proj.res,
    file = str_glue("~/storage/data/archive/muscle/mutation_freq/GSEA/{proj}.xlsx")
  )
}

# Enrichment summary -----
# res <- data.frame(
#   Project = rep(projects, each = 4),
#   TumorStage = rep(stages, times = length(projects)),
#   NumberOfEnrichedGenes = rep(NA, length(projects) * length(stages))
# )
# for (proj in projects) {
#   for (stage in stages) {
#     dat <- openxlsx::read.xlsx(str_glue(
#       "~/storage/data/archive/muscle/mutation_freq/GSEA/{proj}.xlsx"),
#       sheet = stage)
#     if (is.null(dat)) {
#       next()
#     }
#     dat <- paste(dat$LeadingEdge, collapse = ", ")
#     dat <- length(unique(str_split(dat, ", ")[[1]]))
#     res$NumberOfEnrichedGenes[res$Project == proj & res$TumorStage == stage] <- dat
#   }
# }
#
# mutation.summary %>%
#   select(Project, TumorStage, nGenes) %>%
#   inner_join(res, by = c("Project", "TumorStage")) %>%
#   rename(TotalNumberOfGenes = nGenes) %>%
#   data.table::fwrite("~/storage/data/archive/muscle/mutation_freq/gene_number_summary.csv")

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
# [1] reactable_0.1.0.1 fgsea_1.14.0      forcats_0.5.0     stringr_1.4.0     dplyr_0.8.5       purrr_0.3.4       readr_1.3.1       tidyr_1.0.3       tibble_3.0.1
# [10] ggplot2_3.3.0     tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] nlme_3.1-147         fs_1.4.1             lubridate_1.7.8      bit64_0.9-7          progress_1.2.2       httr_1.4.1           tools_4.0.0          backports_1.1.6
# [9] utf8_1.1.4           R6_2.4.1             DBI_1.1.0            BiocGenerics_0.34.0  colorspace_1.4-1     withr_2.2.0          tidyselect_1.0.0     gridExtra_2.3
# [17] prettyunits_1.1.1    curl_4.3             bit_1.1-15.2         compiler_4.0.0       cli_2.0.2            rvest_0.3.5          Biobase_2.48.0       xml2_1.3.2
# [25] scales_1.1.0         askpass_1.1          rappdirs_0.3.1       digest_0.6.25        pkgconfig_2.0.3      htmltools_0.4.0      limma_3.44.1         dbplyr_1.4.3
# [33] htmlwidgets_1.5.1    rlang_0.4.6          readxl_1.3.1         rstudioapi_0.11      RSQLite_2.2.0        generics_0.0.2       jsonlite_1.6.1       BiocParallel_1.22.0
# [41] zip_2.0.4            magrittr_1.5         Matrix_1.2-18        Rcpp_1.0.4.6         munsell_0.5.0        S4Vectors_0.26.0     fansi_0.4.1          lifecycle_0.2.0
# [49] stringi_1.4.6        yaml_2.2.1           org.Hs.eg.db_3.11.1  BiocFileCache_1.12.0 grid_4.0.0           blob_1.2.1           parallel_4.0.0       crayon_1.3.4
# [57] lattice_0.20-41      haven_2.2.0          hms_0.5.3            pillar_1.4.4         biomaRt_2.44.0       stats4_4.0.0         fastmatch_1.1-0      reprex_0.3.0
# [65] XML_3.99-0.3         glue_1.4.0           data.table_1.12.8    modelr_0.1.7         vctrs_0.2.4          cellranger_1.1.0     gtable_0.3.0         openssl_1.4.1
# [73] reactR_0.4.2         assertthat_0.2.1     openxlsx_4.1.5       broom_0.5.6          AnnotationDbi_1.50.0 memoise_1.1.0        IRanges_2.22.1       ellipsis_0.3.0
