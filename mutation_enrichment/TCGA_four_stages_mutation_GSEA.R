library(tidyverse)
library(fgsea)

rm(list = ls())

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

# projects <- annot %>%
#   count(Project, TumorStage) %>%
#   filter(n >= 15 | Project == "TCGA-LIHC") %>%
#   spread(TumorStage, n) %>%
#   select(-Normal) %>%
#   filter(rowSums(is.na(.)) == 0) %>%
#   pull(Project)
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")

annot <- annot %>%
  filter(TumorStage != "Normal") %>%
  filter(Project %in% projects)
stages <- sort(unique(annot$TumorStage))

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

# Gene length information -----
# GetGeneLength <- function(symbols) {
#   mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                            dataset = "hsapiens_gene_ensembl")
#
#   biomaRt::getBM(
#     attributes = c("external_gene_name", "start_position", "end_position"),
#     filters = "external_gene_name",
#     values = symbols, mart = mart
#   ) %>%
#     mutate(GeneLength = end_position - start_position) %>%
#     select(Symbol = external_gene_name, GeneLength)
# }
#
# gene.lengths <- GetGeneLength(unique(mutations$Symbol))
# gene.lengths <- gene.lengths %>%
#   group_by(Symbol) %>%
#   summarise(GeneLength = round(sum(GeneLength)))

# GSEA of mutations -----
gsea.analysis <- function(dat) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pathways, stats = stats, eps = 0, scoreType = "pos") %>%
    arrange(padj, pval)
}

if (!dir.exists("~/storage/data/archive/muscle/mutation_freq/GSEA_stages")) {
  dir.create("~/storage/data/archive/muscle/mutation_freq/GSEA_stages")
}

for (proj in projects) {
  set.seed(1005)
  print(str_glue("Working on {proj}"))
  res.mut <- vector("list", length = length(stages))
  names(res.mut) <- stages

  for (stage in stages) {
    print(str_glue("----{stage}"))
    proj.annot <- annot %>%
      filter(Project == proj & TumorStage == stage)

    proj.mutations <- mutations %>%
      filter(SampleID %in% proj.annot$SampleID) %>%
      count(Symbol, Impact) %>%
      spread(Impact, n) %>%
      replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
      mutate(MutationNum = HIGH + MODERATE + MODIFIER) %>%
      arrange(desc(MutationNum), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
      select(Symbol, MutationNum)

    res.mut[[stage]] <- gsea.analysis(proj.mutations)
  }
  openxlsx::write.xlsx(res.mut, file = str_glue(
    "~/storage/data/archive/muscle/mutation_freq/GSEA_stages/{proj}.xlsx"))
}

# sessionInfo()
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
# [1] fgsea_1.14.0    forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.4     readr_1.3.1     tidyr_1.0.3     tibble_3.0.1    ggplot2_3.3.0   tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] tidyselect_1.0.0    haven_2.2.0         lattice_0.20-41     colorspace_1.4-1    vctrs_0.2.4         generics_0.0.2      htmltools_0.4.0     blob_1.2.1
# [9] rlang_0.4.6         pillar_1.4.4        glue_1.4.0          withr_2.2.0         DBI_1.1.0           BiocParallel_1.22.0 bit64_0.9-7         dbplyr_1.4.3
# [17] modelr_0.1.7        readxl_1.3.1        lifecycle_0.2.0     munsell_0.5.0       gtable_0.3.0        cellranger_1.1.0    rvest_0.3.5         zip_2.0.4
# [25] memoise_1.1.0       parallel_4.0.0      fansi_0.4.1         broom_0.5.6         Rcpp_1.0.4.6        scales_1.1.0        backports_1.1.6     jsonlite_1.6.1
# [33] bit_1.1-15.2        fs_1.4.1            gridExtra_2.3       fastmatch_1.1-0     digest_0.6.25       hms_0.5.3           stringi_1.4.6       openxlsx_4.1.5
# [41] grid_4.0.0          cli_2.0.2           tools_4.0.0         magrittr_1.5        RSQLite_2.2.0       crayon_1.3.4        pkgconfig_2.0.3     ellipsis_0.3.0
# [49] Matrix_1.2-18       data.table_1.12.8   xml2_1.3.2          reprex_0.3.0        lubridate_1.7.8     assertthat_0.2.1    httr_1.4.1          rstudioapi_0.11
# [57] R6_2.4.1            nlme_3.1-147        compiler_4.0.0
