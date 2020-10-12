library(tidyverse)
library(ggpubr)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c(
  "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC","TCGA-KIRP",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")

sample.groups <- data.table::fread(
  "~/storage/data/archive/muscle/compare_random_groups/out_barcode.txt",
  sep = ":", sep2 = ",", header = F, col.names = c("Label", "SampleID")) %>%
  as_tibble() %>%
  mutate(Label = map(Label, ~ str_split(.x, " ")[[1]])) %>%
  mutate(Project = map_chr(Label, ~ .x[1]),
         Iteration = as.integer(map_chr(Label, ~ .x[3])),
         Group = map_chr(Label, ~ .x[5]),
         SampleID = map(SampleID, ~ str_split(.x, ",")[[1]])) %>%
  select(Project, Iteration, Group, SampleID)

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

clinical <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/GDC-PANCAN.basic_phenotype.tsv"
) %>%
  filter(project_id %in% projects) %>%
  filter(
    sample_type %in% c("Solid Tissue Normal", "Primary Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood")
  )

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
  filter(Project %in% projects) %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  filter(TumorStage != "Normal")

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
# See the link above for impact rating
mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
  group_by(Impact) %>%
  summarise(Term = paste(Term, collapse = "|"))
mutations <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/GDC-PANCAN.mutect2_snv.tsv") %>%
  filter(Sample_ID %in% annot$SampleID) %>%
  mutate(Impact = case_when(
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "HIGH"]) ~ "HIGH",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODERATE"]) ~ "MODERATE",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODIFIER"]) ~ "MODIFIER",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "LOW"]) ~ "LOW",
    T ~ "OTHER"  # not present
  )) %>%
  filter(Impact != "LOW")

# Density plot of all vs. microindels (<4bp) -----
mutations %>%
  inner_join(annot, by = c("Sample_ID" = "SampleID")) %>%
  mutate(MutationLength = end - start + 1) %>%
  select(Project, TumorStage, MutationLength) %>%
  filter(MutationLength < 4) %>%
  ggdensity(x = "MutationLength", color = "TumorStage") %>%
  facet(facet.by = "Project", nrow = 3, scales = "free")

# Correct aliases -----
mutations <- mutations %>%
  mutate(MutationLength = end - start + 1) %>%
  filter(MutationLength < 4) %>%
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
  select(SampleID, Symbol = Symbol.y, Impact) %>%
  inner_join(annot, by = "SampleID") %>%
  select(Project, TumorStage, SampleID, Barcode, Impact, Symbol) %>%
  arrange(Project, TumorStage, Impact, Symbol)

data.table::fwrite(mutations, "~/storage/data/archive/muscle/microindel_list.csv")

# Get gene length information -----
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")

gene.lengths <- biomaRt::getBM(
  attributes = c("external_gene_name", "start_position", "end_position"),
  filters = "external_gene_name",
  values = unique(mutations$Symbol), mart = mart
) %>%
  mutate(GeneLength = end_position - start_position) %>%
  select(Symbol = external_gene_name, GeneLength)

gene.lengths <- gene.lengths %>%
  group_by(Symbol) %>%
  summarise(GeneLength = round(sum(GeneLength))) %>%
  ungroup()

# Get sample sizes -----
sample.num <- annot %>%
  count(Project, TumorStage) %>%
  rename(SampeNum = n)


# Sort by mutation number and impact -----
mutations <- mutations %>%
  count(Project, TumorStage, Symbol, Impact) %>%
  inner_join(gene.lengths, by = "Symbol") %>%
  inner_join(sample.num, by = c("Project", "TumorStage")) %>%
  spread(Impact, n) %>%
  replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
  mutate(MutationNum = HIGH + MODERATE + MODIFIER) %>%
  arrange(Project, TumorStage, desc(MutationNum), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
  select(Project, TumorStage, SampeNum, Symbol, MutationNum, GeneLength)

data.table::fwrite(mutations, "~/storage/data/archive/muscle/microindels.csv")

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
# [1] ggpubr_0.3.0    forcats_0.5.0   stringr_1.4.0   dplyr_1.0.0     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0     tibble_3.0.1    ggplot2_3.3.2   tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] Biobase_2.48.0       httr_1.4.1           bit64_0.9-7          jsonlite_1.6.1       carData_3.0-4        modelr_0.1.8         assertthat_0.2.1     askpass_1.1
# [9] BiocFileCache_1.12.0 stats4_4.0.0         blob_1.2.1           cellranger_1.1.0     progress_1.2.2       pillar_1.4.4         RSQLite_2.2.0        backports_1.1.8
# [17] lattice_0.20-41      glue_1.4.1           limma_3.44.3         digest_0.6.25        ggsignif_0.6.0       rvest_0.3.5          colorspace_1.4-1     XML_3.99-0.3
# [25] pkgconfig_2.0.3      broom_0.5.6          biomaRt_2.44.1       haven_2.3.1          scales_1.1.1         openxlsx_4.1.5       rio_0.5.16           openssl_1.4.1
# [33] generics_0.0.2       farver_2.0.3         IRanges_2.22.2       car_3.0-8            ellipsis_0.3.1       withr_2.2.0          BiocGenerics_0.34.0  cli_2.0.2
# [41] magrittr_1.5         crayon_1.3.4         readxl_1.3.1         memoise_1.1.0        fs_1.4.1             fansi_0.4.1          nlme_3.1-148         rstatix_0.6.0
# [49] xml2_1.3.2           foreign_0.8-80       prettyunits_1.1.1    tools_4.0.0          data.table_1.12.8    hms_0.5.3            org.Hs.eg.db_3.11.4  lifecycle_0.2.0
# [57] S4Vectors_0.26.1     munsell_0.5.0        reprex_0.3.0         zip_2.0.4            AnnotationDbi_1.50.0 compiler_4.0.0       rlang_0.4.6          grid_4.0.0
# [65] rstudioapi_0.11      rappdirs_0.3.1       labeling_0.3         gtable_0.3.0         abind_1.4-5          DBI_1.1.0            curl_4.3             R6_2.4.1
# [73] lubridate_1.7.9      bit_1.1-15.2         stringi_1.4.7        parallel_4.0.0       Rcpp_1.0.4.6         vctrs_0.3.1          dbplyr_1.4.4         tidyselect_1.1.0
