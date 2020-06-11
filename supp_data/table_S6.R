library(tidyverse)

rm(list = ls())
options(max.print = 50)

# Selected genes in each protein complex -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-READ", "TCGA-STAD", "TCGA-THCA")
# stages <- str_glue("Stage {c('I', 'II', 'III', 'IV')}")
stages <- c("Early", "Late")
gene.list <- data.table::fread("~/storage/data/archive/muscle/supp_tables/table_S5_genes.csv")

gene.list <- gene.list %>%
  gather(GeneType, Symbol) %>%
  filter(Symbol != "")

# mutation.dir <- "~/storage/data/archive/muscle/mutation_filtered_by_frequency"
# length.cutoff <- 5
# mutations <- NULL
# for (proj in projects) {
#   for (stage in stages) {
#     proj.mutations <- data.table::fread(str_glue("{mutation.dir}/{proj}_{stage}_{length.cutoff}.csv"))
#     mutations <- rbind(mutations, proj.mutations)
#   }
# }

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

res <- mutations %>%
  select(Project, TumorStage, Symbol, MutationNum) %>%
  filter(Symbol %in% gene.list$Symbol) %>%
  group_by(Project, Symbol) %>%
  summarise(MeanMutationNum = sum(MutationNum, na.rm = T) / 2) %>%
  ungroup() %>%
  inner_join(gene.list, by = "Symbol") %>%
  group_by(Project, GeneType) %>%
  top_n(n = 10, wt = MeanMutationNum) %>%
  select(Project, GeneType, Symbol, MeanMutationNum) %>%
  arrange(Project, GeneType, desc(MeanMutationNum))

# Average expression in TPM -----
annot <- data.table::fread(
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
  mutate(tumor_stage = case_when(
    tumor_stage %in% c("Stage I", "Stage II") ~ "Early",
    tumor_stage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal")) %>%
  filter(tumor_stage != "Normal") %>%
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage)

id.map <- data.table::fread(
  "~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv"
) %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id) %>%
  filter(Symbol %in% gene.list$Symbol)

FPKM2TPM <- function(x) {
  1e6 * x / sum(x)
}

dat.tpm <- NULL
for (proj in projects) {
  dat <- data.table::fread(
    str_glue("~/CSBL_shared/RNASeq/TCGA/FPKM/{proj}.FPKM.csv")
  ) %>%
    select(Ensembl, one_of(annot$Barcode[annot$Project == proj])) %>%
    mutate_if(is.numeric, FPKM2TPM) %>%
    mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    select(-Ensembl) %>%
    gather(Barcode, TPM, -Symbol) %>%
    inner_join(annot, by = "Barcode") %>%
    select(Project, TumorStage, Symbol, Barcode, TPM)
  dat.tpm <- bind_rows(dat.tpm, dat)
}

dat.avg <- dat.tpm %>%
  group_by(Project, TumorStage, Symbol) %>%
  summarise(AvgTPM = mean(TPM, na.rm = T)) %>%
  spread(Project, AvgTPM) %>%
  gather(Project, TPM, -c(TumorStage, Symbol)) %>%
  spread(TumorStage, TPM)

res %>%
  inner_join(dat.avg, by = c("Project", "Symbol")) %>%
  select(-MeanMutationNum) %>%
  gather(TumorStage, TPM, -c(Project, GeneType, Symbol)) %>%
  unite(Project, Project, TumorStage) %>%
  spread(Project, TPM) %>%
  arrange(GeneType, Symbol) %>%
  openxlsx::write.xlsx("~/storage/data/archive/muscle/supp_tables/table_S6_early_late.xlsx")
