library(tidyverse)

rm(list = ls())

# Gene list --------------------------------------------------------------------
if (!file.exists("~/storage/data/archive/muscle/myogenesis.csv")) {
muscles <- c("CARDIAC", "SKELETAL", "SMOOTH", "STRIATED",
             "VASCULAR_SMOOTH", "CARDIAC_VASCULAR_SMOOTH")
biological.processes <- c(
  "MUSCLE_CELL_DIFFERENTIATION",
  "MUSCLE_TISSUE_DEVELOPMENT",
  "MUSCLE_ORGAN_DEVELOPMENT",
  "MUSCLE_CELL_PROLIFERATION",
  "MUSCLE_CELL_MIGRATION",
  "MUSCLE_CONTRACTION",
  "MUSCLE_CELL_MIGRATION",
  "MUSCLE_TISSUE_REGENERATION"
)
known.pws <- paste(
  rep(muscles, each = length(biological.processes)),
  rep(biological.processes, times = length(muscles)),
  sep = "_"
)

pathways <- c(
  "REACTOME_MYOGENESIS",
  str_glue("GO_{known.pws}"),
  str_glue("GO_REGULATION_OF_{known.pws}"),
  str_glue("GO_POSITIVE_REGULATION_OF_{known.pws}"),
  str_glue("GO_NEGATIVE_REGULATION_OF_{known.pws}")
)

gene.list <- c(
  fgsea::gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt"),
  fgsea::gmtPathways("~/storage/data/annotation/MSigDB/c5.all.v7.1.symbols.gmt")
) %>%
  enframe() %>%
  rename(Pathway = name, Symbol = value) %>%
  filter(Pathway %in% pathways) %>%
  # Drop REGULATION_* pathways if POSITIVE & NEGATIVE REGULATION are present
  mutate(
    Process = str_replace(Pathway, "^GO_.*OF_(.*)$", "\\1"),
    Regulation = case_when(
      str_detect(Pathway, "POSITIVE") ~ "P",
      str_detect(Pathway, "NEGATIVE") ~ "N",
      str_detect(Pathway, "REGULATION") ~ "R",
      T ~ NA_character_
    )) %>%
  arrange(Process, Regulation) %>%
  group_by(Process) %>%
  mutate(RegulationCount = n()) %>%
  ungroup() %>%
  filter(!(RegulationCount == 3 & Regulation == "R")) %>%
  select(Pathway, Symbol)

data.table::fwrite(gene.list, "~/storage/data/archive/muscle/myogenesis.csv")
} else {
  gene.list <- data.table::fread("~/storage/data/archive/muscle/myogenesis.csv")
}

# Annotation files -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-LIHC")
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
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage)

id.map <- data.table::fread(
  "~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv"
) %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id) %>%
  filter(Symbol %in% gene.list$Symbol)

FPKMtoTPM <- function(x) {
  1e6 * x / sum(x)
}

res <- NULL
for (proj in projects) {
  dat <- data.table::fread(str_glue(
    "~/CSBL_shared/RNASeq/TCGA/FPKM/{proj}.FPKM.csv"
  )) %>%
    mutate_if(is.numeric, FPKMtoTPM) %>%
    mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    select(Symbol, one_of(annot$Barcode[annot$Project == proj])) %>%
    gather(Barcode, TPM, -Symbol) %>%
    inner_join(annot, by = "Barcode") %>%
    group_by(Project, TumorStage, Symbol) %>%
    summarise(MeanTPM = mean(TPM)) %>%
    spread(TumorStage, MeanTPM) %>%
    inner_join(gene.list, by = "Symbol") %>%
    select(Project, Pathway, Symbol, everything())
  res <- bind_rows(res, dat)
}

res <- res %>%
  arrange(Project, Pathway, Symbol)

res %>%
  nest(data = -Project) %>%
  deframe() %>%
  openxlsx::write.xlsx("~/storage/data/archive/muscle/myogenesis_expression.xlsx")

