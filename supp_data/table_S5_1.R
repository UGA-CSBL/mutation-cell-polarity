library(tidyverse)

rm(list = ls())

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-READ", "TCGA-STAD", "TCGA-THCA")

gene.list <- c("ACTA1", "ACTA2", "ACTB", "ACTBL2", "ACTC1", "ACTG1", "ACTG2",
               "ACTL10", "ACTL6A", "ACTL6B", "ACTL7A", "ACTL7B", "ACTL8", "ACTL9",
               "KRT1", "KRT10", "KRT12", "KRT13", "KRT14", "KRT15", "KRT16",
               "KRT17", "KRT18", "KRT19", "KRT2", "KRT20", "KRT222", "KRT23",
               "KRT24", "KRT25", "KRT26", "KRT27", "KRT28", "KRT3", "KRT31",
               "KRT32", "KRT33A", "KRT33B", "KRT34", "KRT35", "KRT36", "KRT37",
               "KRT38", "KRT39", "KRT4", "KRT40", "KRT5", "KRT6A", "KRT6B",
               "KRT6C", "KRT7", "KRT71", "KRT72", "KRT73", "KRT74", "KRT75",
               "KRT76", "KRT77", "KRT78", "KRT79", "KRT8", "KRT80", "KRT81",
               "KRT82", "KRT83", "KRT84", "KRT85", "KRT86", "KRT9", "TUBA1A",
               "TUBA1B", "TUBA1C", "TUBA3C", "TUBA3D", "TUBA3E", "TUBA4A",
               "TUBA4B", "TUBA8", "TUBAL3", "TUBB", "TUBB1", "TUBB2A", "TUBB2B",
               "TUBB3", "TUBB4A", "TUBB4B", "TUBB6", "TUBB8", "TUBD1", "TUBE1",
               "TUBG1", "TUBG2")

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
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage)

id.map <- data.table::fread(
  "~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv"
) %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id) %>%
  filter(Symbol %in% gene.list)

# Calculate average TPM -----
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
  ungroup() %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  spread(TumorStage, AvgTPM) %>%
  select(Project, Symbol, Normal, Early, Late)

projects <- str_remove(projects, "^TCGA-")
res <- vector("list", length(projects))
names(res) <- projects

for (proj in projects) {
  res[[proj]] <- dat.avg %>%
    filter(Project == proj) %>%
    select(-Project)
}

openxlsx::write.xlsx(res, file = "~/storage/data/archive/muscle/supp_tables/Table S5_1.xlsx")
