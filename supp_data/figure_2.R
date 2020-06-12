library(tidyverse)
library(ggpubr)

rm(list = ls())

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
polarity.genes <- tibble(
  Complex = c(
    rep("PAR", 5), rep("Crumbs", 5), rep("Scribble", 7),
    rep("Planar complex", 11), rep("Stabilizer", 2), rep("Destabilizer", 3)
  ),
  Symbol = c(
    "PARD3", "PARD3B", "PARD6B", "PRKCI", "PRKCZ",
    "CRB1", "CRB2", "CRB3", "MPP5", "LIN7A",
    "SCRIB", "MARK2", "LLGL1", "LLGL2", "DLG1", "DLG2", "DLG3",
    "VANGL1", "VANGL2", "PRICKLE1", "PRICKLE2", "CELSR1", "CELSR2", "CELSR3", "ANKRD6", "FZD2", "FZD3", "FZD6",
    "VANGL1", "VANGL2",
    "DVL1", "DVL2", "DVL3"
  )
)

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
  filter(Symbol %in% polarity.genes$Symbol)

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

projects <- str_remove(projects, "^TCGA-")
gene.list <- dat.tpm %>%
  group_by(Project, Symbol) %>%
  summarise(TPM = mean(TPM, na.rm = T)) %>%
  ungroup() %>%
  inner_join(polarity.genes, by = "Symbol") %>%
  group_by(Project, Complex) %>%
  top_n(1, desc(TPM)) %>%
  ungroup() %>%
  select(Project, Complex, Symbol) %>%
  mutate(Project = str_remove(Project, "^TCGA-"))

dat.avg <- dat.tpm %>%
  group_by(Project, TumorStage, Symbol) %>%
  summarise(AvgTPM = mean(TPM, na.rm = T)) %>%
  ungroup() %>%
  mutate(Project = str_remove(Project, "^TCGA-"),
         TumorStage = factor(TumorStage, levels = c("Normal", "Early", "Late"))) %>%
  inner_join(polarity.genes, by = "Symbol")

p <- dat.avg %>%
  inner_join(gene.list, by = c("Project", "Complex", "Symbol")) %>%
  mutate(
    TumorStage = case_when(
      TumorStage == "Early" ~ "Early Stage",
      TumorStage == "Late" ~ "Advanced Stage",
      T ~ "Normal"
    ),
    TumorStage = factor(TumorStage, levels = c("Normal", "Early Stage", "Advanced Stage")),
    Complex = factor(Complex, levels = c("PAR", "Crumbs", "Scribble", "Planar complex", "Stabilizer", "Destabilizer"))
  ) %>%
  ggbarplot(x = "Complex", y = "AvgTPM",
            fill = "TumorStage", position = position_dodge(0.8)) %>%
  facet(facet.by = "Project", nrow = 3, scales = "free") %>%
  ggpar(xlab = "Protein complex",
        ylab = "Average expression (TPM)",
        legend.title = "",
        font.main = c(18, "bold", "black"),
        font.legend = 16, font.x = 16, font.y = 16,
        x.text.angle = 30,
        palette = "jco")

ggsave(str_glue("~/storage/data/archive/muscle/supp_tables/figure_2/Figure 2.tiff"),
       p, width = 9, height = 11, units = "in", dpi = 200)
