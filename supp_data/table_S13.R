# Table S13: Up-regulated pathways out of the 78 in early vs. advanced stages across nine cancer types.
library(tidyverse)
library(fgsea)

rm(list = ls())
set.seed(1005)
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
two.stages <- c("Early", "Late")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Ensembl = ensembl_gene_id, Symbol = external_gene_name)

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)
pws <- openxlsx::read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                           colNames = F)$X3
pws <- pws[2:length(pws)]
pathways <- pathways[pws]

# DEA GSEA results -----
dea.gsea <- NULL
for (proj in projects) {
  for (stage in two.stages) {
    dat <- openxlsx::read.xlsx(str_glue("~/storage/data/archive/muscle/DEA_GSEA_early_late/{proj}_up.xlsx"), sheet = stage) %>%
      mutate(Project = proj, Stage = stage, Regulation = "Up")
    dea.gsea <- bind_rows(dea.gsea, dat)

    dat <- openxlsx::read.xlsx(str_glue("~/storage/data/archive/muscle/DEA_GSEA_early_late/{proj}_down.xlsx"), sheet = stage) %>%
      mutate(Project = proj, Stage = stage, Regulation = "Down")
    dea.gsea <- bind_rows(dea.gsea, dat)
  }
}

dea.gsea <- dea.gsea %>%
  as_tibble() %>%
  filter(pval <= 0.05) %>%
  filter(pathway %in% pws)

# Pathway FC -----
pathways %>%
  enframe() %>%
  select(Pathway = name) %>%
  left_join(
    dea.gsea %>%
      select(Project, Stage, Regulation, Pathway = pathway) %>%
      mutate(Project = str_remove(Project, "^TCGA-")) %>%
      filter(Regulation == "Up") %>%
      mutate(Stage = case_when(
        Stage == "Late" ~ "Advanced",
        T ~ Stage
      )) %>%
      select(Project, Stage, Pathway),
    by = "Pathway") %>%
  mutate(Regulation = "x") %>%
  unite(Project, Project, Stage, sep = " ") %>%
  pivot_wider(names_from = "Project", values_from = "Regulation") %>%
  select(Pathway, one_of(
    paste(rep(str_remove(projects, "^TCGA-"), each = 2),
          rep(c("Early", "Advanced"), times = length(projects)),
          sep = " ")
  )) %>%
  openxlsx::write.xlsx("~/storage/data/archive/muscle/supp_tables/Table_S13.xlsx")
