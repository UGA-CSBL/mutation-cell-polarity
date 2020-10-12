library(tidyverse)
library(fgsea)

rm(list = ls())

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

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

mut.gsea.dir <- "~/storage/data/archive/muscle/selected_mutations/GSEA_early_late/GSEA_summary"
mut.gsea <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- openxlsx::read.xlsx(str_glue("{mut.gsea.dir}/{proj}.xlsx"), sheet = stage) %>%
      mutate(Project = proj, TumorStage = stage) %>%
      select(Project, TumorStage, everything())
    mut.gsea <- bind_rows(mut.gsea, dat)
  }
}
mut.gsea <- mut.gsea %>%
  filter(Pathway %in% pws) %>%
  as_tibble()

res <- mut.gsea %>%
  count(Pathway, Project) %>%
  rename(Type = n) %>%
  count(Pathway, Type) %>%
  pivot_wider(names_from = "Type", values_from = "n") %>%
  replace_na(list(`1` = 0, `2` = 0)) %>%
  rename(
    `# Cancer types where it is enriched in both stages` = `2`,
    `# Cancer types where it is enriched in one of the stages` = `1`
  ) %>%
  mutate(Pathway = factor(Pathway, levels = pws)) %>%
  arrange(Pathway)

data.table::fwrite(res, file = "~/storage/data/archive/muscle/supp_tables/Table S12.csv")
