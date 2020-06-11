library(tidyverse)
library(fgsea)

rm(list = ls())
set.seed(1005)
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
two.stages <- c("Early", "Late")
four.stages <- paste("Stage", c("I", "II", "III", "IV"), sep = " ")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Ensembl = ensembl_gene_id, Symbol = external_gene_name)

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)
# pws <- openxlsx::read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
#                            colNames = F)$X3
# pws <- pws[2:length(pws)]
# pathways <- pathways[pws]

# DEA GSEA results -----
# dea.gsea <- NULL
# for (proj in projects) {
#   for (stage in four.stages) {
#     dat <- openxlsx::read.xlsx(str_glue("~/storage/data/archive/muscle/DEA_GSEA/{proj}_up.xlsx"), sheet = stage) %>%
#       mutate(Project = proj, Stage = stage, Regulation = "Up")
#     dea.gsea <- bind_rows(dea.gsea, dat)
#
#     dat <- openxlsx::read.xlsx(str_glue("~/storage/data/archive/muscle/DEA_GSEA/{proj}_down.xlsx"), sheet = stage) %>%
#       mutate(Project = proj, Stage = stage, Regulation = "Down")
#     dea.gsea <- bind_rows(dea.gsea, dat)
#   }
# }

# # DEA results -----
dea.four.res <- data.table::fread("~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv") %>%
  as_tibble() %>%
  filter(project %in% projects) %>%
  filter(Symbol %in% unique(unlist(pathways))) %>%
  mutate(comparison = str_replace(comparison, "([^_]+)_vs_N", "Stage \\1")) %>%
  select(Project = project, TumorStage = comparison, Symbol,
         log2FC = log2_fold_change, p_adjusted)

dea.two.res <- NULL
for (proj in projects) {
  for (stage in two.stages) {
    dat <- data.table::fread(str_glue(
      "~/storage/data/TCGA/DEA/EarlyLateVsNormal/{proj}_{tolower(stage)}_vs_normal.csv")) %>%
      rename(Ensembl = V1) %>%
      mutate(Project = proj,
             TumorStage = stage,
             Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
      inner_join(id.map, by = "Ensembl") %>%
      select(-Ensembl) %>%
      select(Project, TumorStage, Symbol, log2FC = logFC, p_adjusted = adj.P.Val)
    dea.two.res <- bind_rows(dea.two.res, dat)
  }
}

# Pathway FC -----
pw.fc <- pathways %>%
  enframe() %>%
  select(Pathway = name, Symbol = value) %>%
  unnest(Symbol) %>%
  # Calculate pathway "fold change"
  inner_join(dea.four.res, by = "Symbol") %>%
  mutate(p_adjusted = -log10(p_adjusted)) %>%
  group_by(Project, TumorStage, Pathway) %>%
  mutate(Score = (p_adjusted * log2FC) / sum(p_adjusted)) %>%
  ungroup() %>%
  filter(abs(Score) >= 0.2)

pw.fc.score <- pw.fc %>%
  group_by(Project, TumorStage, Pathway) %>%
  summarise(PathwayRegulation = sum(Score, na.rm = T)) %>%
  ungroup() %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  pivot_wider(names_from = TumorStage, values_from = PathwayRegulation)

pw.fc.score[is.na(pw.fc.score)] <- 0

# ans <- tibble(Cutoff = seq(0, 0.2, 0.01)) %>%
#   mutate(
#     Data = map(Cutoff, ~ {
#       pw.fc.score <- pw.fc %>%
#         filter(abs(Score) >= .x) %>%
#         group_by(Project, TumorStage, Pathway) %>%
#         summarise(PathwayRegulation = sum(Score, na.rm = T)) %>%
#         ungroup() %>%
#         mutate(Project = str_remove(Project, "^TCGA-")) %>%
#         pivot_wider(names_from = TumorStage, values_from = PathwayRegulation)
#
#       pw.fc.score[is.na(pw.fc.score)] <- 0
#       list(
#         `Stage I only` = pw.fc.score %>%
#           filter(`Stage I` > 0 & `Stage II` <= 0 & `Stage III` <= 0 & `Stage IV` <= 0) %>%
#           count(Project) %>%
#           arrange(desc(n)) %>%
#           unite(Project, Project, n, sep = ":") %>%
#           pull(Project) %>%
#           paste(collapse = ", "),
#         `Stage IV only` = pw.fc.score %>%
#           filter(`Stage I` <= 0 & `Stage II` <= 0 & `Stage III` <= 0 & `Stage IV` > 0) %>%
#           count(Project) %>%
#           arrange(desc(n)) %>%
#           unite(Project, Project, n, sep = ":") %>%
#           pull(Project) %>%
#           paste(collapse = ", "),
#         `Early up` = pw.fc.score %>%
#           filter(`Stage I` <= 0 & `Stage II` <= 0 & `Stage III` <= 0 & `Stage IV` > 0) %>%
#           count(Project) %>%
#           arrange(desc(n)) %>%
#           unite(Project, Project, n, sep = ":") %>%
#           pull(Project) %>%
#           paste(collapse = ", "),
#         `Late up` = pw.fc.score %>%
#           filter(`Stage I` <= 0 & `Stage II` <= 0 & `Stage III` > 0 & `Stage IV` > 0) %>%
#           count(Project) %>%
#           arrange(desc(n)) %>%
#           unite(Project, Project, n, sep = ":") %>%
#           pull(Project) %>%
#           paste(collapse = ", ")
#       )
#     })
#   )
#
# ans %>%
#   mutate(Data = map(Data, ~ tibble(Type = names(.x), Count = unlist(unname(.x))))) %>%
#   unnest(Data) %>%
#   spread(Type, Count) %>%
#   select(`Score cutoff` = Cutoff,
#          `Only stage I up` = `Stage I only`,
#          `Only stage IV up` = `Stage IV only`,
#          `Only Early stages up` = `Early up`,
#          `Only Advanced stages up` = `Late up`) %>%
#   data.table::fwrite("~/different_score_cutoffs_pathway_count.csv")
#
pw.fc.score %>%
  filter(`Stage I` < 0 & `Stage II` < 0 & `Stage III` < 0) %>%
  filter(`Stage IV` > 0) %>%
  arrange(Project, desc(`Stage IV`)) %>%
data.table::fwrite("~/storage/data/archive/muscle/supp_tables/up_in_stageIV.csv")

pw.fc.score %>%
  filter(`Stage I` > 0) %>%
  filter(`Stage II` < 0 & `Stage III` < 0 & `Stage IV` < 0) %>%
  arrange(Project, desc(`Stage IV`)) %>%
  data.table::fwrite("~/storage/data/archive/muscle/supp_tables/up_in_stageI.csv")
