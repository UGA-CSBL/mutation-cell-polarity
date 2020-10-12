# Table S8 is what used to be Table 4. Table S9 was Table S7.
# Table S8: The number of enriched pathways and the differential expressions of the mutated genes.
# Table S9: Fold-changes in gene-expressions in 78 selected pathways in two stages across nine cancer types.
library(tidyverse)
library(openxlsx)
library(fgsea)

rm(list = ls())
set.seed(1005)
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
# stages <- str_glue("Stage {c('I', 'II', 'III', 'IV')}")
stages <- c("Early", "Late")

pws <- read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                 colNames = F)$X3
pws <- pws[2:length(pws)]

# Load mutation data -----
data.dir <- "~/storage/data/archive/muscle/selected_mutations"
mutation.dir <- str_glue("{data.dir}/mutation_by_selection/adjpval/1e-3")
mutations <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- data.table::fread(str_glue("{mutation.dir}/{proj}_{stage}.csv")) %>%
      select(-V1) %>%
      mutate(TumorStage = stage)
    mutations <- bind_rows(mutations, dat)
  }
}

mutations <- mutations %>%
  select(Project, TumorStage, Symbol, n = MutationNum)

# Load mutation GSEA and DEA results -----
mut.dir <- "~/storage/data/archive/muscle/selected_mutations/GSEA_early_late/GSEA_summary"

mut.gsea <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- read.xlsx(str_glue("{mut.dir}/{proj}.xlsx"), sheet = stage) %>%
      mutate(Project = proj, TumorStage = stage) %>%
      select(Project, TumorStage, everything())
    mut.gsea <- bind_rows(mut.gsea, dat)
  }
}

mut.gsea <- mut.gsea %>%
  rename(pval = `p-value`)

# dea.res <- data.table::fread("~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")

# dea.res <- dea.res %>%
#   as_tibble() %>%
#   mutate(TumorStage = paste0("Stage ", str_remove(comparison, "_.+$"))) %>%
#   filter(!is.na(p_adjusted)) %>%
#   filter(p_adjusted != 0) %>%
#   # filter(p_adjusted <= 0.01) %>%
# select(Project = project, TumorStage, Symbol, log2FC = log2_fold_change)

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

dea.res <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- data.table::fread(str_glue(
      "~/storage/data/TCGA/DEA/EarlyLateVsNormal/{proj}_{tolower(stage)}_vs_normal.csv")) %>%
      rename(Ensembl = V1) %>%
      mutate(Project = proj,
             Stage = stage,
             Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
      inner_join(id.map, by = "Ensembl") %>%
      select(-Ensembl) %>%
      select(Project, Stage, Symbol, log2FC = logFC, p_adjusted = adj.P.Val)
    dea.res <- bind_rows(dea.res, dat)
  }
}

dea.res <- dea.res %>%
  as_tibble() %>%
  filter(Symbol %in% mutations$Symbol)

# Run GSEA using DEA data for mutated genes -----
# gsea.analysis <- function(dat, pathways) {
#   stats <- dat %>%
#     deframe()
#
#   fgsea(pathways = pathways, stats = stats, eps = 0,
#         minSize = 15, maxSize = 500, scoreType = "std") %>%
#     arrange(padj, pval)
# }
#
# pathways <- c(
#   gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
#   gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
# )
#
# if (!file.exists("~/storage/data/archive/muscle/DEA_GSEA_early_late/all.csv")) {
#
#   dea.gsea <- NULL
#   for (proj in projects) {
#     print(str_glue("Working on {proj}"))
#     for (stage in stages) {
#       print(str_glue("---{stage}"))
#
#       stat <- dea.res %>%
#         filter(Project == proj & TumorStage == stage) %>%
#         distinct(Symbol, Score)
#
#       proj.dea.gsea <- gsea.analysis(stat, pathways) %>%
#         mutate(Project = proj, Stage = stage) %>%
#         select(Project, Stage, everything())
#
#       dea.gsea <- bind_rows(dea.gsea, proj.dea.gsea)
#     }
#   }
#   data.table::fwrite(dea.gsea, "~/storage/data/archive/muscle/DEA_GSEA_early_late/all.csv")
# } else {
#   dea.gsea <- data.table::fread("~/storage/data/archive/muscle/DEA_GSEA_early_late/all.csv")
# }
#
# dea.gsea.reg <- dea.gsea %>%
#   filter(pathway %in% pws) %>%
#   mutate(NES = case_when(
#     NES > 0 ~ "Up",
#     NES < 0 ~ "Down",
#     T ~ "No"
#   )) %>%
#   filter(NES != "No") %>%
#   select(Project, TumorStage = Stage, Pathway = pathway, NES)


# Up-/Down-regulation of pathways with p-value <= 0.2 -----
pw.mut.count <- mut.gsea %>%
  filter(!is.na(pval)) %>%
  filter(pval <= 0.2) %>%
  filter(Pathway %in% pws) %>%
  count(Project, TumorStage) %>%
  spread(TumorStage, n)

pw.mut.count[is.na(pw.mut.count)] <- 0
write.xlsx(pw.mut.count, "~/storage/data/archive/muscle/supp_tables/TableS8_mut.xlsx")

pw.fc <- mut.gsea %>%
  as_tibble() %>%
  filter(Pathway %in% pws) %>%
  filter(!is.na(pval)) %>%
  filter(pval <= 0.2) %>%
  # select(Project, TumorStage, Pathway = pathway, leadingEdge) %>%

  select(Project, TumorStage, Pathway, PathwayGenes) %>%
  mutate(PathwayGenes = map(PathwayGenes, ~ str_split(.x, ", ")[[1]])) %>%
  unnest(PathwayGenes) %>%
  rename(Symbol = PathwayGenes) %>%
  # Only look for genes that have mutations
  inner_join(mutations, by = c("Project", "TumorStage", "Symbol")) %>%
  select(-n) %>%
  # Calculate pathway "fold change"
  inner_join(dea.res, by = c("Project", "TumorStage" = "Stage", "Symbol")) %>%
  group_by(Project, TumorStage, Pathway) %>%
  mutate(p_adjusted = -log10(p_adjusted)) %>%
  mutate(Score = (p_adjusted * log2FC) / sum(p_adjusted)) %>%
  ungroup() %>%
  filter(abs(Score) >= 0.2) %>%
  group_by(Project, TumorStage, Pathway) %>%
  summarise(PathwayRegulation = sum(Score, na.rm = T)) %>%
  ungroup()

pw.fc.count <- pw.fc %>%
  mutate(PathwayRegulation = case_when(
    PathwayRegulation > 0 ~ "Up",
    PathwayRegulation < 0 ~ "Down",
    T ~ "None"
  )) %>%
  count(Project, TumorStage, PathwayRegulation) %>%
  unite(TumorStage, TumorStage, PathwayRegulation) %>%
  spread(TumorStage, n) %>%
  select(Project, Early_Up, Early_Down, Late_Up, Late_Down) %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  rename(`Cancer Type` = Project)

pw.fc.count[is.na(pw.fc.count)] <- 0
write.xlsx(pw.fc.count, "~/storage/data/archive/muscle/supp_tables/TableS8_DE.xlsx")

pw.fc %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  mutate(Pathway = factor(Pathway, levels = pws)) %>%
  unite(TumorStage, Project, TumorStage) %>%
  spread(TumorStage, PathwayRegulation) %>%
  write.xlsx("~/storage/data/archive/muscle/supp_tables/Table_S9.xlsx")

# pw.fc %>%
#   mutate(Project = str_remove(Project, "^TCGA-")) %>%
#   mutate(Pathway = factor(Pathway, levels = pws)) %>%
#   group_by(Project, Pathway) %>%
#   summarise(PathwayRegulation = sum(PathwayRegulation, na.rm = T)) %>%
#   ungroup() %>%
#   spread(Project, PathwayRegulation) %>%
#   write.xlsx("~/storage/data/archive/muscle/supp_tables/Table S9_early_late.xlsx")
