library(tidyverse)
library(openxlsx)
library(fgsea)

rm(list = ls())
set.seed(1005)
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")

pws <- read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                 colNames = F)$X3
pws <- pws[2:length(pws)]

# Load mutation data -----
data.dir <- "~/storage/data/archive/muscle/selected_mutations"
mutation.dir <- str_glue("{data.dir}/mutation_by_selection/adjpval/1e-3")
stages <- c("Early", "Late")
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
  select(Project, TumorStage, Symbol, n = MutationNum) %>%
  group_by(Project, Symbol) %>%
  summarise(n = sum(n)) %>%
  ungroup()

# Load mutation GSEA and DEA results -----
mut.dir <- "~/storage/data/archive/muscle/selected_mutations/GSEA_all_samples/GSEA_summary"

mut.gsea <- NULL
for (proj in projects) {
    dat <- data.table::fread(str_glue("{mut.dir}/{proj}.csv")) %>%
      mutate(Project = proj) %>%
      select(Project, everything())
    mut.gsea <- bind_rows(mut.gsea, dat)
}

mut.gsea <- mut.gsea %>%
  rename(pval = `p-value`)

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

dea.res <- NULL
for (proj in projects) {
    dat <- data.table::fread(str_glue(
      "~/storage/data/TCGA/DEA/TumorVsNormal/{proj}.csv")) %>%
      rename(Ensembl = V1) %>%
      mutate(Project = proj,
             Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
      inner_join(id.map, by = "Ensembl") %>%
      select(-Ensembl) %>%
      select(Project, Symbol, log2FC = log2FoldChange, p_adjusted = padj)
    dea.res <- bind_rows(dea.res, dat)
}

dea.res <- dea.res %>%
  as_tibble() %>%
  filter(Symbol %in% mutations$Symbol) %>%
  filter(!is.na(p_adjusted))

# Up-/Down-regulation of pathways with p-value <= 0.2 -----
pw.mut.count <- mut.gsea %>%
  filter(!is.na(pval)) %>%
  filter(pval <= 0.2) %>%
  filter(Pathway %in% pws) %>%
  count(Project)

pw.mut.count[is.na(pw.mut.count)] <- 0

pw.fc <- mut.gsea %>%
  as_tibble() %>%
  filter(Pathway %in% pws) %>%
  filter(!is.na(pval)) %>%
  filter(pval <= 0.2) %>%
  select(Project, Pathway, PathwayGenes) %>%
  mutate(PathwayGenes = map(PathwayGenes, ~ str_split(.x, "\\|")[[1]])) %>%
  unnest(PathwayGenes) %>%
  rename(Symbol = PathwayGenes) %>%
  # Only look for genes that have mutations
  inner_join(mutations, by = c("Project", "Symbol")) %>%
  distinct(Project, Pathway, Symbol, .keep_all = T) %>%
  # Calculate pathway "fold change"
  inner_join(dea.res, by = c("Project", "Symbol")) %>%
  group_by(Project, Pathway) %>%
  mutate(p_adjusted = -log10(p_adjusted)) %>%
  mutate(Score = (p_adjusted * log2FC) / sum(p_adjusted)) %>%
  ungroup() %>%
  filter(abs(Score) >= 0.2) %>%
  group_by(Project, Pathway) %>%
  summarise(PathwayRegulation = sum(Score, na.rm = T)) %>%
  ungroup()

pw.fc %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  mutate(Pathway = factor(Pathway, levels = pws)) %>%
  spread(Project, PathwayRegulation) %>%
  write.xlsx("~/storage/data/archive/muscle/supp_tables/Table_2_all_samples.xlsx")
