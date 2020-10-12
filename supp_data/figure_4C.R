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
  select(Project, TumorStage, Symbol, MutationNum)

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

# Cancer type
# Total number of selected genes by mutation
total_mut_genes <- mutations %>%
  inner_join(
    pathways %>%
      enframe() %>%
      filter(name %in% pws) %>%
      unnest(value) %>%
      rename(Pathway = name, Symbol = value),
    by = "Symbol") %>%
  distinct(Project, Symbol) %>%
  count(Project)

# Number of unique mutated genes that enrich pathways in early stage samples
early_gene_num <- mut.gsea %>%
  filter(TumorStage == "Early") %>%
  mutate(LeadingEdge = map(LeadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
  unnest(LeadingEdge) %>%
  group_by(Project) %>%
  summarise(col3 = length(unique(LeadingEdge)))

# Number of unique mutated genes that enrich pathways in advanced stage samples
advanced_gene_num <- mut.gsea %>%
  filter(TumorStage == "Late") %>%
  mutate(LeadingEdge = map(LeadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
  unnest(LeadingEdge) %>%
  group_by(Project) %>%
  summarise(col4 = length(unique(LeadingEdge)))

# Number of unique mutated genes that enrich pathways over cancer samples
all_gene_num <- mut.gsea %>%
  mutate(LeadingEdge = map(LeadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
  unnest(LeadingEdge) %>%
  group_by(Project) %>%
  summarise(col5 = length(unique(LeadingEdge)))

res <- total_mut_genes %>%
  inner_join(early_gene_num, by = "Project") %>%
  inner_join(advanced_gene_num, by = "Project") %>%
  inner_join(all_gene_num, by = "Project") %>%
  mutate(col6 = col5 / n,
         col6 = 100 * round(col6, 4),
         col5 = str_glue("{col5} ({col6}%)")) %>%
  select(-col6) %>%
  rename(`Cancer Type` = "Project",
         `Total # of selected genes by mutation` = "n",
         `# of unique mutated genes that enrich pathways in early stage samples` = "col3",
         `# of unique mutated genes that enrich pathways in advanced stage samples` = "col4",
         `# of unique mutated genes that enrich pathways over cancer samples` = "col5"
  )

data.table::fwrite(res, file = "~/storage/data/archive/muscle/supp_tables/Figure_4C.csv")
