library(tidyverse)
library(clusterProfiler)

rm(list = ls())
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

pathways <- bind_rows(
  read.gmt("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  read.gmt("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)  # term, gene

gene2pw <- pathways %>%
  group_by(gene) %>%
  summarise(Pathways = paste(term, collapse = ", ")) %>%
  ungroup() %>%
  mutate(Pathways = map(Pathways, ~ str_split(.x, ", ")[[1]])) %>%
  rename(Symbol = gene)

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id) %>%
  filter(Symbol %in% unique(pathways$gene))

# Run GO enrichment for each subgroup -----
go.enrich <- function(genes, pathways, pws) {
  if (length(pathways) == 0) {
    return(NULL)
  }
  set.seed(1005)
  res <- enricher(genes, TERM2GENE = pws[pws$term %in% pathways, ])
  res@result
}


res <- vector("list", length(projects))
names(res) <- projects

for (project in projects) {
  res[[project]] <- openxlsx::read.xlsx(
    "~/storage/data/archive/muscle/umbrella_pathway_gene_summary.xlsx",
    sheet = project)
}

res <- map(res, as_tibble)

for (project in projects) {
  df <- res[[project]] %>%
    pivot_longer(names_to = "GeneGroup",
                 values_to = "Symbols",
                 -c(Pathway, Stage, PathwayExpressionRegulation)) %>%
    mutate(Symbols = map(Symbols, ~ str_split(.x, ", ")[[1]])) %>%
    mutate(RelatedPathways = map(Symbols, ~ gene2pw %>%
                                   filter(Symbol %in% .x) %>%
                                   unnest(Pathways) %>%
                                   distinct(Pathways) %>%
                                   pull(Pathways))) %>%
    head() %>%
    mutate(
      EnrichGO = map2(Symbols, RelatedPathways, ~ go.enrich(.x, .y, pathways))
    )

  save(df, file = str_glue(
    "~/storage/data/archive/muscle/umbrella_pathway/{project}.RData"))
}
