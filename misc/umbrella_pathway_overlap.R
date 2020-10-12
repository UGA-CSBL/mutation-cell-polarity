library(tidyverse)

rm(list = ls())
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

pws <- openxlsx::read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                           colNames = T) %>%
  fill(Umbrella) %>%
  filter(Umbrella == "Development pathway of multiple tissue types") %>%
  pull(Name.of.GO.pathways)

pathways <- c(
  fgsea::gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  fgsea::gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
) %>%
  enframe() %>%
  select(Pathway = name, Genes = value) %>%
  filter(Pathway %in% pws)

# gene2pw <- pathways %>%
#   group_by(gene) %>%
#   summarise(Pathways = paste(term, collapse = ", ")) %>%
#   ungroup() %>%
#   mutate(Pathways = map(Pathways, ~ str_split(.x, ", ")[[1]])) %>%
#   rename(Symbol = gene)

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id) %>%
  filter(Symbol %in% unique(pathways$gene))

# Get GSEA results of DE genes -----
dea.gsea <- NULL
for (project in projects) {
  for (stage in stages) {
    df <- openxlsx::read.xlsx(str_glue(
      "~/storage/data/archive/muscle/DEA_GSEA_early_late/{project}_up.xlsx"),
      sheet = stage) %>%
      mutate(Project = project, Stage = stage, PathwayRegulation = "Up")
    dea.gsea <- bind_rows(dea.gsea, df)
    df <- openxlsx::read.xlsx(str_glue(
      "~/storage/data/archive/muscle/DEA_GSEA_early_late/{project}_down.xlsx"),
      sheet = stage) %>%
      mutate(Project = project, Stage = stage, PathwayRegulation = "Down")
    dea.gsea <- bind_rows(dea.gsea, df)
  }
}

dea.gsea <- dea.gsea %>%
  as_tibble() %>%
  relocate(Project, Stage, PathwayRegulation)

# For each pathway in the umbrella pathway, we focus on their leading edges (LEu)
# check what other pathways these genes enrich. We check the intersection of LEu
# with the leading edges of all other enriched pathways (LEo).
# We find the overlap's significance using Fisher's exact test.
test.overlap <- function(A.size, B.size, intersection.size, union.size) {
  # Find number of all unique gene symbols 15539
  # dea.gsea %>%
  #   select(leadingEdge) %>%
  #   mutate(leadingEdge = map(leadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
  #   unnest(leadingEdge) %>%
  #   distinct()

  contingency.table <- matrix(c(
    15539 - union.size,
    B.size - intersection.size,
    A.size - intersection.size,
    intersection.size
  ), ncol = 2)

  rownames(contingency.table) <- c('notB', 'inB')
  colnames(contingency.table) <- c('notA', 'inA')

  res <- fisher.test(contingency.table, alternative = "greater", conf.int = F)
  res$p.value
}

res <- NULL
for (pw in pws) {
  LEu <- dea.gsea %>%
    filter(pathway == pw) %>%
    mutate(leadingEdge = map(leadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
    select(Project, Stage, PathwayRegulation, leadingEdgeUmbrella = leadingEdge)

  LEo <- dea.gsea %>%
    filter(padj <= 0.01) %>%
    mutate(leadingEdge = map(leadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
    select(Project, Stage, PathwayRegulation,
           Pathway = pathway, padj, NES, leadingEdge)

  pw.res <- LEo %>%
    inner_join(LEu, by = c("Project", "Stage", "PathwayRegulation")) %>%
    mutate(
      GeneOverlap = map2(leadingEdge, leadingEdgeUmbrella, ~ intersect(.x, .y)),
      GeneUnion = map2(leadingEdge, leadingEdgeUmbrella, ~ union(.x, .y)),
      OverlapSize = map_int(GeneOverlap, length),
      UnionSize = map_int(GeneUnion, length),
      LeuSize = map_int(leadingEdgeUmbrella, length),
      LeoSize = map_int(leadingEdge, length),
    ) %>%
    mutate(
      p_overlap = pmap_dbl(
        list(A.size = LeoSize, B.size = LeuSize,
             intersection.size = OverlapSize, union.size = UnionSize),
        test.overlap)) %>%
    filter(p_overlap <= 0.01) %>%
    mutate(PathwayA = pw) %>%
    select(Project, Stage, PathwayRegulation, PathwayA, PathwayB = Pathway, padj,
           PathwayASize = LeuSize, PathwayBSize = LeoSize,
           OverlapSize, OverlapPValue = p_overlap,
           leadingEdgeA = leadingEdgeUmbrella, leadingEdgeB = leadingEdge)

  res <- bind_rows(res, pw.res)
}

data.table::fwrite(res, file = "~/storage/data/archive/muscle/umbrella_development_pathway_overlap.csv")

