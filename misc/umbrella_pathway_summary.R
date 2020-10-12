# 对于每个specific pathway在每一肿瘤中, 给出 其上下调的状态，基因总数，
# 上调的基因数，下调的基因数，未变的基因数，被选择的突变基因数，其中有多少在上调
# 的基因中，有多少在下调的基因中，多少在未变的基因中。有了这些数据后，我们需要
# 一张直观的图把这些显示出来。我们需要研究出一个好的方法，直观地将这些信息都表示
# 出来。先以tissue development 的14个 pathways 为例。

library(tidyverse)
library(fgsea)

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
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
) %>%
  enframe() %>%
  select(Pathway = name, Genes = value) %>%
  filter(Pathway %in% pws) %>%
  mutate(PathwaySize = map_int(Genes, length))

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id) %>%
  filter(Symbol %in% unique(unlist(pathways$Genes)))

# Differential expression results -----
dea.res <- NULL
for (project in projects) {
  df <- data.table::fread(str_glue(
    "~/storage/data/TCGA/DEA/EarlyLateVsNormal/{project}_early_vs_normal.csv")) %>%
    select(Ensembl = V1, logFC, padj = adj.P.Val) %>%
    filter(padj <= 0.01) %>%
    mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    select(-Ensembl) %>%
    distinct(Symbol, .keep_all = T) %>%
    mutate(Project = project, Stage = "Early")
  dea.res <- bind_rows(dea.res, df)

  df <- data.table::fread(str_glue(
    "~/storage/data/TCGA/DEA/EarlyLateVsNormal/{project}_late_vs_normal.csv")) %>%
    select(Ensembl = V1, logFC, padj = adj.P.Val) %>%
    filter(padj <= 0.01) %>%
    mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    select(-Ensembl) %>%
    distinct(Symbol, .keep_all = T) %>%
    mutate(Project = project, Stage = "Late")
  dea.res <- bind_rows(dea.res, df)
}

dea.res <- dea.res %>%
  mutate(Regulation = case_when(
    logFC > 0 ~ "Up",
    logFC < 0 ~ "Down",
    T ~ "None"
  )) %>%
  select(Project, Stage, Symbol, Regulation)

# Mutation data -----
data.dir <- "~/storage/data/archive/muscle/selected_mutations"
mutation.dir <- str_glue("{data.dir}/mutation_by_selection/adjpval/1e-3")
mutations <- NULL
for (project in projects) {
  for (stage in stages) {
    dat <- data.table::fread(str_glue("{mutation.dir}/{project}_{stage}.csv")) %>%
      select(-V1) %>%
      mutate(TumorStage = stage)
    mutations <- bind_rows(mutations, dat)
  }
}

mutations <- mutations %>%
  select(Project, Stage = TumorStage, Symbol) %>%
  filter(Symbol %in% unique(unlist(pathways$Genes)))

# Summary results -----
res <- vector("list", length(projects))
names(res) <- projects

for (project in projects) {
  ## Pathway up/down regulation -----
  exp.gsea.up.early <- openxlsx::read.xlsx(str_glue(
    "~/storage/data/archive/muscle/DEA_GSEA_early_late/{project}_up.xlsx"),
    sheet = "Early") %>%
    filter(pathway %in% pws) %>%
    filter(padj <= 0.01) %>%
    pull(pathway)

  exp.gsea.up.late <- openxlsx::read.xlsx(str_glue(
    "~/storage/data/archive/muscle/DEA_GSEA_early_late/{project}_up.xlsx"),
    sheet = "Late") %>%
    filter(pathway %in% pws) %>%
    filter(padj <= 0.01) %>%
    pull(pathway)

  exp.gsea.down.early <- openxlsx::read.xlsx(str_glue(
    "~/storage/data/archive/muscle/DEA_GSEA_early_late/{project}_down.xlsx"),
    sheet = "Early") %>%
    filter(pathway %in% pws) %>%
    filter(padj <= 0.01) %>%
    pull(pathway)

  exp.gsea.down.late <- openxlsx::read.xlsx(str_glue(
    "~/storage/data/archive/muscle/DEA_GSEA_early_late/{project}_down.xlsx"),
    sheet = "Late") %>%
    filter(pathway %in% pws) %>%
    filter(padj <= 0.01) %>%
    pull(pathway)

  exp.gsea.early <- pathways %>%
    select(Pathway) %>%
    mutate(
      Stage = "Early",
      PathwayExpressionRegulation = case_when(
        Pathway %in% exp.gsea.up.early & Pathway %in% exp.gsea.down.early ~ "Up & Down",
        Pathway %in% exp.gsea.up.early ~ "Up",
        Pathway %in% exp.gsea.down.early ~ "Down",
        T ~ "None"
      )
    )

  exp.gsea.late <- pathways %>%
    select(Pathway) %>%
    mutate(
      Stage = "Late",
      PathwayExpressionRegulation = case_when(
        Pathway %in% exp.gsea.up.late & Pathway %in% exp.gsea.down.late ~ "Up & Down",
        Pathway %in% exp.gsea.up.late ~ "Up",
        Pathway %in% exp.gsea.down.late ~ "Down",
        T ~ "None"
      )
    )

  exp.gsea <- bind_rows(exp.gsea.early, exp.gsea.late)

  rm(exp.gsea.early, exp.gsea.late,
     exp.gsea.up.early, exp.gsea.up.late, exp.gsea.down.early, exp.gsea.down.late)


  ## Number of up-regulated genes, down-regulated genes and not DE genes -----
  proj.de.gene <- pathways %>%
    mutate(PathwayGenes = Genes) %>%
    unnest(Genes) %>%
    inner_join(
      dea.res %>%
        filter(Project == project) %>%
        select(-Project),
      by = c("Genes" = "Symbol")) %>%
    nest(DEGenes = c(Genes)) %>%
    mutate(DEGenes = map(DEGenes, ~ .x$Genes)) %>%
    select(-PathwaySize) %>%
    pivot_wider(names_from = Regulation, values_from = DEGenes) %>%
    mutate(NotDE = pmap(list(.all = PathwayGenes, .x = Up, .y = Down),
                        function(.all, .x, .y) .all[!.all %in% .x & !.all %in% .y]))

  ## Selected mutated genes -----
  proj.mut.gene <- pathways %>%
    unnest(Genes) %>%
    inner_join(
      mutations %>%
        filter(Project == project) %>%
        select(-Project),
      by = c("Genes" = "Symbol")) %>%
    nest(MutGenes = c(Genes)) %>%
    mutate(MutGenes = map(MutGenes, ~ .x$Genes)) %>%
    select(-PathwaySize)

  ## Mutated & Up, Mutated & Down, Mutated & not DE -----
  proj.res <- proj.mut.gene %>%
    inner_join(proj.de.gene, by = c("Pathway", "Stage")) %>%
    mutate(
      MutUp = map2(MutGenes, Up, ~ intersect(.x, .y)),
      MutDown = map2(MutGenes, Down, ~ intersect(.x, .y)),
      MutNotDE = map2(MutGenes, NotDE, ~ intersect(.x, .y))
    )

  ## Final summary -----
  res[[project]] <- exp.gsea %>%
    inner_join(proj.res, by = c("Pathway", "Stage")) %>%
    # mutate_at(4:ncol(.), function(x) map_int(x, length)) %>%
    select(Pathway, Stage, PathwayExpressionRegulation,
           PathwayGenes = PathwayGenes,
           UpRegulatedGenes = Up,
           DownRegulatedGenes = Down,
           NotDEGenes = NotDE,
           MutatedGenes = MutGenes,
           MutatedAndUp = MutUp,
           MutatedAndDown = MutDown,
           MutatedAndNotDE = MutNotDE
    )
}

openxlsx::write.xlsx(res, file = "~/storage/data/archive/muscle/umbrella_pathway_gene_summary.xlsx")
