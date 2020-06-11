library(tidyverse)
library(fgsea)

rm(list = ls())

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

# Correlation data from Zheng An -----
dat.cor <- NULL
for (project in projects) {
  for (stage in stages) {
    dat <- openxlsx::read.xlsx(
      "~/storage/data/archive/muscle/supp_tables/mutation_cor/mutation_corr.xlsx",
      sheet = str_glue("{project}{stage}")) %>%
      rename(Symbol = 1) %>%
      mutate(Project = project, Stage = stage)
    dat.cor <- bind_rows(dat.cor, dat)
  }
}

dat.cor <- dat.cor %>%
  as_tibble() %>%
  filter(Symbol != "Mutation_Number") %>%
  filter(!is.na(pval)) %>%
  select(Project, Stage, Symbol, RSq = Rsquare, pval) %>%
  filter(RSq > 0)

# GSEA -----
gsea.analysis <- function(dat, pathways, score.type = "pos") {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pathways, stats = stats, eps = 0,
        minSize = 15, scoreType = score.type) %>%
    arrange(padj, pval)
}

res <- vector("list", length(projects) * length(stages))
names(res) <- paste(
  rep(projects, each = length(stages)),
  rep(stages, times = length(projects)),
  sep = "_"
)
for (project in projects) {
  print(str_glue("Working on {project}"))
  set.seed(1005)

  for (stage in stages) {
    print(str_glue("---{stage}"))

    stat <- dat.cor %>%
      filter(Project == project & Stage == stage) %>%
      mutate(Score = -log10(pval) * RSq) %>%
      select(Symbol, Score) %>%
      arrange(desc(Score))

    res[[str_glue("{project}_{stage}")]] <- gsea.analysis(stat, pathways, score.type = "pos") %>%
      select(Pathway = pathway, `p-value` = pval, NES, size, LeadingEdge = leadingEdge)
  }
}
openxlsx::write.xlsx(res, "~/storage/data/archive/muscle/supp_tables/mutation_cor/mutation_cor_GSEA.xlsx")
