library(tidyverse)
library(fgsea)

rm(list = ls())

# GSEA results (by cancer type and cancer stage) -----
gsea.dir <- "~/storage/data/archive/muscle/mutation_freq/GSEA"
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- paste0("Stage ", c("I", "II", "III", "IV"))

# DEA results -----
id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Ensembl = ensembl_gene_id, Symbol = external_gene_name)

if (file.exists("~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")) {
  dea.res <- data.table::fread("~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")
} else {
  con <- DBI::dbConnect(RSQLite::SQLite(), "~/CSBL_shared/RNASeq/TCGA/DEA/results.db")
  dea.res <- tbl(con, "stages_vs_normal") %>%
    filter(project %in% projects) %>%
    filter(comparison %in% c("I_vs_N", "II_vs_N", "III_vs_N", "IV_vs_N")) %>%
    collect() %>%
    mutate(ensembl = str_sub(ensembl, 1, 15)) %>%  # remove Ensembl version
    inner_join(id.map, by = c("ensembl" = "Ensembl")) %>%
    select(-ensembl) %>%
    select(Symbol, everything()) %>%
    filter(!is.na(p_adjusted) & p_adjusted != 0)
  DBI::dbDisconnect(con)

  for (stage in str_remove(stages, "^Stage ")) {
    lihc.dea <- data.table::fread(str_glue("~/storage/data/TCGA/DEA/csv/TCGA-LIHC_{stage}_vs_N.csv")) %>%
      rename(ensembl = V1) %>%
      mutate(project = "TCGA-LIHC", comparison = paste0(stage, "_vs_N"),
             ensembl = str_sub(ensembl, 1, 15)) %>%  # remove Ensembl version
      inner_join(id.map, by = c("ensembl" = "Ensembl")) %>%
      select(-ensembl) %>%
      select(Symbol, project, comparison, base_mean = baseMean,
             log2_fold_change = log2FoldChange, p_value = pvalue, p_adjusted = padj) %>%
      filter(!is.na(p_adjusted) & p_adjusted != 0)
    dea.res <- bind_rows(dea.res, lihc.dea)
  }

  data.table::fwrite(dea.res, file = "~/storage/data/TCGA/DEA/9_cancer_types_staged_DEA.csv")
}

# GSEA using DEA results -----
pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

gsea.analysis <- function(dat, pathways, score.type) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pathways, stats = stats, eps = 0, scoreType = score.type) %>%
    arrange(padj, pval)
}

for (proj in projects) {
  set.seed(1005)
  print(str_glue("Working on {proj}"))
  res.up <- vector("list", length(stages))
  res.down <- vector("list", length(stages))
  names(res.up) <- stages
  names(res.down) <- stages
  for (stage in stages) {
    print(str_glue("---{stage}"))
    comparison <- paste0(str_remove(stage, "^Stage "), "_vs_N")

    stat <- dea.res %>%
      filter(project == proj & comparison == comparison) %>%
      filter(!is.na(p_adjusted) & p_adjusted != 0) %>%
      mutate(Score = -log10(p_adjusted) * sign(log2_fold_change)) %>%
      select(Symbol, Score)

    stat.up <- stat %>%
      filter(Score > 0) %>%
      arrange(desc(Score))
    stat.down <- stat %>%
      filter(Score < 0) %>%
      arrange(desc(Score))

    res.up[[stage]] <- gsea.analysis(stat.up, pathways, score.type = "pos")
    res.down[[stage]] <- gsea.analysis(stat.down, pathways, score.type = "neg")
  }
  openxlsx::write.xlsx(res.up, file = str_glue("~/storage/data/archive/muscle/DEA_GSEA/{proj}_up.xlsx"))
  openxlsx::write.xlsx(res.down, file = str_glue("~/storage/data/archive/muscle/DEA_GSEA/{proj}_down.xlsx"))
}
