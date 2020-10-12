library(tidyverse)
library(ggpubr)
library(fgsea)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
) %>%
  enframe() %>%
  select(Pathway = name, Symbol = value) %>%
  mutate(Pathway = str_remove(Pathway, "^REACTOME_|GO_"))

# Load mutation data ------------------------------
data.dir <- "~/storage/data/archive/muscle/selected_mutations"

mutation.dirs <- c(
  str_glue("{data.dir}/mutation_by_selection/adjpval/1e-2"),
  str_glue("{data.dir}/mutation_by_selection/adjpval/1e-3"),
  str_glue("{data.dir}/mutation_by_selection/adjpval/1e-4"),
  str_glue("{data.dir}/mutation_by_selection/adjpval/1e-5"),
  str_glue("{data.dir}/mutation_by_selection/pval/1e-3"),
  str_glue("{data.dir}/mutation_by_selection/pval/1e-4"),
  # str_glue("{data.dir}/mutation_by_selection/pval/1e-5"),
  str_glue("{data.dir}/mutation_by_selection/unfiltered"),
  str_glue("{data.dir}/mutation_by_selection/unfiltered_with_adjpval"),
  str_glue("{data.dir}/mutation_by_chance")
)
res.dirs <- mutation.dirs %>%
  str_replace(str_glue("^({data.dir})/(.*)$"), "\\1/GSEA_early_late/\\2")

for (i in seq(length(mutation.dirs))) {
  mutation.dir <- mutation.dirs[i]
  res.dir <- res.dirs[i]

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

  for (proj in projects) {
    for (stage in stages) {
      res.type <- res.dir %>%
        str_remove(str_glue("^{data.dir}/GSEA_early_late/mutation_by_")) %>%
        str_replace_all("/", "_")
      plt.dir <- str_glue("{data.dir}/mutation_density/{proj}/{stage}")
      if (!dir.exists(plt.dir)) {dir.create(plt.dir, recursive = T)}

      proj.mutations <- mutations %>%
        filter(Project == proj & TumorStage == stage) %>%
        select(Symbol, MutationNum)

      gsea.res <- openxlsx::read.xlsx(str_glue("{res.dir}/{proj}.xlsx"), sheet = stage)
      gsea.res <- gsea.res %>%
        as_tibble() %>%
        slice_min(pval, n = 5) %>%
        select(Pathway = pathway, leadingEdge) %>%
        mutate(Pathway = str_remove(Pathway, "^REACTOME_|GO_"))
      pw.level <- gsea.res$Pathway
      gsea.res <- gsea.res %>%
        inner_join(pathways, by = "Pathway") %>%
        unnest(Symbol) %>%
        left_join(proj.mutations, by = "Symbol") %>%
        replace_na(list(MutationNum = 0)) %>%
        mutate(Pathway = factor(Pathway, levels = pw.level))

      p <- ggdensity(gsea.res, x = "MutationNum",
                     color = "Pathway", fill = "Pathway", palette = "jco",
                     alpha = 0.7, rug = T,
                     title = str_glue("{proj} {stage}"), xlab = "# mutations") %>%
        facet("Pathway", scales = "free") +
        guides(colour = guide_legend(nrow = 3)) +
        theme(axis.text = element_text(size = 14, colour = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              legend.text = element_text(size = 9))

      ggsave(filename = str_glue("{plt.dir}/{res.type}.png"), p,
             width = 12, height = 10, unit = "in", dpi = "retina")
    }
  }
}
