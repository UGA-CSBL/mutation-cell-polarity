library(tidyverse)
library(openxlsx)
rm(list = ls())

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

res <- NULL
for (project in projects) {
  data.dir <- "~/storage/data/archive/muscle/selected_mutations/GSEA_early_late/GSEA_summary/"
  for (stage in stages) {
    gsea.res <- openxlsx::read.xlsx(str_glue(
      "{data.dir}/{project}.xlsx"
    ), sheet = stage) %>%
      filter(`p-value` < 0.05) %>%
      mutate(Project = project, TumorStage = stage) %>%
      select(Project, TumorStage, everything())

    if (nrow(gsea.res) == 0) {
      next()
    }
    res <- bind_rows(res, gsea.res)
  }
}

res <- res %>%
  as_tibble() %>%
  select(Project, TumorStage, Pathway, `p-value`, `Adjusted p-value` = `Adjusted.p-value`,
         nPathwayGenes, nMutatedGenes, nPathwayMutations) %>%
  mutate(TumorStage = case_when(TumorStage == "Late" ~ "Advanded",
                                T ~ TumorStage)) %>%
  unite(Project, Project, TumorStage, sep = " ") %>%
  group_by(Project) %>%
  nest() %>%
  deframe()

wb.sheets <- names(res)
highlight.style <- createStyle(fontColour = "#000000", bgFill = "#B59562")
wb <- createWorkbook()

for (st in wb.sheets) {
  addWorksheet(wb, st)
  writeData(wb, st, res[[st]])
  conditionalFormatting(wb, st,
                        cols = 1, rows = 1:nrow(res[[st]]),
                        rule = "$B1 < 0.01", style = highlight.style)
}

saveWorkbook(wb, "~/storage/data/archive/muscle/supp_tables/Table_S21.xlsx", TRUE)
