# Table S3: Numbers of mutations in cytoskeletal genes in the two stages across the nine cancer types.
library(tidyverse)

rm(list = ls())

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
# stages <- str_glue("Stage {c('I', 'II', 'III', 'IV')}")
stages <- c("Early", "Late")
gene.list <- data.table::fread("~/storage/data/archive/muscle/supp_tables/table_S5_genes.csv")

gene.list <- gene.list %>%
  gather(GeneType, Symbol) %>%
  filter(Symbol != "")

# mutation.dir <- "~/storage/data/archive/muscle/mutation_filtered_by_frequency"
# length.cutoff <- 5
# mutations <- NULL
# for (proj in projects) {
#   for (stage in stages) {
#     proj.mutations <- data.table::fread(str_glue("{mutation.dir}/{proj}_{stage}_{length.cutoff}.csv"))
#     mutations <- bind_rows(mutations, proj.mutations)
#   }
# }

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

res <- mutations %>%
  select(Project, TumorStage, Symbol, MutationNum) %>%
  unite(Project, Project, TumorStage) %>%
  spread(Project, MutationNum) %>%
  right_join(gene.list, by = "Symbol") %>%
  select(GeneType, Symbol, everything()) %>%
  arrange(GeneType, Symbol)

res[is.na(res)] <- 0

# openxlsx::write.xlsx(res, "~/storage/data/archive/muscle/supp_tables/table_S5.xlsx")
openxlsx::write.xlsx(res, "~/storage/data/archive/muscle/supp_tables/Table_S3.xlsx")
