library(tidyverse)

rm(list = ls())
options(max.print = 50)

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")
pws <- read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                 colNames = F)$X3
pws <- pws[2:length(pws)]

pathways <- c(
  fgsea::gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  fgsea::gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

# Load mutation data -----
mut.dir <- "~/storage/data/archive/muscle/mutation_filtered_by_frequency/combined-stage"
mutations <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- data.table::fread(str_glue("{mut.dir}/{proj}_{stage}.csv")) %>%
      select(-V1) %>%
      mutate(TumorStage = stage)
    mutations <- bind_rows(mutations, dat)
  }
}
mutations <- mutations %>%
  as_tibble()

# Load mutation GSEA and DEA results -----
mut.gsea.dir <- "~/storage/data/archive/muscle/mutation_freq/GSEA_early_late/GSEA_summary"
mut.gsea <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- read.xlsx(str_glue("{mut.gsea.dir}/{proj}.xlsx"), sheet = stage) %>%
      mutate(Project = proj, TumorStage = stage) %>%
      select(Project, TumorStage, everything())
    mut.gsea <- bind_rows(mut.gsea, dat)
  }
}
mut.gsea <- mut.gsea %>%
  filter(Pathway %in% pws) %>%
  as_tibble()

# 每个pathway 在每个肿瘤期, 每个肿瘤，所有肿瘤分别对应的mutated gene list, 基因个数及突变个数 -----
gsea.sum <- mut.gsea %>%
  select(Project, TumorStage, Pathway,
         nPathwayGenes, nLeadingEdge, LeadingEdge) %>%
  mutate(LeadingEdge = map(LeadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
  unnest(LeadingEdge) %>%
  inner_join(mutations %>% select(Project, TumorStage, Symbol, MutationNum),
             by = c("Project", "TumorStage", "LeadingEdge" = "Symbol")) %>%
  mutate(Project = str_remove(Project, "^TCGA-"))

gsea.sum %>%
  group_by(Project, TumorStage, Pathway) %>%
  top_n(5, MutationNum) %>%
  ungroup() %>%
  unite(LeadingEdge, LeadingEdge, MutationNum, sep = ":") %>%
  group_by(Project, TumorStage, Pathway) %>%
  summarise(LeadingEdge = paste(LeadingEdge, collapse = ", ")) %>%
  ungroup() %>%
  spread(TumorStage, LeadingEdge) %>%
  data.table::fwrite("~/selecting_genes.csv")

pw.sum.stage <- gsea.sum %>%
  unite(LeadingEdge, LeadingEdge, MutationNum, sep = ":") %>%
  group_by_at(vars(-LeadingEdge)) %>%
  summarise(LeadingEdge = paste(LeadingEdge, collapse = ", ")) %>%
  rename(`Cancer Type` = Project)

pw.sum.tumor <- gsea.sum %>%
  group_by(Project, Pathway, nPathwayGenes, LeadingEdge) %>%
  summarise(MutationNum = sum(MutationNum)) %>%
  ungroup() %>%
  group_by(Project, Pathway) %>%
  arrange(Pathway, desc(MutationNum)) %>%
  mutate(nLeadingEdge = n()) %>%
  ungroup() %>%
  unite(LeadingEdge, LeadingEdge, MutationNum, sep = ":") %>%
  group_by(Project, Pathway, nPathwayGenes, nLeadingEdge) %>%
  summarise(LeadingEdge = paste(LeadingEdge, collapse = ", ")) %>%
  rename(`Cancer Type` = Project)

pw.sum.all <- gsea.sum %>%
  group_by(Pathway, nPathwayGenes, LeadingEdge) %>%
  summarise(MutationNum = sum(MutationNum)) %>%
  ungroup() %>%
  group_by(Pathway) %>%
  mutate(nLeadingEdge = n()) %>%
  arrange(Pathway, desc(MutationNum)) %>%
  ungroup() %>%
  unite(LeadingEdge, LeadingEdge, MutationNum, sep = ":") %>%
  group_by(Pathway, nPathwayGenes, nLeadingEdge) %>%
  summarise(LeadingEdge = paste(LeadingEdge, collapse = ", "))

openxlsx::write.xlsx(list(
  `Each cancer stage` = pw.sum.stage,
  `Each cancer type` = pw.sum.tumor,
  `Each pathway` = pw.sum.all
), file = "~/storage/data/archive/muscle/supp_tables/summary_stats/pathway_mutation_nums.xlsx")


# 在每一肿瘤中，哪些pathways 只在早期出现，哪些只在晚期出现，哪些一致都出现；
pw.early.late <- mut.gsea %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  filter(`p-value` <= 0.05) %>%
  select(Project, TumorStage, Pathway, nLeadingEdge) %>%
  spread(TumorStage, nLeadingEdge) %>%
  mutate(Occurance = case_when(
    (!is.na(Early)) & (!is.na(Late)) ~ "Both",
    (is.na(Early)) & (!is.na(Late)) ~ "Late",
    (!is.na(Early)) & (is.na(Late)) ~ "Early",
    T ~ "None"
  ))

pw.early.late %>%
  group_by(Project, Occurance) %>%
  summarise(Pathway = paste(Pathway, collapse = ", ")) %>%
  ungroup() %>%
  spread(Occurance, Pathway) %>%
  rename(`Cancer Type` = Project) %>%
  data.table::fwrite("~/storage/data/archive/muscle/supp_tables/summary_stats/pathway_pval0.05_early_late_both.csv")

# 4) 在不同肿瘤中，有哪些pathways 倾向于在绝大多数的肿瘤早期出现；哪些倾向于在绝大多数的肿瘤晚期出现，及一直处出现在绝大多数肿瘤两期中。
pw.early.late %>%
  count(Pathway, Occurance) %>%
  spread(Occurance, n) %>%
  arrange(desc(Both), desc(Early), desc(Late)) %>%
  data.table::fwrite("~/storage/data/archive/muscle/supp_tables/summary_stats/pathway_pval0.05_occurance_stages.csv")

# 5) 有哪些mutated pathways 只出现在某一种肿瘤中？
pw.early.late %>%
  group_by(Pathway) %>%
  summarise(Project = paste(Project, collapse = ", "), n = n()) %>%
  filter(n == 1) %>%
  arrange(Project) %>%
  select(Pathway, `Cancer Type` = Project) %>%
  data.table::fwrite("~/storage/data/archive/muscle/supp_tables/summary_stats/pathway_pval0.05_only_one_cancer.csv")

# For each cancer type, -----
## Total number of selected genes by mutation
mut.gsea %>%
  filter(`p-value` <= 0.05) %>%
  select(Project, TumorStage, Pathway, PathwayGenes) %>%
  mutate(PathwayGenes = map(PathwayGenes, ~ str_split(.x, ", ")[[1]])) %>%
  unnest(PathwayGenes) %>%
  inner_join(mutations %>% select(Project, TumorStage, Symbol, MutationNum),
             by = c("Project", "TumorStage", "PathwayGenes" = "Symbol")) %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  distinct(Project, PathwayGenes) %>%
  count(Project) %>%
  rename(TotalNumSelectedMutatedGenes = n) %>%
  ## number of unique mutated genes that enrich pathways in early/advanced stage samples
  inner_join(mut.gsea %>%
               filter(`p-value` <= 0.05) %>%
               select(Project, TumorStage, Pathway,
                      nPathwayGenes, nLeadingEdge, LeadingEdge) %>%
               mutate(LeadingEdge = map(LeadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
               unnest(LeadingEdge) %>%
               inner_join(mutations %>% select(Project, TumorStage, Symbol, MutationNum),
                          by = c("Project", "TumorStage", "LeadingEdge" = "Symbol")) %>%
               mutate(Project = str_remove(Project, "^TCGA-")) %>%
               distinct(Project, TumorStage, LeadingEdge) %>%
               count(Project, TumorStage) %>%
               spread(TumorStage, n), by = "Project") %>%
  ## number of unique mutated genes that enrich pathways in all cancer samples
  inner_join(mut.gsea %>%
               filter(`p-value` <= 0.05) %>%
               select(Project, TumorStage, Pathway,
                      nPathwayGenes, nLeadingEdge, LeadingEdge) %>%
               mutate(LeadingEdge = map(LeadingEdge, ~ str_split(.x, ", ")[[1]])) %>%
               unnest(LeadingEdge) %>%
               inner_join(mutations %>% select(Project, TumorStage, Symbol, MutationNum),
                          by = c("Project", "TumorStage", "LeadingEdge" = "Symbol")) %>%
               mutate(Project = str_remove(Project, "^TCGA-")) %>%
               distinct(Project, LeadingEdge) %>%
               count(Project) %>%
               rename(TotalNumMutatedGenesAllSamples = n),
             by = "Project") %>%
  data.table::fwrite("~/storage/data/archive/muscle/supp_tables/summary_stats/Table 5.csv")
