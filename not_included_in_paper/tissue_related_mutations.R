library(tidyverse)

rm(list = ls())

mut.dir <- "~/storage/data/archive/muscle/mutation_freq/GSEA_stages"
projects <- list.files(mut.dir)
projects <- str_remove(projects, "\\.xlsx$")
stages <- paste0("Stage ", c("I", "II", "III", "IV"))

tissues <- data.table::fread("~/storage/data/archive/muscle/GTEx_GO_term.csv")

res <- NULL
for (proj in projects) {
  for (stage in stages) {
    dat <- openxlsx::read.xlsx(str_glue("{mut.dir}/{proj}.xlsx"), sheet = stage) %>%
      as_tibble() %>%
      filter(pval <= 0.05)

    pws <- dat$pathway
    dat.tissue <- map(tissues$Regex, ~ pws[str_detect(pws, .x)])
    names(dat.tissue) <- tissues$TissueType
    dat.tissue <- dat.tissue %>%
      enframe() %>%
      rename(Tissue = name, Pathway = value) %>%
      unnest(Pathway) %>%
      count(Tissue) %>%
      mutate(Project = proj, Stage = stage)

    res <- bind_rows(res, dat.tissue)
  }
}

res %>%
  mutate(Stage = str_remove(Stage, "^Stage ")) %>%
  unite(Project, Project, Stage) %>%
  spread(Project, n) %>%
  gather(Project, n, -Tissue) %>%
  mutate(n = case_when(!is.na(n) ~ 1, T ~ 0)) %>%
  group_by(Tissue) %>%
  summarise(PMutated = sum(n) / n()) %>%
  arrange(desc(PMutated)) %>%
  data.table::fwrite("~/storage/data/archive/muscle/mutation_freq/tissue_freq.csv")

# Correlation with inositol marker genes -----
res <- data.table::fread("~/storage/data/archive/muscle/mutation_freq/tissue_freq.csv")
annot <- data.table::fread("~/CSBL_shared/RNASeq/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") %>%
  filter(SMAFRZE == "RNASEQ") %>%
  group_by(SMTS) %>%
  summarise(Barcode = paste(SAMPID, collapse = ", ")) %>%
  ungroup() %>%
  mutate(Barcode = map(Barcode, ~ str_split(.x, ", ")[[1]]))
dat.exp <- data.table::fread("~/CSBL_shared/RNASeq/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2)

annot <- data.table::fread("~/CSBL_shared/RNASeq/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") %>%
  filter(SMAFRZE == "RNASEQ") %>%
  distinct(SMTS, SMTSD)
