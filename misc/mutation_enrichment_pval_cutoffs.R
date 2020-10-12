library(tidyverse)
library(fgsea)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

clinical <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/GDC-PANCAN.basic_phenotype.tsv"
) %>%
  filter(project_id %in% projects) %>%
  filter(
    sample_type %in% c("Solid Tissue Normal", "Primary Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood")
  )

fpkm.annot <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/fpkm_annot.csv"
) %>%
  filter(project %in% projects) %>%
  mutate(tumor_stage = case_when(
    sample_type == "Solid Tissue Normal" ~ "Normal",
    str_detect(tumor_stage, "(^i$)|(\\si[abc]?$)|(1)") ~ "Stage I",
    str_detect(tumor_stage, "(^ii$)|(\\si{2}[abc]?$)|(2)") ~ "Stage II",
    str_detect(tumor_stage, "(^iii$)|(\\si{3}[abc]?$)|(3)") ~ "Stage III",
    str_detect(tumor_stage, "(^iv$)|(\\siv[abc]?$)|(4)") ~ "Stage IV",
    T ~ "Unknown"
  )) %>%
  filter(tumor_stage != "Unknown") %>%
  select(Project = project, Barcode = barcode, TumorStage = tumor_stage) %>%
  mutate(SampleID = str_sub(Barcode, 1, 16))

annot <- clinical %>%
  select(Project = project_id, SampleID = sample) %>%
  inner_join(fpkm.annot, by = c("Project", "SampleID")) %>%
  distinct(SampleID, .keep_all = T)

## Stages I-IV to early & late
annot <- annot %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  filter(TumorStage != "Normal")

rm(clinical, fpkm.annot)

# Load mutation data and run GSEA ------------------------------
gsea.analysis <- function(dat, pws) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pws, stats = stats, eps = 0, scoreType = "pos") %>%
    arrange(padj, pval)
}

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
    # mutate(MutationNum = 1000 * MutationNum / GeneLength) %>%
    select(Project, TumorStage, Symbol, MutationNum)

  if (!dir.exists(res.dir)) {
    dir.create(res.dir, recursive = T)
  }

  print(str_glue("Working on {mutation.dir}"))

  for (proj in projects) {
    print(str_glue("---Working on {proj}"))
    res.mut <- vector("list", length = length(stages))
    names(res.mut) <- stages

    for (stage in stages) {
      set.seed(1005)
      proj.annot <- annot %>%
        filter(Project == proj & TumorStage == stage)

      proj.mutations <- mutations %>%
        filter(Project == proj & TumorStage == stage) %>%
        select(Symbol, MutationNum) %>%
        arrange(desc(MutationNum))

      res.mut[[stage]] <- gsea.analysis(proj.mutations, pathways) %>%
        filter(pval <= 0.2)
    }
    openxlsx::write.xlsx(res.mut, file = str_glue("{res.dir}/{proj}.xlsx"))
  }
}
