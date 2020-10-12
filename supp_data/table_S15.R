# Table S15: Pathway enriched by mutated gene in stage 1 and stage 4, respectively, across nine cancer types.
library(tidyverse)
library(fgsea)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- paste("Stage", c("I", "IV"), sep = " ")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

pathways <- c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)

clinical <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/GDC-PANCAN.basic_phenotype.tsv"
) %>%
  filter(str_starts(project_id, "TCGA-")) %>%
  filter(
    sample_type %in% c("Solid Tissue Normal", "Primary Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood")
  )

fpkm.annot <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/fpkm_annot.csv"
) %>%
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
  filter(TumorStage != "Normal") %>%
  filter(Project %in% projects)

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
  group_by(Impact) %>%
  summarise(Term = paste(Term, collapse = "|")) %>%
  ungroup()

mutations <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/DNASeq/GDC-PANCAN.mutect2_snv.tsv"
) %>%
  filter(Sample_ID %in% annot$SampleID) %>%
  mutate(Impact = case_when(
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "HIGH"]) ~ "HIGH",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODERATE"]) ~ "MODERATE",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODIFIER"]) ~ "MODIFIER",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "LOW"]) ~ "LOW",
    T ~ "OTHER"  # not present
  )) %>%
  filter(Impact != "LOW") %>%
  select(SampleID = Sample_ID, Symbol = gene, Impact)

transform.alias <- function(genes) {
  bind_cols(
    Alias = genes,
    Symbol = limma::alias2SymbolTable(genes)
  )
}

mutated.genes <- unique(mutations$Symbol)
mutated.genes <- transform.alias(mutated.genes) %>%
  filter(!is.na(Symbol))

mutations <- mutations %>%
  inner_join(mutated.genes, by = c("Symbol" = "Alias")) %>%
  select(SampleID, Symbol = Symbol.y, Impact) %>%
  inner_join(annot, by = "SampleID")

mutations <- mutations %>%
  count(Project, TumorStage, Symbol, Impact) %>%
  pivot_wider(names_from = Impact, values_from = n) %>%
  replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
  mutate(MutationNum = HIGH + MODERATE + MODIFIER)

# GSEA of mutations -----
gsea.analysis <- function(dat, pws) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pws, stats = stats, eps = 0, scoreType = "pos") %>%
    arrange(padj, pval)
}

res.dir <- "~/storage/data/archive/muscle/supp_tables/Table_S15"
if (!dir.exists(res.dir)) {
  dir.create(res.dir, recursive = T)
}

res <- NULL
for (project in projects) {
  print(str_glue("Working on {project}"))

  for (stage in stages) {
    set.seed(1005)
    proj.mutations <- mutations %>%
      filter(Project == project & TumorStage == stage) %>%
      arrange(desc(MutationNum), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
      select(Symbol, MutationNum)

    res.mut <- gsea.analysis(proj.mutations, pws = pathways) %>%
      mutate(Project = project, TumorStage = stage) %>%
      select(Project, TumorStage, everything())

    data.table::fwrite(res.mut, str_glue("{res.dir}/{project}_{stage}.csv"))

    res <- bind_rows(res, res.mut)
  }
}

res %>%
  as_tibble() %>%
  filter(pval <= 0.2) %>%
  unite(Project, Project, TumorStage, sep = " ") %>%
  group_by(Project) %>%
  nest() %>%
  deframe() %>%
  openxlsx::write.xlsx(str_glue("~/storage/data/archive/muscle/supp_tables/Table_S15.xlsx"))
