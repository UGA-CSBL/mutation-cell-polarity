# Of all the expressed pathways, how many of the mutated genes are captured?

# Section 3的30个pathways中，有% 的mutations 是出现在上调的specific pathways中。我是想看看，是否被选择突变主要出现在 上调的pathway 中？谢谢。
# 我最关心的问题是：mutations 是否是上调的一个后果，即 如果一个pathway 不上调，它与cell polarity相关的基因也不（需要）被选择有mutation。
# 举一个例子，对于Tissue development pathway, mutations 是否基本都出现在上调的specific development pathways, 下调的部分基本与mutated 基因都富集 于tissue growth and morphogenesis的部分。
# I noted that some mutations may take place in expressed but not necessarily up-regulated pathways. So you may want to classify pathways into three types, i.e., up-, down- and unchanged when analyzing mutations.

library(tidyverse)
library(fgsea)

rm(list = ls())

# Setup -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Ensembl = ensembl_gene_id, Symbol = external_gene_name)

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
  distinct(SampleID, .keep_all = T) %>%
  filter(Project %in% projects) %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  filter(TumorStage != "Normal")

rm(clinical, fpkm.annot)

# Calculate mean TPM -----
FPKMtoTPM <- function(x) {
  return(1e6 * x / sum(x))
}

exp.res <- NULL
for (project in projects) {
  df <- data.table::fread(
    str_glue("/home/yi/CSBL_shared/RNASeq/TCGA/FPKM/{project}.FPKM.csv")
  ) %>%
    mutate_if(is.numeric, FPKMtoTPM) %>%
    # Ensembl ID to gene symbol
    mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    select(one_of(c("Symbol", annot$Barcode[annot$Project == project]))) %>%
    pivot_longer(-Symbol, names_to = "Barcode", values_to = "TPM") %>%
    group_by(Symbol) %>%
    summarise(MeanTPM = mean(TPM, na.rm = T)) %>%
    ungroup() %>%
    mutate(Project = project) %>%
    select(Project, Symbol, MeanTPM)
  exp.res <- bind_rows(exp.res, df)
}

exp.res <- exp.res %>%
  as_tibble() %>%
  arrange(Project, Symbol) %>%
  filter(MeanTPM > 10)

# Pathway expression -----
gsea.analysis <- function(dat, pws) {
  stats <- dat %>%
    deframe()

  fgsea(pathways = pws,
        stats = stats,
        eps = 0,
        scoreType = "pos"
  ) %>%
    arrange(padj, pval)
}

gsea.res <- NULL
for (project in projects) {
  set.seed(1005)
  dat <- exp.res %>%
    filter(Project == project) %>%
    select(Symbol, MeanTPM) %>%
    arrange(desc(MeanTPM))

  dat <- gsea.analysis(dat, pathways) %>%
    mutate(Project = project) %>%
    relocate(Project)

  gsea.res <- bind_rows(gsea.res, dat)
}

pws <- gsea.res %>%
  as_tibble() %>%
  filter(padj <= 0.2) %>%
  distinct(Project, pathway) %>%
  nest(data = pathway) %>%
  mutate(data = map(data, ~ .x$pathway)) %>%
  deframe()

expressed.pw.genes <- map(pws, ~ unique(unlist(pathways[.x]))) %>%
  enframe() %>%
  select(Project = name, GenesInExpressedPathways = value)

# Mutations -----
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

selected.mutations <- mutations %>%
  distinct(Project, Symbol)

mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
  group_by(Impact) %>%
  summarise(Term = paste(Term, collapse = "|"))

all.mutations <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/GDC-PANCAN.mutect2_snv.tsv") %>%
  filter(Sample_ID %in% annot$SampleID) %>%
  mutate(Impact = case_when(
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "HIGH"]) ~ "HIGH",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODERATE"]) ~ "MODERATE",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODIFIER"]) ~ "MODIFIER",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "LOW"]) ~ "LOW",
    T ~ "OTHER"  # not present
  )) %>%
  filter(Impact != "LOW")

# Correct aliases -----
all.mutations <- all.mutations %>%
  select(SampleID = Sample_ID, Symbol = gene, Impact)

transform.alias <- function(genes) {
  bind_cols(
    Alias = genes,
    Symbol = limma::alias2SymbolTable(genes)
  )
}

mutated.genes <- unique(all.mutations$Symbol)
mutated.genes <- transform.alias(mutated.genes) %>%
  filter(!is.na(Symbol))

all.mutations <- all.mutations %>%
  inner_join(mutated.genes, by = c("Symbol" = "Alias")) %>%
  select(SampleID, Symbol = Symbol.y, Impact) %>%
  inner_join(annot, by = "SampleID") %>%
  distinct(Project, Symbol)

# Are there mutations not in expressed pathways -----
selected.mutations %>%
  nest(data = Symbol) %>%
  mutate(data = map(data, ~ .x$Symbol)) %>%
  rename(SelectedMutatedGenes = data) %>%
  inner_join(
    all.mutations %>%
      nest(data = Symbol) %>%
      mutate(data = map(data, ~ .x$Symbol)) %>%
      rename(AllMutatedGenes = data),
    by = "Project"
  ) %>%
  inner_join(expressed.pw.genes, by = "Project") %>%
  mutate(
    ExpressedAndSelectedMutatedGenes = map2(SelectedMutatedGenes, GenesInExpressedPathways, ~ .x[.x %in% .y]),
    ExpressedAndAllMutatedGenes = map2(AllMutatedGenes, GenesInExpressedPathways, ~ .x[.x %in% .y])
  ) %>%
  mutate_at(2:6, ~ map_int(., length)) %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  openxlsx::write.xlsx("~/mutated_and_expressed_genes.xlsx")
