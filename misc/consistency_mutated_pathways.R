library(tidyverse)
library(fgsea)

rm(list = ls())

# Load annotation data ------------------------------
projects <- c(
  "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC","TCGA-KIRP",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")

sample.groups <- data.table::fread(
  "~/storage/data/archive/muscle/compare_random_groups/out_barcode.txt",
  sep = ":", sep2 = ",", header = F, col.names = c("Label", "SampleID")) %>%
  as_tibble() %>%
  mutate(Label = map(Label, ~ str_split(.x, " ")[[1]])) %>%
  mutate(Project = map_chr(Label, ~ .x[1]),
         Iteration = as.integer(map_chr(Label, ~ .x[3])),
         Group = map_chr(Label, ~ .x[5]),
         SampleID = map(SampleID, ~ str_split(.x, ",")[[1]])) %>%
  select(Project, Iteration, Group, SampleID)

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
  distinct(SampleID, .keep_all = T) %>%
  filter(Project %in% projects) %>%
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  )) %>%
  filter(TumorStage != "Normal")

rm(clinical, fpkm.annot)

pws <- openxlsx::read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                 colNames = F)$X3
pws <- pws[2:length(pws)]

# Load mutation data ------------------------------
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
# See the link above for impact rating
# if (file.exists("~/storage/data/archive/muscle/compare_random_groups/full_mutations_list.csv")) {
#   mutations <- data.table::fread("~/storage/data/archive/muscle/compare_random_groups/full_mutations_list.csv")
# } else {
#   mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
#     group_by(Impact) %>%
#     summarise(Term = paste(Term, collapse = "|"))
#   mutations <- data.table::fread(
#     "~/CSBL_shared/UCSC_Xena/DNASeq/GDC-PANCAN.mutect2_snv.tsv"
#   ) %>%
#     filter(Sample_ID %in% annot$SampleID) %>%
#     mutate(Impact = case_when(
#       str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "HIGH"]) ~ "HIGH",
#       str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODERATE"]) ~ "MODERATE",
#       str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODIFIER"]) ~ "MODIFIER",
#       str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "LOW"]) ~ "LOW",
#       T ~ "OTHER"  # not present
#     )) %>%
#     filter(Impact != "LOW") %>%
#     select(SampleID = Sample_ID, Symbol = gene, Impact)
#
#   transform.alias <- function(genes) {
#     bind_cols(
#       Alias = genes,
#       Symbol = limma::alias2SymbolTable(genes)
#     )
#   }
#
#   mutated.genes <- unique(mutations$Symbol)
#   mutated.genes <- transform.alias(mutated.genes) %>%
#     filter(!is.na(Symbol))
#
#   mutations <- mutations %>%
#     inner_join(mutated.genes, by = c("Symbol" = "Alias")) %>%
#     select(SampleID, Symbol = Symbol.y, Impact) %>%
#     inner_join(annot, by = "SampleID")
#
#   data.table::fwrite(mutations, "~/storage/data/archive/muscle/compare_random_groups/full_mutations_list.csv")
# }
mutations <- data.table::fread("~/storage/data/archive/muscle/microindel_list.csv")

# GSEA of mutations -----
gsea.analysis <- function(dat, pws) {
  stats <- dat %>%
    deframe()

  fgseaMultilevel(pathways = pws, stats = stats,
                   eps = 0, scoreType = "pos") %>%
    arrange(padj, pval)
}

for (project in projects) {
  res <- NULL
  proj.mutations <- mutations %>%
    filter(Project == project)

  print(str_glue("Working on {project}...\n"))
  pb <- txtProgressBar(0, 99, style = 3)
  for (i in 0:99) {
    set.seed(1005)
    # Group A
    sample.ids <- sample.groups %>%
      filter(Project == project & Iteration == i & Group == "A") %>%
      pull(SampleID) %>%
      unlist()
    dat <- proj.mutations %>%
      filter(SampleID %in% sample.ids) %>%
      count(Symbol, Impact) %>%
      spread(Impact, n) %>%
      replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
      mutate(MutationNum = HIGH + MODERATE + MODIFIER) %>%
      arrange(desc(MutationNum), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
      select(Symbol, MutationNum) %>%
      filter(MutationNum > 1) %>%
      gsea.analysis(pws = pathways) %>%
      mutate(Project = project, Iteration = i, Group = "A") %>%
      filter(pval <= 0.2) %>%
      select(Project, Iteration, Group, everything())
    res <- bind_rows(res, dat)

    # Group B
    sample.ids <- sample.groups %>%
      filter(Project == project & Iteration == i & Group == "B") %>%
      pull(SampleID) %>%
      unlist()
    dat <- proj.mutations %>%
      filter(SampleID %in% sample.ids) %>%
      count(Symbol, Impact) %>%
      spread(Impact, n) %>%
      replace_na(list(HIGH = 0, MODERATE = 0, MODIFIER = 0)) %>%
      mutate(MutationNum = HIGH + MODERATE + MODIFIER) %>%
      arrange(desc(MutationNum), desc(HIGH), desc(MODERATE), desc(MODIFIER), Symbol) %>%
      select(Symbol, MutationNum) %>%
      filter(MutationNum > 1) %>%
      gsea.analysis(pws = pathways) %>%
      mutate(Project = project, Iteration = i, Group = "B") %>%
      filter(pval <= 0.2) %>%
      select(Project, Iteration, Group, everything())
    res <- bind_rows(res, dat)
    setTxtProgressBar(pb, i)
  }
  data.table::fwrite(res, file = str_glue(
    "~/storage/data/archive/muscle/compare_random_groups/{project}.csv"))
}

# Calculate avg concordant pathways -----
res <- NULL
for (project in projects) {
  dat <- data.table::fread(str_glue(
    "~/storage/data/archive/muscle/compare_random_groups/{project}.csv"))
  dat <- dat %>%
    filter(pval <= 0.05) %>%
    select(Iteration, Group, pathway) %>%
    filter(pathway %in% pws) %>%
    group_by(Iteration, Group) %>%
    nest(data = c(pathway)) %>%
    ungroup() %>%
    mutate(data = map(data, ~ .x$pathway)) %>%
    pivot_wider(names_from = Group, values_from = data) %>%
    # mutate(IntersectPathways = map2(A, B, ~ intersect(.x, .y)),
    #        UnionPathways = map2(A, B, ~ union(.x, .y))) %>%
    # mutate(
    #   Project = project,
    #   ConsistentRate = map2_dbl(IntersectPathways, UnionPathways,
    #                              ~ length(.x) / length(.y))) %>%
    mutate(IntersectPathways = map2(A, B, ~ intersect(.x, .y))) %>%
    mutate(
      Project = project,
      ConsistentRate = map_dbl(IntersectPathways, ~ length(.x) / length(pws))) %>%
    select(Project, Iteration, ConsistentRate)

  res <- bind_rows(res, dat)
}

data.table::fwrite(res, file = "~/storage/data/archive/muscle/compare_random_groups/pathway_consistency.csv")
