library(tidyverse)
library(ggpubr)

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

# Load mutation data ------------------------------
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
# See the link above for impact rating
mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
  group_by(Impact) %>%
  summarise(Term = paste(Term, collapse = "|"))
mutations <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/GDC-PANCAN.mutect2_snv.tsv") %>%
  filter(Sample_ID %in% annot$SampleID) %>%
  mutate(Impact = case_when(
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "HIGH"]) ~ "HIGH",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODERATE"]) ~ "MODERATE",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "MODIFIER"]) ~ "MODIFIER",
    str_detect(effect, mutation.impact$Term[mutation.impact$Impact == "LOW"]) ~ "LOW",
    T ~ "OTHER"  # not present
  )) %>%
  filter(Impact != "LOW") %>%
  relocate(SampleID = Sample_ID, Symbol = gene, Impact)

# Correct aliases -----
transform.alias <- function(genes) {
  bind_cols(
    Alias = genes,
    Symbol = limma::alias2SymbolTable(genes)
  )
}

mutated.genes <- unique(mutations$gene)
mutated.genes <- transform.alias(mutated.genes) %>%
  filter(!is.na(Symbol))

mutations <- mutations %>%
  inner_join(mutated.genes, by = c("gene" = "Alias")) %>%
  select(-gene) %>%
  rename(SampleID = Sample_ID) %>%
  inner_join(annot, by = "SampleID") %>%
  relocate(Project, TumorStage, Impact, Symbol, SampleID, Barcode) %>%
  arrange(Project, TumorStage, Impact, Symbol)

mutation.freq <- mutations %>%

  # Only keep microindels
  mutate(MutationLength = end - start + 1) %>%
  filter(MutationLength < 4) %>%

  # Exclude hypermutators (more than 500 mutations)
  count(Project, TumorStage, SampleID) %>%
  # filter(n < 500) %>%
  mutate(n = n / 30) %>%  # The human exome is ~30Mb
  # group_by(Project, TumorStage) %>%
  # filter(n < quantile(n, 0.75) + 1.5*IQR(n)) %>%
  # ungroup() %>%

  group_by(Project, TumorStage) %>%
  summarise(MutationFreq = mean(n)) %>%
  ungroup()

min(mutation.freq$MutationFreq)

max(mutation.freq$MutationFreq)

mutation.freq %>%
  pivot_wider(names_from = TumorStage, values_from = MutationFreq) %>%
  mutate(Project = str_remove(Project, "^TCGA-")) %>%
  rename(`Cancer Type` = Project, Advanced = Late) %>%
  data.table::fwrite("~/storage/data/archive/muscle/mutation_frequency.csv")
