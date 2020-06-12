library(tidyverse)
library(ggpubr)
library(foreach)
library(doSNOW)

rm(list = ls())
set.seed(1005)

# Load annotation data ------------------------------
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")

id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

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
  filter(Project %in% projects) %>%
  distinct(SampleID, .keep_all = T) %>%
  mutate(TumorStage = case_when(
    TumorStage =="Normal" ~ "Normal",
    T ~ "Tumor"
  ))

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
# See the link above for impact rating
mutation.impact <- data.table::fread("~/CSBL_shared/UCSC_Xena/DNASeq/mutation_impact.csv") %>%
  group_by(Impact) %>%
  summarise(Term = paste(Term, collapse = "|"))
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
  select(SampleID, Symbol = Symbol.y, Impact)

mut.summary <- mutations %>%
  inner_join(annot, by = "SampleID") %>%
  count(Project, TumorStage, Symbol) %>%
  inner_join(count(annot, Project, TumorStage), by = c("Project", "TumorStage")) %>%
  rename(MutationNum = n.x, SampleNum = n.y) %>%
  mutate(NormMutations = MutationNum / SampleNum) %>%
  select(Project, TumorStage, Symbol, NormMutations)

# Average expression in TPM -----
FPKM2TPM <- function(x) {
  1e6 * x / sum(x)
}

cl <- makeCluster(length(projects), outfile = "")
registerDoSNOW(cl)

dat.tpm <- foreach(proj = projects, .combine = bind_rows,
                   .packages = c("dplyr", "stringr", "tidyr")) %dopar%
  {
    print(proj)
    data.table::fread(
      str_glue("~/CSBL_shared/RNASeq/TCGA/FPKM/{proj}.FPKM.csv")
    ) %>%
      select(Ensembl, one_of(annot$Barcode[annot$Project == proj])) %>%
      mutate_if(is.numeric, FPKM2TPM) %>%
      mutate(Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
      inner_join(id.map, by = "Ensembl") %>%
      select(-Ensembl) %>%
      filter(Symbol %in% unique(mutations$Symbol)) %>%
      gather(Barcode, TPM, -Symbol) %>%
      inner_join(annot, by = "Barcode") %>%
      select(Project, TumorStage, Symbol, Barcode, TPM)

  }
stopCluster(cl)

dat.tpm.avg <- dat.tpm %>%
  group_by(Project, TumorStage, Symbol) %>%
  summarise(AvgTPM = mean(TPM, na.rm = T)) %>%
  ungroup()

# Plot number of mutations vs. average expression -----
# For each stage (early and advanced) of a cancer type, plot the relationship
# between the normalized number of mutations per gene vs. the average expession
# level of the gene as follows, where the normalized number of mutations per gene
# is the total number of mutations in a gene divided by the number of genomes:
# the x-axis is for normalized number of genes; you can bin genes toegther with
# similar numbers of mutations, and the y-axis is for the average expressions.
# For example, if you have a bin containing k genes with similar numbers of
# mutations, then plot the average expession level for each gene across all the k
# genes like a boxplot for the bin.
dat <- dat.tpm.avg %>%

  # Add in mutation data
  inner_join(mut.summary, by = c("Project", "TumorStage", "Symbol")) %>%

  # Keep expression data of normal samples
  # pivot_wider(names_from = "TumorStage", values_from = "AvgTPM") %>%
  # inner_join(dat.tpm.avg %>%
  #              filter(TumorStage == "Normal") %>%
  #              select(Project, Symbol, Normal = AvgTPM),
  #            by = c("Project", "Symbol")) %>%
  # pivot_longer(c(Tumor, Normal), names_to = "TumorStage", values_to = "AvgTPM") %>%
  # filter(!is.na(AvgTPM)) %>%
  select(Project, TumorStage, Symbol, NormMutations, AvgTPM) %>%

  # quantile-based binning
  group_by(Project) %>%
  mutate(NormMutations = cut(NormMutations,
                             breaks = unname(quantile(NormMutations, probs = c(0, 0.95, 1))),
                             include.lowest = T)) %>%
  ungroup() %>%

  # Tweak for plotting
  mutate(
    AvgTPM = log(AvgTPM + 0.1),
    # TumorStage = case_when(
    #   TumorStage == "Late" ~ "Advanced",
    #   T ~ TumorStage
    # ),
    # TumorStage = factor(TumorStage, levels = c("Normal", "Early", "Advanced")),
    TumorStage = factor(case_when(
      TumorStage == "Normal" ~ "Normal",
      T ~ "Tumor"
    )),
    Project = str_remove(Project, "^TCGA-")
  )

dat.sample.num <- count(dat, Project, TumorStage, NormMutations)

p <- ggboxplot(dat, x = "NormMutations", y = "AvgTPM",
               # group = "TumorStage", color = "TumorStage"
               ) %>%
  facet(facet.by = "Project", nrow = 3, scales = "free") %>%
  ggpar(xlab = "Normalized # of mutations",
        ylab = "Average expression in ln(TPM + 0.1)",
        legend.title = "",
        font.main = c(18, "bold", "black"),
        font.legend = 16, font.x = 16, font.y = 16,
        x.text.angle = 45,
        palette = "jco") +
  geom_text(dat = dat.sample.num, aes(x = NormMutations, y = Inf, label = n), vjust = 1)

ggsave(str_glue("~/storage/data/archive/muscle/supp_tables/Figure S1.tiff"),
       p, width = 12, height = 13, units = "in", dpi = 150)
