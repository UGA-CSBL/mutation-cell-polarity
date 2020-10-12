library(tidyverse)
library(fgsea)
library(doSNOW)

rm(list = ls())
set.seed(1005)

# Load annotation data ------------------------------
id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  select(Symbol = external_gene_name, Ensembl = ensembl_gene_id)

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

clinical <- data.table::fread(
  "~/CSBL_shared/UCSC_Xena/GDC-PANCAN.basic_phenotype.tsv"
) %>%
  filter(project_id %in% projects) %>%
  filter(sample_type %in% c("Solid Tissue Normal",
                            "Primary Tumor",
                            "Primary Blood Derived Cancer - Peripheral Blood"))

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
  mutate(TumorStage = case_when(
    TumorStage %in% c("Stage I", "Stage II") ~ "Early",
    TumorStage %in% c("Stage III", "Stage IV") ~ "Late",
    T ~ "Normal"
  ))

rm(clinical, fpkm.annot)

# Load mutation data ------------------------------
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

mutations <- mutations %>%
  select(Project, TumorStage, Symbol, n = MutationNum) %>%
  group_by(Project, Symbol) %>%
  summarise(n = sum(n)) %>%
  ungroup()

if (!dir.exists(str_glue("{data.dir}/GSEA_all_samples/GSEA_summary"))) {
  dir.create(str_glue("{data.dir}/GSEA_all_samples/GSEA_summary", recursive = T))
}

if (!file.exists(str_glue("{data.dir}/GSEA_all_samples/GSEA_summary/mutations.xlsx"))) {
  mutation.summary <- mutations %>%
    group_by(Project) %>%
    summarise(nGenes = n(),
              nMutations = sum(n)) %>%
    ungroup() %>%
    inner_join(count(annot, Project),
               by = c("Project")) %>%
    rename(SampleNum = n)
  openxlsx::write.xlsx(list(Mutations = mutations, Summary = mutation.summary),
                       file = str_glue("{data.dir}/GSEA_all_samples/GSEA_summary/mutations.xlsx"))
}

# DEA results -----
dea.res <- NULL
for (proj in projects) {
  dat <- data.table::fread(str_glue(
    "~/storage/data/TCGA/DEA/TumorVsNormal/{proj}.csv")) %>%
    rename(Ensembl = V1) %>%
    mutate(Project = proj,
           Ensembl = str_remove(Ensembl, "\\.\\d+$")) %>%
    inner_join(id.map, by = "Ensembl") %>%
    select(-Ensembl) %>%
    select(Project, Symbol, log2FC = log2FoldChange, p_adjusted = padj)
  dea.res <- bind_rows(dea.res, dat)
}

# Frequency of mutations in interested pathways -----
## Note that genes with only 1 mutation was filtered out during GSEA
## So the `size` doesn't match the actual `nPathwayGenes`
pathways <- enframe(c(
  gmtPathways("~/storage/data/annotation/MSigDB/c5.bp.v7.1.symbols.gmt"),
  gmtPathways("~/storage/data/annotation/MSigDB/c2.cp.reactome.v7.1.symbols.gmt")
)) %>%
  rename(Pathway = name, Symbols = value) %>%
  mutate(nPathwayGenes = map_int(Symbols, length))

cl <- makeCluster(length(projects) + 1, outfile = "")
registerDoSNOW(cl)
foreach(project = projects,
        .packages = c("tibble", "tidyr", "purrr", "dplyr", "stringr", "openxlsx")) %dopar%
  {
    data.dir <- "~/storage/data/archive/muscle/selected_mutations"
    proj.res <- vector("list", length(stages))
    names(proj.res) <- stages
    print(str_glue("Working on project {project}"))
    gsea.res <- data.table::fread(str_glue(
      "{data.dir}/GSEA_all_samples/filtered_mutations_all_samples_GSEA.csv"
    )) %>%
      as_tibble() %>%
      filter(Project == project & pval <= 0.2)

    if (nrow(gsea.res) == 0) {
      next()
    }
    proj.mutations <- mutations %>%
      filter(Project == project) %>%
      select(Symbol, n)

    proj.dea <- dea.res %>%
      filter(Project == project) %>%
      select(Symbol, log2FC, p_adjusted)

    proj.res <- gsea.res %>%
      inner_join(pathways, by = c("pathway" = "Pathway")) %>%
      mutate(leadingEdge = str_split(leadingEdge, ", "),
             nLeadingEdge = map_int(leadingEdge, length),
             PrMutatedGenes = size / nPathwayGenes,
             PrLeadingEdge = nLeadingEdge / nPathwayGenes) %>%
      mutate(
        nPathwayMutations = map_int(Symbols, ~ proj.mutations %>%
                                      filter(Symbol %in% .x) %>%
                                      pull(n) %>% sum()),
        nLeadingEdgeMutations = map_int(leadingEdge, ~ proj.mutations %>%
                                          filter(Symbol %in% .x) %>%
                                          pull(n) %>% sum())
      ) %>%
      select(Pathway = pathway, `p-value` = pval, `Adjusted p-value` = padj,
             nPathwayGenes, nMutatedGenes = size, PrMutatedGenes,
             nLeadingEdge, PrLeadingEdge,
             nPathwayMutations, nLeadingEdgeMutations,
             LeadingEdge = leadingEdge, PathwayGenes = Symbols
      )

    data.table::fwrite(
      proj.res,
      file = str_glue("{data.dir}/GSEA_all_samples/GSEA_summary/{project}.csv")
    )
  }
stopCluster(cl)
