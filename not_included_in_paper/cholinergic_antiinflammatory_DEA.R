library(tidyverse)
library(data.table)

rm(list = ls())
options(datatable.print.nrows = 10)
source("~/storage/R_script/snippets/get_TCGA_DE.R")

projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-LIHC",
              "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
raw.gene.list <- fread("~/storage/data/archive/muscle/CAP_genelist.csv",
                       header = F)$V1
raw.gene.list <- unique(raw.gene.list)
raw.gene.list[raw.gene.list == "AP2A"] <- "TFAP2A"
raw.gene.list[raw.gene.list == "ACHA7"] <- "CHRNA7"
raw.gene.list[raw.gene.list == "SC6A2"] <- "SLC6A2"

# Get DE genes for cholinergic anti-inflammatory pathway (CAP) genes
res.cap <- get.DE(raw.gene.list, projects = projects, stages = T)

# Get TFs in different tissue types
tf.path <- "~/CSBL_shared/transcription_factor/FANTOM5_individual_networks/394_individual_networks"
tissues <- c(
  "breast_adult",
  "breast_carcinoma_cell_line",
  "colon_adult",
  "colon_carcinoma_cell_line",
  "kidney_adult",
  "embryonic_kidney_cell_line",
  "lung_adult",
  "lung_adenocarcinoma_cell_line",
  "large_cell_lung_carcinoma_cell_line",
  "stomach_fetal",
  "thyroid_adult",
  "thyroid_carcinoma_cell_line"
)

tf.res <- NULL
for (tissue in tissues) {
  dat <- fread(str_glue("{tf.path}/{tissue}.txt.gz"))
  dat <- dat[, .(TF = V1, Symbol = V2, edge_weight = V3)][
    Symbol %in% res.cap$Symbol]
  dat$tissue_type <- tissue
  tf.res <- bind_rows(tf.res, dat)
}

tf.res[edge_weight >= 0.3][, .(.N), tissue_type]

res.tf <- get.DE(unique(tf.res[edge_weight >= 0.3, Symbol]), projects = projects)

