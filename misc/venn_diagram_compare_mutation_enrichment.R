library(tidyverse)
library(openxlsx)
library(VennDiagram)

rm(list = ls())

# Mutation GSEA results -----
projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
              "TCGA-LIHC", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
stages <- c("Early", "Late")

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

if (!file.exists("~/pval_mutation_GSEA.csv")) {
  res <- NULL
  for (res.dir in res.dirs) {
    for (project in projects) {
      for (stage in stages) {
        dat <- read.xlsx(str_glue("{res.dir}/{project}.xlsx"), sheet = stage) %>%
          mutate(
            Project = project,
            Stage = stage,
            Cutoff = str_remove(res.dir, str_glue(
              "^{data.dir}/GSEA_early_late/mutation_by_"))
          )
        res <- bind_rows(res, dat)
      }
    }
  }
  data.table::fwrite(res, file = "~/pval_mutation_GSEA.csv")
} else {
  res <- data.table::fread("~/pval_mutation_GSEA.csv")
}

res.ori <- NULL
ori.data.dir <- "~/storage/data/archive/muscle/mutation_freq/GSEA_early_late"
for (project in projects) {
  for (stage in stages) {
    dat <- read.xlsx(str_glue("{ori.data.dir}/{project}.xlsx"), sheet = stage) %>%
      mutate(Project = project, Stage = stage, Cutoff = "Original")
    res.ori <- bind_rows(res.ori, dat)
  }
}

res <- as_tibble(res)
res.ori <- as_tibble(res.ori)

dat <- res %>%
  bind_rows(res.ori) %>%
  # filter(padj <= 0.2) %>%
  group_by(Project, Stage, Cutoff) %>%
  summarise(Pathway = paste(pathway, collapse = "|")) %>%
  ungroup() %>%
  mutate(Pathway = map(Pathway, ~ str_split(.x, "\\|")[[1]]))

# 79 pathways -----
pws <- read.xlsx("~/storage/data/archive/muscle/supp_tables/Table S6.xlsx",
                 colNames = F)$X3
pws <- pws[2:length(pws)]

# Venn diagrams -----
myCol <- ggsci::pal_jco()(3)
for (project in projects) {
  for (stage in stages) {
    if (!dir.exists(str_glue("{data.dir}/venn/{project}/{stage}"))) {
      dir.create(str_glue("{data.dir}/venn/{project}/{stage}"), recursive = T)
    }

    x.all <- dat %>%
      filter(Project == project & Stage == stage) %>%
      select(Cutoff, Pathway)

    x.ori <- x.all$Pathway[x.all$Cutoff == "Original"][[1]]
    x.all <- x.all %>%
      filter(Cutoff != "Original")

    if (nrow(x.all) == 0) {
      next
    }

    for (i in seq(nrow(x.all))) {
      category.name <- str_replace_all(x.all$Cutoff[i], "/", "_")
      venn.diagram(
        x = list(x.ori, x.all$Pathway[[i]], pws),
        category.names = c("Original", category.name, "78 Pathways"),
        filename = str_glue("{data.dir}/venn/{project}/{stage}/{category.name}.png"),

        # Output features
        imagetype="png" ,
        resolution = 300,

        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,

        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",

        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.fontfamily = "sans",
      )

    }
  }
}

# Clean up logs
for (f in (list.files(str_glue("{data.dir}/venn/"),
                      pattern = "\\.log$",
                      full.names = T,
                      recursive = T))) {
  file.remove(f)
}
