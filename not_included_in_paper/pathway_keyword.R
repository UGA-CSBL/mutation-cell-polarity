library(tidyverse)

rm(list = ls())

gsea.dir <- "~/storage/data/archive/muscle/mutation_freq/GSEA_stages"
projects <- "TCGA-KIRP"
stages <- paste0("Stage ", c("I", "II", "III", "IV"))

keywords <- data.table::fread("~/storage/data/archive/muscle/polarity_dependent_bio_functions.csv")

for (proj in projects) {
  res <- vector("list", length(stages))
  names(res) <- stages
  for (stage in stages) {
    dat <- openxlsx::read.xlsx(str_glue("{gsea.dir}/{proj}.xlsx"), sheet = stage) %>%
      as_tibble() %>%
      filter(pval <= 0.2)

    pws <- dat$pathway
    dat.kw <- map(keywords$Regex, ~ pws[str_detect(pws, .x)])
    names(dat.kw) <- keywords$Keyword
    dat.kw <- dat.kw %>%
      enframe() %>%
      rename(BioFunction = name, Pathway = value) %>%
      unnest(Pathway) %>%
      group_by(Pathway) %>%
      summarise(BioFunction = paste(BioFunction, collapse = ", ")) %>%
      ungroup()

    dat.res <- dat %>%
      left_join(dat.kw, by = c("pathway" = "Pathway")) %>%
      select(Pathway = pathway, BioFunction, everything()) %>%
      arrange(pval, BioFunction)
    res[[stage]] <- dat.res
  }

  openxlsx::write.xlsx(res, file = str_glue(
    "~/storage/data/archive/muscle/polarity_dependent_bio_func/{proj}.xlsx"
  ))
}
