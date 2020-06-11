library(tidyverse)
library(ggpubr)

rm(list = ls())

dat <- data.table::fread("~/storage/data/archive/muscle/supp_tables/umbrella_pathways.csv")

ggscatter(dat, x = "up_cancer_types", y = "down_cancer_types") %>%
  ggpar(xlab = "# of cancer types having a pathway up-regulated",
        ylab = "# of cancer types having a pathway down-regulated",
        legend.title = "",
        font.main = c(18, "bold", "black"),
        font.legend = 16, font.x = 16, font.y = 16) +
  scale_y_continuous(breaks=0:9) +
  scale_x_continuous(breaks=0:9)

ggsave("~/storage/data/archive/muscle/supp_tables/Figure S33.tiff",
       width = 7, height = 7, units = "in", dpi = "print")
