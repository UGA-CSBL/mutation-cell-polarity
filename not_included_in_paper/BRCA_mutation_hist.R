mutations %>%
  as_tibble() %>%
  inner_join(annot, by = "SampleID") %>%
  select(SampleID, Symbol, MutationType, TumorStage) %>%
  count(TumorStage, SampleID) %>%
  arrange(n) %>%
  gghistogram(x = "n", bins = 40, palette = "jco",
              xlab = "Number of Mutations Per Sample") %>%
  facet("TumorStage", scales = "free")

ggsave("~/num_of_mutations_per_sample_staged.png", width = 10, height = 8)

mutations %>%
  as_tibble() %>%
  inner_join(annot, by = "SampleID") %>%
  select(SampleID, Symbol, MutationType, TumorStage) %>%
  count(TumorStage, SampleID) %>%
  filter(n < 1000) %>%
  gghistogram(x = "n", bins = 40, palette = "jco",
              xlab = "Number of Mutations Per Sample") %>%
  facet("TumorStage", scales = "free")

ggsave("~/num_of_mutations_per_sample_staged_leq1000.png", width = 10, height = 8)

mutations %>%
  as_tibble() %>%
  inner_join(annot, by = "SampleID") %>%
  select(SampleID, Symbol, MutationType, TumorStage) %>%
  count(TumorStage, SampleID) %>%
  filter(n < 1000) %>%
  gghistogram(x = "n", bins = 40, palette = "jco",
              xlab = "Number of Mutations Per Sample")

ggsave("~/num_of_mutations_per_sample_leq1000.png", width = 8, height = 6)
