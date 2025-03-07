library(recount3)
library(tidyverse)

##########
## GTEx ##
##########

tissues <- read_tsv(
  "../data/gtex/tissueInfo.tsv",
  col_types = cols(
    tissueSiteDetail = "c",
    tissueSiteDetailAbbr = "c",
    .default = "-"
  )
) |>
  rename(dataset = tissueSiteDetailAbbr)

proj_info <- available_projects() |>
  filter(file_source == "gtex")

samples <- proj_info |>
  reframe(
    tibble(project, organism, file_source, project_home, project_type, n_samples) |>
      create_rse() |>
      colData() |>
      as_tibble() |>
      select(
        sample = gtex.sampid,
        tissueSiteDetail = gtex.smtsd,
        url = BigWigURL
      ),
    .by = project
  ) |>
  as_tibble() |>
  left_join(tissues, by = "tissueSiteDetail", relationship = "many-to-one") |>
  select(sample, dataset, url) |>
  arrange(dataset, sample)

write_tsv(samples, "../data/gtex/bigwigs.tsv.gz")

##########
## TCGA ##
##########

proj_info <- available_projects() |>
  filter(file_source == "tcga")

samples <- proj_info |>
  reframe(
    tibble(project, organism, file_source, project_home, project_type, n_samples) |>
      create_rse() |>
      colData() |>
      as_tibble() |>
      select(
        sample = tcga.tcga_barcode,
        dataset = study,
        url = BigWigURL
      ),
    .by = project
  ) |>
  as_tibble() |>
  select(sample, dataset, url) |>
  arrange(dataset, sample)

write_tsv(samples, "../data/tcga/bigwigs.tsv.gz")
