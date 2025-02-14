library(recount3)
library(tidyverse)

proj_info <- available_projects()

samples_brain <- proj_info |>
    filter(project == "BRAIN") |>
    create_rse() |>
    colData() |>
    as_tibble()

samples_blood <- proj_info |>
    filter(project == "BLOOD") |>
    create_rse() |>
    colData() |>
    as_tibble()

set.seed(20250213)

samples <- bind_rows(
    samples_brain |>
        select(gtex.sampid, gtex.smts, BigWigURL) |>
        slice_sample(n = 8),
    samples_blood |>
        select(gtex.sampid, gtex.smts, BigWigURL) |>
        slice_sample(n = 8),    
) |>
    rename(
        sample = gtex.sampid,
        tissue = gtex.smts,
        url = BigWigURL
    )

write_tsv(samples, "samples.tsv")
