library(recount3)
library(tidyverse)

tissues <- read_tsv(
    "data/tissueInfo.tsv",
    col_types = cols(
        tissueSiteDetail = "c",
        tissueSiteDetailAbbr = "c",
        .default = "-"
    )
) |>
    rename(tissue = tissueSiteDetailAbbr)

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
    select(sample, tissue, url) |>
    arrange(tissue, sample)

write_tsv(samples, "data/bigwigs.tsv.gz")

# Also filter to only include a smaller dataset of 5 selected tissues

gtex5 <- read_lines("data/tissues.gtex5.txt")

samples |>
    filter(tissue %in% gtex5) |>
    write_tsv("data/bigwigs.gtex5.tsv.gz")
