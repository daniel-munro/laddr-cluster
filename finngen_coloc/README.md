# FinnGen Credible-Set Colocalization

This workflow performs credible-set-based colocalization between DDP/QTL SuSiE
credible sets and FinnGen R12 SuSiE fine-mapping summary files.

The analysis is downstream of `finemap/`: `finemap/` creates the QTL credible
sets, while this directory compares those credible sets against FinnGen disease
and trait credible sets.

## Inputs

Configured in `config.yaml`:

- `qtl_susie`: QTL SuSiE table, currently
  `../finemap/output/cis_susie.all_tissues.tsv.gz`.
- `finngen_traits`: explicit list of FinnGen trait IDs, one per line, currently
  `../data/finngen_finemap/traits.txt`.
- `finngen_summary_dir`: directory containing FinnGen summary files, currently
  `../data/finngen_finemap/summary`.

For every trait ID listed in `finngen_traits`, the workflow expects:

- `finngen_R12_<trait>.SUSIE_extend.snp.filter.tsv`
- `finngen_R12_<trait>.SUSIE_extend.cred.summary.tsv`

The Snakefile expands these expected paths from `traits.txt`, so missing
FinnGen files fail during normal Snakemake input checking.

## Run

From the repository root:

```bash
snakemake -s finngen_coloc/Snakefile --cores 1
```

The main script streams FinnGen data one trait at a time. The QTL table is
loaded once, and progress messages report QTL loading, QTL metadata
precomputation, and per-trait processing status.

## Variant Matching

QTL and FinnGen variants are normalized to exact `chrom:pos:ref:alt` keys.

- QTL variants are parsed from IDs like `chr1_64764_C_T_b38`.
- FinnGen variants use `chromosome`, `position`, `allele1`, and `allele2`.

Only exact allele matches are used. There is no strand flip, allele swap, or
liftover rescue in this version.

## Colocalization Metrics

Each output row is a QTL credible set x FinnGen credible set pair with at least
the configured minimum number of shared variants.

- `clpp`: eCAVIAR-style colocalization posterior probability, computed as
  `sum(qtl_pip * finngen_cs_specific_prob)` across shared variants.
- `clpa`: FinnGen-style colocalization posterior agreement, computed as
  `sum(min(qtl_pip, finngen_cs_specific_prob))` across shared variants.
- `cs_overlap`: number of shared variants between the QTL and FinnGen credible
  sets after exact variant matching.

By default, rows are restricted to FinnGen credible sets with `good_cs == true`
and at least one shared variant.

## Outputs

All outputs are written to `finngen_coloc/output/`.

### `coloc.clpp.tsv.gz`

Primary colocalization table using the default FinnGen credible-set assignment:

- FinnGen credible set column: `cs`
- FinnGen posterior column: `cs_specific_prob`

This is the main table for downstream review. It contains all QTL CS x FinnGen
CS pairs passing the configured overlap and FinnGen credible-set filters, not
only rows passing the significance thresholds.

Important columns:

- `tissue`, `phenotype_id`, `phenotype_class`: QTL phenotype metadata.
- `qtl_cs_id`, `qtl_cs_size`, `qtl_top_variant`, `qtl_top_pip`: QTL credible-set
  metadata.
- `finngen_trait`, `finngen_region`, `finngen_cs`, `finngen_cs_size`: FinnGen
  credible-set identifiers and size.
- `finngen_top_variant`, `finngen_top_prob`: top FinnGen variant within the
  matched credible set by FinnGen CS-specific probability.
- `cs_overlap`, `clpp`, `clpa`: overlap and colocalization metrics.
- `top_shared_variant`, `top_shared_qtl_pip`, `top_shared_finngen_prob`,
  `top_shared_clpp`: shared variant contributing the largest single-variant
  CLPP term.
- `finngen_cs_log10bf`, `finngen_cs_avg_r2`, `finngen_cs_min_r2`,
  `finngen_low_purity`, `finngen_good_cs`: FinnGen credible-set quality
  annotations from `*.cred.summary.tsv`.
- `finngen_lead_*`: FinnGen lead variant annotations from the credible-set
  summary file.

### `coloc.clpp_99.tsv.gz`

Sensitivity colocalization table using FinnGen 99% credible sets:

- FinnGen credible set column: `cs_99`
- FinnGen posterior column: `cs_specific_prob_99`

The schema matches `coloc.clpp.tsv.gz`. This table is intended as a sensitivity
analysis for pairs whose support changes when FinnGen 99% credible sets are
used.

### `coloc.significant.tsv.gz`

Thresholded subset of `coloc.clpp.tsv.gz`.

Rows are included when they pass either configured threshold:

- `clpp >= clpp_threshold`
- `clpa >= clpa_threshold`

The default thresholds are:

- `clpp_threshold: 0.01`
- `clpa_threshold: 0.5`

This file uses only the primary FinnGen credible sets, not the 99% sensitivity
sets.

### `summary_by_trait.tsv`

Trait-level summary of the primary colocalization table.

Columns:

- `finngen_trait`: FinnGen trait ID.
- `n_coloc_pairs`: number of QTL CS x FinnGen CS pairs in `coloc.clpp.tsv.gz`.
- `n_phenotypes`: number of distinct QTL phenotypes with at least one pair.
- `n_tissues`: number of tissues represented.
- `max_clpp`, `max_clpa`: maximum primary colocalization scores for the trait.

### `summary_by_phenotype_class.tsv`

Phenotype-class summary of the primary colocalization table.

`phenotype_class` is the prefix before the first `:` in `phenotype_id`. Phenotype
IDs without `:` are assigned to `ddp`.

Columns:

- `phenotype_class`: QTL phenotype class.
- `n_coloc_pairs`: number of QTL CS x FinnGen CS pairs in `coloc.clpp.tsv.gz`.
- `n_phenotypes`: number of distinct QTL phenotypes with at least one pair.
- `n_traits`: number of FinnGen traits represented.
- `max_clpp`, `max_clpa`: maximum primary colocalization scores for the class.

### `coloc_examples.annotated.tsv`

Selected colocalization examples from `coloc_examples.tsv`, joined back to
`coloc.clpp.tsv.gz`.

This table has one row per selected example and includes the selected
`tissue`, `phenotype_id`, `finngen_trait`, `qtl_cs_id`, `finngen_region`, and
`finngen_cs`, plus the corresponding CLPP/CLPA values, top shared variant, and
FinnGen credible-set annotations.

### `coloc_examples.qtl_region.tsv.gz`

Regional xQTL nominal association statistics for the selected examples.

The workflow subsets the GTEx residual phenotype BED to the selected examples,
replaces each temporary BED row's coordinates with that example's
`finngen_region`, runs tensorQTL `cis_nominal` for the relevant tissues, and
writes one combined table. This makes the QTL and FinnGen regional statistics
use the same plotting window. If the same QTL phenotype is selected for more
than one FinnGen region, the temporary BED uses unique internal phenotype IDs
and maps them back to the original `phenotype_id` in the output.

Important columns:

- `example_id`, `tissue`, `phenotype_id`, `qtl_cs_id`: selected example
  identifiers.
- `tensorqtl_phenotype_id`: internal phenotype ID used in the temporary BED
  and tensorQTL output.
- `variant_id`, `variant_key`: tensorQTL variant ID and normalized
  `chrom:pos:ref:alt` key.
- `pval_nominal`, `slope`, `slope_se`, `af`: xQTL nominal association
  statistics.
- `in_qtl_cs`: whether the variant is in the selected QTL credible set.
- `qtl_pip`, `qtl_af`: QTL SuSiE PIP and allele frequency for variants in the
  selected credible set.

### `coloc_examples.finngen_region.tsv.gz`

Regional FinnGen summary statistics for the selected examples, queried from
tabix-indexed files in `data/finngen/summary_stats`.

Important columns:

- `example_id`, `finngen_trait`, `finngen_region`: selected example
  identifiers.
- `#chrom`, `pos`, `ref`, `alt`, `variant_key`: FinnGen variant coordinates and
  normalized key.
- `pval`, `mlogp`, `beta`, `sebeta`, `af_alt`: FinnGen regional association
  statistics.
- `in_finngen_cs`: whether the variant is in the selected FinnGen credible set.
- `finngen_cs_specific_prob`: FinnGen CS-specific posterior probability for
  variants in the selected credible set.
- `finngen_finemap_p`, `finngen_finemap_beta`, `finngen_finemap_se`: fine-map
  SNP-file statistics for variants in the selected credible set.

## Configuration

Current defaults:

```yaml
qtl_phenotype_prefixes: []
require_finngen_good_cs: true
min_cs_overlap: 1
clpp_threshold: 0.01
clpa_threshold: 0.5
```

`qtl_phenotype_prefixes: []` means all QTL phenotypes in the SuSiE table are
included.

## Notes

- `CLPP` is sensitive to credible-set size because it sums products of posterior
  probabilities over shared variants.
- `CLPA` is included as FinnGen's size-independent agreement metric.
- Candidate examples of DDP discoveries missed by other quantification methods
  should be selected downstream after reviewing `coloc.clpp.tsv.gz` and
  `coloc.significant.tsv.gz`.
