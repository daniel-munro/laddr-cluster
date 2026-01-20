# laddr-cluster

Cluster-facing workflows for the LaDDR (latent data-driven RNA phenotyping) project. This submodule contains Snakemake pipelines, setup scripts, and utilities used to train latent RNA models, run genetic analyses (xQTL/xTWAS), and package outputs for the manuscript and data release. It is designed to run on an HPC cluster with shared data directories and large memory jobs.

## How this fits into the parent repo

From the top-level `laddr` repo:

- `cluster/` (this submodule) runs heavy compute workflows and produces outputs.
- `data/`, `ref/`, and `pantry/` (outside this submodule) supply inputs and shared tooling.
- `analyses/` and `scripts/` consume outputs for exploratory analyses and for the paper.

Most pipelines here write outputs in-place under their own subdirectories, then the `cluster/repo` scripts assemble the final data repository for release.

## Core dependencies and assumptions

- `laddr` CLI is available (from the `LaDDR` package).
- Snakemake is used throughout; rules assume large-memory nodes and long runtimes.
- `Pantry` is expected at `~/tools/Pantry` in several scripts.
- External tools used by pipelines include: `tensorQTL`, `FUSION`, `STAR`, `seqtk`, `plink`, `bgzip`, `tabix`, `Rscript`, and `python3`.

## Directory map and what each pipeline does

- `phenos/`: Train latent RNA models and generate latent phenotypes from coverage.
  - Each dataset directory (`gtex-full`, `gtextcga-full`, `geuvadis-full`, etc.) contains a `Snakefile`, `config.yaml`, and `coverage_manifest.tsv`.
  - The pipeline runs `laddr binning`, `laddr coverage`, `laddr fit`, and `laddr transform`.
  - `phenos/scripts/setup.sh` shows how manifests and project directories were initialized.
- `pheast/`: QTL/TWAS workflows built on the Pantry `pheast` template.
  - `scripts/setup.sh` copies `~/tools/Pantry/pheast` into local runs and injects LaDDR-specific steps.
  - `scripts/qtl.smk` defines cis/trans tensorQTL mapping; `scripts/twas.smk` handles FUSION model building.
- `twas/`: Aggregates TWAS across traits using FUSION.
  - Consumes Pheast outputs and GWAS sumstats, then assembles per-trait and combined hit tables.
  - `scripts/FUSION.assoc_test.R` and `scripts/trait_corrs.py` handle association tests and trait correlations.
- `compare/`: Hyperparameter sweeps for binning/model settings.
  - Runs multiple `trial*` directories to compare binning strategies and PCA/FPCA model settings.
  - Produces latent phenotypes and QTL-ready BEDs for comparison.
- `held_out/`: Leave-one-modality-out ablations.
  - Generates latent phenotypes with a single Pantry modality held out, then combines with the remaining modalities for cross-modality QTL runs.
  - Modalities are: `alt_polyA`, `alt_TSS`, `expression`, `isoforms`, `splicing`, `stability`.
- `prune_anno/`: Annotation pruning experiments.
  - `scripts/prune_annotations.py` removes non-canonical isoforms in steps (e.g., 20% at a time).
  - `latent-*` and `pantry-*` directories run LaDDR/Pantry on pruned GTFs to test isoform sparsity sensitivity.
- `seqsim/`: Sequencing depth/read-length simulations.
  - Trims/subsamples FASTQs, aligns with STAR, generates bigWigs, then phenotypes and QTL inputs.
  - Used to quantify robustness to read length and depth.
- `repo/`: Data-release packaging helpers.
  - `prepare_repo.sh` copies QTL/phenotype/TWAS outputs into a publishable repository layout.
  - `repo/README.md` documents the final released data repo.
- `twas_example/`: Minimal example Snakemake workflow and example phenotype tables.

## Typical workflow (high level)

1. **Prepare coverage inputs** (bigWigs/FASTQs) under `../data/`.
2. **Run LaDDR phenotyping** in `phenos/<project>` to produce latent phenotypes.
3. **Run Pheast** in `pheast/<project>` to generate covariates, QTLs, and TWAS weights.
4. **Run TWAS aggregation** in `twas/` to compute GWAS associations and trait correlations.
5. **Assemble release repository** with `repo/prepare_repo.sh`.

Each step is a separate Snakemake project and can be run independently once inputs are in place.
