# Phenotyping Runtime Benchmark

This benchmark reruns phenotyping on the fixed 100-sample GEUVADIS subset in `../seqsim/data/samples.txt` and reports:

- shared `alignment`: STAR alignment, BAM shrinking, and BAM indexing
- shared `bamCoverage`: bigWig generation from the shared BAMs
- Pantry modality-group runtimes:
  - `alt_TSS/alt_polyA`
  - `expression/isoforms`
  - `splicing`
  - `stability`
  - `Pantry total`
- LaDDR runtimes:
  - `laddr setup`
  - `laddr binning`
  - `laddr coverage`
  - `laddr fit`
  - `laddr transform`
  - `LaDDR total`

The LaDDR benchmark does not use optional residualization against Pantry phenotype BEDs.

Excluded from totals:

- STAR genome index generation
- Pantry template copy/setup

The alignment step reuses the existing STAR index at `../seqsim/data/bam/reference/star_index_75/`.

Outputs:

- raw benchmark files under `benchmarks/`
- summary tables at `results/benchmark_summary.tsv` and `results/benchmark_summary.csv`
