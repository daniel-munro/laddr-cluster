# Phenotyping Runtime Benchmark

This benchmark reruns phenotyping on the fixed 100-sample GEUVADIS subset in `../seqsim/data/samples.txt`.

Procedure:

1. Run `bash scripts/setup.sh` from `benchmark/`.
2. Edit `pantry/config.yml` and `laddr/config.yaml` if needed.
3. Run Pantry from `benchmark/pantry/`. This includes BAM generation.
4. Run LaDDR from `benchmark/laddr/`. This uses BAM files from `../pantry/intermediate/bam/` to generate bigWigs and then runs the LaDDR pipeline.
5. Summarize benchmarks from `benchmark/` with:
   `python3 scripts/summarize_benchmarks.py --samples 100 --pantry-threads 16 --bigwig-threads 16 --output-tsv benchmarks/benchmark_summary.tsv`

Reported components:

- shared `alignment`: STAR alignment and BAM indexing
- Pantry modality-group runtimes:
  - `bam shrinking`
  - `alt_TSS/alt_polyA`
  - `expression/isoforms`
  - `splicing`
  - `stability`
  - `Pantry total`
- LaDDR runtimes:
  - `bamCoverage`
  - `laddr setup`
  - `laddr binning`
  - `laddr coverage`
  - `laddr fit`
  - `laddr transform`
  - `LaDDR total`

The LaDDR benchmark does not use optional residualization against Pantry phenotype BEDs.

Outputs:

- raw benchmark files under `benchmarks/`
- summary table at `benchmarks/benchmark_summary.tsv`
