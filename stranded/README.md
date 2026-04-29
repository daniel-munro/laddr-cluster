# LaDDR Stranded Coverage Test

This analysis compares LaDDR output from the same SRA RNA-seq samples using:

- `runs/unstranded`: conventional unstranded bigWigs
- `runs/stranded`: plus/minus bigWigs as listed in `data/coverage_manifest_stranded.tsv`
- `runs/stranded_swapped`: the same plus/minus bigWigs with manifest columns swapped

The workflow starts from paired FASTQs in `data/fastq`, optionally subsets each pair with `seqtk`, aligns with STAR, builds unstranded and strand-filtered bigWigs, then runs LaDDR separately for each analysis.

Run from this directory:

```bash
snakemake -j 8 --resources mem_mb=128000
```

The subset size is controlled by `subset_read_pairs` in `config.yaml`. Set it to `0` or `null` to use full FASTQs.

The current STAR index and GTF are the small chr1 0-2 Mb test reference from the Pantry example, so this is configured as a pipeline/code-path test rather than a genome-wide analysis.
