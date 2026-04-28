# Held-out rDDP recapture colocalization

This workflow fine-maps GEUVADIS baseline and held-out cross-modality
conditionally independent QTLs with SuSiE, then computes credible-set overlap
between original KDP modality QTLs and `latent_residual` rDDP QTLs from either
the baseline control or the matching held-out run.

Run from the repository root with:

```bash
snakemake -s held_out_coloc/Snakefile
```

Primary outputs:

- `output/qtl_coloc.all.tsv.gz`: all same-gene credible-set overlaps.
- `output/qtl_coloc.significant.tsv.gz`: overlaps passing `CLPP >= 0.01` or
  `CLPA >= 0.5`.
- `output/summary_by_modality.tsv`: per held-out modality counts of original
  credible sets colocalized in the baseline control, held-out run, and
  held-out-only subset.
