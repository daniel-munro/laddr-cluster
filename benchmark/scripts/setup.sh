#!/usr/bin/env bash

set -euo pipefail

laddr_dir=laddr

rsync -a ~/tools/Pantry/phenotyping/ pantry/ \
    --exclude input \
    --exclude intermediate \
    --exclude output \
    --exclude .snakemake

mkdir -p pantry/intermediate/reference
ln -s ../../../../seqsim/data/bam/reference/star_index_75 pantry/intermediate/reference/star_index_75

## Then edit config.yml and add benchmark arguments to Snakefile

laddr --init --template snakemake laddr
awk '{ print "benchmark\t" $1 "\t" $1 ".bw" }' ../seqsim/data/samples.txt > laddr/coverage_manifest.tsv

cat > laddr/config.yaml <<'EOF'
input:
  gtf: ../../../ref/human-ensembl/Homo_sapiens.GRCh38.113.chr.chrom.gtf
  coverage:
    method: manifest
    directory: covg_bigwig
    manifest: coverage_manifest.tsv
  min_samples_expressed: 0.5

binning:
  use_existing: false
  batch_size: 40
  max_bin_width: 1024
  adaptive:
    max_samples: 256
    bins_per_gene: 256

model:
  use_existing: false
  var_explained_max: 0.8
  n_pcs_max: 16
EOF

## Then edit configs as needed and run:
##   cd pantry && snakemake ...
##   cd ../laddr && snakemake ...
