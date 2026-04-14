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

## Then edit config.yaml and add benchmark arguments to Snakefile
## Also add bigWig generation to Snakefile
