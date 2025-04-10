set -e

python scripts/prune_annotations.py \
    --gtf ../../ref/human/gencode.v47.annotation.gtf.gz \
    --output_prefix data/gencode.v47.annotation \
    --remove_percent 20 \
    --steps 3

for percent in 40 60 80; do
    gzip data/gencode.v47.annotation.pruned_${percent}.gtf
    latent-rna init latent-${percent} --template snakemake
    grep '^BRNCTXB' ../phenos/gtex5-full/coverage_manifest.tsv > latent-${percent}/coverage_manifest.tsv
done

## Edit each config file to use the new gtf file

## Run `latent-rna setup` in each directory

