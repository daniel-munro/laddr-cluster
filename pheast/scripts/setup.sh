set -e

group=$1
tissue=$2
gtf=../../ref/human/gencode.v47.annotation.gtf.gz


projdir="$group-$tissue"

# Confirm projdir does not already exist
if [ -d "$projdir" ]; then
    echo "Error: $projdir already exists"
    exit 1
fi
cp -r ~/tools/Pantry/pheast $projdir
cp scripts/qtl.smk $projdir/steps/

# (Remove example data, edit config)
rm -r $projdir/input

echo "Processing $tissue"

mkdir -p $projdir/input/phenotypes_original

python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
    --type latent \
    --input ../phenos/$group/phenotypes/latent_phenos.$tissue.tsv.gz \
    --ref_anno $gtf \
    --output $projdir/input/phenotypes_original/latent_full.unnorm.bed

python3 ~/tools/Pantry/phenotyping/scripts/normalize_phenotypes.py \
    --input $projdir/input/phenotypes_original/latent_full.unnorm.bed \
    --output $projdir/input/phenotypes_original/latent_full.bed
bgzip $projdir/input/phenotypes_original/latent_full.bed

# Convert to individual IDs and filter samples
python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
    --indir $projdir/input/phenotypes_original \
    --out $projdir/input/phenotypes \
    --map ../data/gtex/sample_individual_map.tsv \
    --individuals ../../pantry/GTEx/geno/ids.txt

tabix -p bed $projdir/input/phenotypes/latent_full.bed.gz

## Get latent phenotype groups
zcat $projdir/input/phenotypes/latent_full.bed.gz \
    | tail -n +2 \
    | cut -f4 \
    | awk '{{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }}' \
    > $projdir/input/phenotypes/latent_full.phenotype_groups.txt

## Get samples with genotypes
zcat $projdir/input/phenotypes/latent_full.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > $projdir/input/samples.txt

# (Run Pheast using Snakemake)
