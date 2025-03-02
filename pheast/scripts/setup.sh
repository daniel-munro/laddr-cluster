set -e

GROUP=$1
TISSUE=$2
GTF=../../ref/human/gencode.v47.annotation.gtf.gz


PROJECT_DIR="$GROUP-$TISSUE"

# Confirm PROJECT_DIR does not already exist
if [ -d "$PROJECT_DIR" ]; then
    echo "Error: $PROJECT_DIR already exists"
    exit 1
fi
cp -r ~/tools/Pantry/pheast $PROJECT_DIR
cp scripts/qtl.smk $PROJECT_DIR/steps/

# (Remove example data, edit config)
rm -r $PROJECT_DIR/input

echo "Processing $TISSUE"

mkdir -p $PROJECT_DIR/input/phenotypes_original

python3 ~/tools/Pantry/phenotyping/scripts/assemble_bed.py \
    --type latent \
    --input ../phenos/$GROUP/phenotypes/latent_phenos.$TISSUE.tsv.gz \
    --ref_anno $GTF \
    --output $PROJECT_DIR/input/phenotypes_original/latent_full.bed
bgzip $PROJECT_DIR/input/phenotypes_original/latent_full.bed

# Convert to individual IDs and filter samples
python3 ~/tools/Pantry/pheast/scripts/prepare_phenotypes.py \
    --indir $PROJECT_DIR/input/phenotypes_original \
    --out $PROJECT_DIR/input/phenotypes \
    --map ../data/gtex/sample_individual_map.tsv \
    --individuals ../../pantry/GTEx/geno/ids.txt \

mv -i $PROJECT_DIR/input/phenotypes_original/latent_full.bed.gz $PROJECT_DIR/input/phenotypes/latent_full.preprocessed.bed.gz
rmdir $PROJECT_DIR/input/phenotypes_original

tabix -p bed $PROJECT_DIR/input/phenotypes/latent_full.bed.gz

## Get latent phenotype groups
zcat $PROJECT_DIR/input/phenotypes/latent_full.bed.gz \
    | tail -n +2 \
    | cut -f4 \
    | awk '{{ g=$1; sub(/:.*$/, "", g); print $1 "\t" g }}' \
    > $PROJECT_DIR/input/phenotypes/latent_full.phenotype_groups.txt

## Get samples with genotypes
zcat $PROJECT_DIR/input/phenotypes/latent_full.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > $PROJECT_DIR/input/samples.txt

# (Run Pheast using Snakemake)
